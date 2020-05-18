#!/usr/bin/env Rscript
# @Author: kt16
# @Date:   2019-12-04 16:14:42
# @Last Modified by:   kt16
# @Last Modified time: 2020-02-16 10:11:25
# @Notes: simple R script to call clones using bcRep's structure using changeo object

options(warn=-1)
if(!require(optparse)){
  install.packages("optparse")
  suppressMessages(library(optparse))
}

option_list = list(
  make_option(c("-c", "--changeo"), type = "character", default = NULL,
    help = "name of original changeo dataframe", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
    help = "name of output file for new changeo dataframe", metavar = "character"),
  make_option(c("-i", "--identity"), type = "numeric", default = 0.85,
    help = "similiarity threshold", metavar = "character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$changeo)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call. = TRUE)}
  if (is.null(opt$output)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (output 2 (processed) file).n", call. = TRUE)}
    if (is.null(opt$identity)){
      print_help(opt_parser)
      stop("At least one argument must be supplied (output 3 (new changeo) file).n", call. = TRUE)}


### program start
      require(parallel)      
      if(!require(doParallel)){
        install.packages("doParallel")
        library(doParallel)
      }     
      if(!require(dplyr)){
        install.packages("dplyr")
        library(dplyr)
      }
      if(!require(seqinr)){
        install.packages("seqinr")
        library(seqinr)
      }

      findClones <- function(changeo, identity = opt$identity, useJ = TRUE, dispD = TRUE,
        dispSeqID = TRUE, dispCDR3aa = TRUE, dispCDR3nt = TRUE, dispJunctionFr.ratio = TRUE,
        dispJunctionFr.list = TRUE, dispFunctionality.ratio = TRUE, dispFunctionality.list = TRUE,
        dispTotalSeq = TRUE, nrCores = detectCores()){

    # only cluster using IGH
        changeo <- changeo %>% filter(LOCUS == 'IGH' & FUNCTIONAL == 'TRUE' & IN_FRAME == 'TRUE' & STOP == 'FALSE')

        changeo$JUNCTION_AA <- NA
        for (x in 1:nrow(changeo)){
          changeo$JUNCTION_AA[x] = paste0(translate(s2c(changeo$JUNCTION[x])), collapse ='')
        }    

        cl <- makeCluster(nrCores)
        registerDoParallel(cl)

        V <- unlist(unique(apply(data.frame(changeo$V_CALL), 1, function(x){strsplit(x, split = " |,|;|[*]")[[1]]})))
        V <- unique(V[grep("V", V)])
        J <- unlist(unique(apply(data.frame(changeo$J_CALL), 1, function(x){strsplit(x, split = " |,|;|[*]")[[1]]})))
        J <- unique(J[grep("J", J)])

        tempout<-vector()
        i<-NULL
        clonelist<-foreach(i=1:length(V)) %dopar% {
          if(useJ==T){
            for(j in J){
              changeo.sub <- changeo[intersect(grep(paste(V[i],"[!/*]",sep=""),changeo$V_CALL,perl=T), grep(paste(j,"[!/*]",sep=""),changeo$J_CALL,perl=T)), c('JUNCTION_AA', 'JUNCTION', 'V_CALL','J_CALL','D_CALL','FUNCTIONAL', 'SEQUENCE_ID', 'IN_FRAME', 'SEQUENCE_INPUT')]

              CDR3length <- apply(data.frame(changeo.sub$JUNCTION_AA),1,function(x){nchar(x)})
              uniqueCDR3length <- sort(as.numeric(unique(CDR3length)))
              if(length(uniqueCDR3length) == 1){ ### only 1 CDR3 length
              temp<-vector()
              uniqueCDR3.sub<-unique(changeo.sub$JUNCTION_AA)
              a.dist<-adist(unique(uniqueCDR3.sub))
              a.dist[upper.tri(a.dist,diag = T)]<-NA
              tr<-floor(uniqueCDR3length*(1-identity))
              for(a in nrow(a.dist):1){
               if(length(which(a.dist[a,]<=tr))>0){
                if(a==nrow(a.dist) || sum(apply(data.frame(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),1,function(x){length(grep(gsub("[*]","-",x),gsub("[*]","-",temp),perl=T))}))<length(intersect(c(which(a.dist[a,]<=tr),a),which(!is.na(a.dist[a,]))))){
                  temp<-c(temp,do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")))

                  changeo.new<-changeo.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",changeo.sub$JUNCTION_AA),perl=T),]

                  tempout <-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")), # shared CDR3 seq.
                   uniqueCDR3length, # CDR3 length
                   length(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), # number_shared_CDR3
                   nrow(changeo.new), # number all sequences, belonging to clone
                   do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",changeo.new$JUNCTION_AA),perl=T))})), sep=", ")), # sequence count
                   V[i], # V_gene
                   do.call(paste, c(as.list(unique(changeo.new$V_CALL)), sep=", ")), # V_gene & allele
                   if(useJ==T){j}else{do.call(paste, c(as.list(unique(changeo.new$J_CALL)), sep=", "))}, # J_gene
                   do.call(paste, c(as.list(unique(changeo.new$J_CALL)), sep=", ")), # J_gene & allele

                   if(dispD==T){do.call(paste, c(as.list(unique(changeo.new$D_CALL)), sep=", "))},
                   if(dispCDR3aa==T){do.call(paste, c(as.list(changeo.new$JUNCTION_AA), sep=", "))},
                   if(dispCDR3nt==T){do.call(paste, c(as.list(changeo.new$JUNCTION), sep=", "))},
                   if(dispFunctionality.list==T){do.call(paste, c(as.list(changeo.new$FUNCTIONAL), sep=", "))},
                   if(dispFunctionality.ratio==T){length(grep("^T",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
                   if(dispFunctionality.ratio==T){length(grep("^F",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
                   if(dispFunctionality.ratio==T){length(grep("^T|^F", changeo.new$FUNCTIONAL,invert=T))/nrow(changeo.new)},
                   if(dispJunctionFr.list==T){do.call(paste, c(as.list(changeo.new$IN_FRAME), sep=", "))},
                   if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="T"))/nrow(changeo.new)},
                   if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="F"))/nrow(changeo.new)},
                   if(dispJunctionFr.ratio==T){length(grep("^T|^F",changeo.new$IN_FRAME,invert=T))/nrow(changeo.new)},
                   if(dispSeqID==T){do.call(paste, c(as.list(changeo.new$SEQUENCE_ID), sep=", "))},
                   if(dispTotalSeq==T){do.call(paste, c(as.list(unique(changeo.new$SEQUENCE_INPUT)), sep=", "))}))
                }
              }else if(length(which(a.dist[,a]<=tr))==0 && length(grep(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep=""), changeo.sub$JUNCTION_AA))>1){
                temp<-c(temp,uniqueCDR3.sub[a])

                changeo.new<-changeo.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",changeo.sub$JUNCTION_AA),perl=T),]

                tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[a]), sep=", ")), # shared CDR3 seq.
                 uniqueCDR3length, # CDR3 length
                 length(uniqueCDR3.sub[a]), # number_shared_CDR3
                 nrow(changeo.new), # number all sequences, belonging to clone
                 do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",changeo.new$JUNCTION_AA),perl=T))})), sep=", ")), # sequence count
                 V[i], # V_gene
                 do.call(paste, c(as.list(unique(changeo.new$V_CALL)), sep=", ")), # V_gene & allele
                 if(useJ==T && length(j)>0){do.call(paste, c(as.list(unique(j)), sep=", "))}else if(useJ==T && length(j)==0){"no J"}, # J_gene
                 do.call(paste, c(as.list(unique(changeo.new$J_CALL)), sep=", ")), # J_gene & allele
                 if(dispD==T){do.call(paste, c(as.list(unique(changeo.new$D_CALL)), sep=", "))},
                 if(dispCDR3aa==T){do.call(paste, c(as.list(changeo.new$JUNCTION_AA), sep=", "))},
                 if(dispCDR3nt==T){do.call(paste, c(as.list(changeo.new$JUNCTION), sep=", "))},
                 if(dispFunctionality.list==T){do.call(paste, c(as.list(changeo.new$FUNCTIONAL), sep=", "))},
                 if(dispFunctionality.ratio==T){length(grep("^T",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
                 if(dispFunctionality.ratio==T){length(grep("^F",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
                 if(dispFunctionality.ratio==T){length(grep("^T|^F",changeo.new$FUNCTIONAL,invert=T))/nrow(changeo.new)},
                 if(dispJunctionFr.list==T){do.call(paste, c(as.list(changeo.new$IN_FRAME), sep=", "))},
                 if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="T"))/nrow(changeo.new)},
                 if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="F"))/nrow(changeo.new)},
                 if(dispJunctionFr.ratio==T){length(grep("^T|^F",changeo.new$IN_FRAME,invert=T))/nrow(changeo.new)},
                 if(dispSeqID==T){do.call(paste, c(as.list(changeo.new$SEQUENCE_ID), sep=", "))},
                 if(dispTotalSeq==T){do.call(paste, c(as.list(unique(changeo.new$SEQUENCE_INPUT)), sep=", "))}))
              }
            }
          }else{ ### several CDR3 length
          for(l in uniqueCDR3length){
            temp<-vector()
            uniqueCDR3.sub<-unique(changeo.sub$JUNCTION_AA[which(CDR3length==l)])
            a.dist<-adist(unique(uniqueCDR3.sub))
            a.dist[upper.tri(a.dist,diag = T)]<-NA
            tr<-floor(l*(1-identity))
            for(a in nrow(a.dist):1){
              if(length(which(a.dist[a,]<=tr))>0){
                if(a==nrow(a.dist) || sum(apply(data.frame(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),1,function(x){length(grep(gsub("[*]","-",x),gsub("[*]","-",temp),perl=T))}))<length(intersect(c(which(a.dist[a,]<=tr),a),which(!is.na(a.dist[a,]))))){
                  temp<-c(temp,do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")))

                  changeo.new<-changeo.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",changeo.sub$JUNCTION_AA),perl=T),]
                  
                  tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")), # shared CDR3 seq.
                   l, # CDR3 length
                   length(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), # number_shared_CDR3
                   nrow(changeo.new), # number all sequences, belonging to clone
                   do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",changeo.new$JUNCTION_AA),perl=T))})), sep=", ")), # sequence count
                   V[i], # V_gene
                   do.call(paste, c(as.list(unique(changeo.new$V_CALL)), sep=", ")), # V_gene & allele
                   if(useJ==T && length(j)>0){do.call(paste, c(as.list(unique(j)), sep=", "))}else if(useJ==T && length(j)==0){"no J"}, # J_gene
                   do.call(paste, c(as.list(unique(changeo.new$J_CALL)), sep=", ")), # J_gene & allele

                   if(dispD==T){do.call(paste, c(as.list(unique(changeo.new$D_CALL)), sep=", "))},
                   if(dispCDR3aa==T){do.call(paste, c(as.list(changeo.new$JUNCTION_AA), sep=", "))},
                   if(dispCDR3nt==T){do.call(paste, c(as.list(changeo.new$JUNCTION), sep=", "))},
                   if(dispFunctionality.list==T){do.call(paste, c(as.list(changeo.new$FUNCTIONAL), sep=", "))},
                   if(dispFunctionality.ratio==T){length(grep("^T",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
                   if(dispFunctionality.ratio==T){length(grep("^F",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
                   if(dispFunctionality.ratio==T){length(grep("^T|^F",changeo.new$FUNCTIONAL,invert=T))/nrow(changeo.new)},
                   if(dispJunctionFr.list==T){do.call(paste, c(as.list(changeo.new$IN_FRAME), sep=", "))},
                   if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="T"))/nrow(changeo.new)},
                   if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="F"))/nrow(changeo.new)},
                   if(dispJunctionFr.ratio==T){length(grep("^T|^F",changeo.new$IN_FRAME,invert=T))/nrow(changeo.new)},
                   if(dispSeqID==T){do.call(paste, c(as.list(changeo.new$SEQUENCE_ID), sep=", "))},
                   if(dispTotalSeq==T){do.call(paste, c(as.list(unique(changeo.new$SEQUENCE_INPUT)), sep=", "))}))
                }
              }else if(length(which(a.dist[,a]<=tr))==0 && length(grep(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep=""), changeo.sub$JUNCTION_AA))>1){
                temp<-c(temp,uniqueCDR3.sub[a])

                changeo.new<-changeo.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",changeo.sub$JUNCTION_AA),perl=T),]
                
                tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[a]), sep=", ")), # shared CDR3 seq.
                 l, # CDR3 length
                 length(uniqueCDR3.sub[a]), # number_shared_CDR3
                 nrow(changeo.new), # number all sequences, belonging to clone
                 do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",changeo.new$JUNCTION_AA),perl=T))})), sep=", ")), # sequence count
                 V[i], # V_gene
                 do.call(paste, c(as.list(unique(changeo.new$V_CALL)), sep=", ")), # V_gene & allele
                 if(useJ==T && length(j)>0){do.call(paste, c(as.list(unique(j)), sep=", "))}else if(useJ==T && length(j)==0){"no J"}, # J_gene
                 do.call(paste, c(as.list(unique(changeo.new$J_CALL)), sep=", ")), # J_gene & allele

                 if(dispD==T){do.call(paste, c(as.list(unique(changeo.new$D_CALL)), sep=", "))},
                 if(dispCDR3aa==T){do.call(paste, c(as.list(changeo.new$JUNCTION_AA), sep=", "))},
                 if(dispCDR3nt==T){do.call(paste, c(as.list(changeo.new$JUNCTION), sep=", "))},
                 if(dispFunctionality.list==T){do.call(paste, c(as.list(changeo.new$FUNCTIONAL), sep=", "))},
                 if(dispFunctionality.ratio==T){length(grep("^T",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
                 if(dispFunctionality.ratio==T){length(grep("^F",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
                 if(dispFunctionality.ratio==T){length(grep("^T|^F",changeo.new$FUNCTIONAL,invert=T))/nrow(changeo.new)},
                 if(dispJunctionFr.list==T){do.call(paste, c(as.list(changeo.new$IN_FRAME), sep=", "))},
                 if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="T"))/nrow(changeo.new)},
                 if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="F"))/nrow(changeo.new)},
                 if(dispJunctionFr.ratio==T){length(grep("^T|^F",changeo.new$IN_FRAME,invert=T))/nrow(changeo.new)},
                 if(dispSeqID==T){do.call(paste, c(as.list(changeo.new$SEQUENCE_ID), sep=", "))},
                 if(dispTotalSeq==T){do.call(paste, c(as.list(unique(changeo.new$SEQUENCE_INPUT)), sep=", "))}))
              }
            }
          }
        }
      }
    } else { ### useJ=F
    changeo.sub<-changeo[grep(paste(V[i],"[!/*]",sep=""),changeo$V_CALL,perl=T),c('JUNCTION_AA', 'JUNCTION', 'V_CALL','J_CALL','D_CALL','FUNCTIONAL','SEQUENCE_ID', 'IN_FRAME', 'SEQUENCE_INPUT')]

    CDR3length<-apply(data.frame(changeo.sub$JUNCTION_AA),1,function(x){nchar(x)})
    uniqueCDR3length<-sort(as.numeric(unique(CDR3length)))
    if(length(uniqueCDR3length)==1){ ### only 1 CDR3 length
    temp<-vector()
    uniqueCDR3.sub<-unique(changeo.sub$JUNCTION_AA)
    a.dist<-adist(unique(uniqueCDR3.sub))
    a.dist[upper.tri(a.dist,diag = T)]<-NA
    tr<-floor(uniqueCDR3length*(1-identity))
    for(a in nrow(a.dist):1){
      if(length(which(a.dist[a,]<=tr))>0){
        if(a==nrow(a.dist) || sum(apply(data.frame(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),1,function(x){length(grep(gsub("[*]","-",x),gsub("[*]","-",temp),perl=T))}))<length(intersect(c(which(a.dist[a,]<=tr),a),which(!is.na(a.dist[a,]))))){
          temp<-c(temp,do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")))

          changeo.new<-changeo.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",changeo.sub$JUNCTION_AA),perl=T),]

          Jgene.uni<-unlist(lapply(changeo.new$J_CALL,function(x){strsplit(x, split=" |[*]|,|S")[[1]]}))
          Jgene.uni<-unique(Jgene.uni[grep("J",Jgene.uni)])

          tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")), # shared CDR3 seq.
           uniqueCDR3length, # CDR3 length
           length(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), # number_shared_CDR3
           nrow(changeo.new), # number all sequences, belonging to clone
           do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",changeo.new$JUNCTION_AA),perl=T))})), sep=", ")), # sequence count
           V[i], # V_gene
           do.call(paste, c(as.list(unique(changeo.new$V_CALL)), sep=", ")), # V_gene & allele
           if(length(Jgene.uni)>0){do.call(paste, c(as.list(Jgene.uni), sep=", "))}else{"No J"}, # J_gene
           do.call(paste, c(as.list(unique(changeo.new$J_CALL)), sep=", ")), # J_gene & allele

           if(dispD==T){do.call(paste, c(as.list(unique(changeo.new$D_CALL)), sep=", "))},
           if(dispCDR3aa==T){do.call(paste, c(as.list(changeo.new$JUNCTION_AA), sep=", "))},
           if(dispCDR3nt==T){do.call(paste, c(as.list(changeo.new$JUNCTION), sep=", "))},
           if(dispFunctionality.list==T){do.call(paste, c(as.list(changeo.new$FUNCTIONAL), sep=", "))},
           if(dispFunctionality.ratio==T){length(grep("^T",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
           if(dispFunctionality.ratio==T){length(grep("^F",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
           if(dispFunctionality.ratio==T){length(grep("^T|^F",changeo.new$FUNCTIONAL,invert=T))/nrow(changeo.new)},
           if(dispJunctionFr.list==T){do.call(paste, c(as.list(changeo.new$IN_FRAME), sep=", "))},
           if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="T"))/nrow(changeo.new)},
           if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="F"))/nrow(changeo.new)},
           if(dispJunctionFr.ratio==T){length(grep("^T|^F",changeo.new$IN_FRAME,invert=T))/nrow(changeo.new)},
           if(dispSeqID==T){do.call(paste, c(as.list(changeo.new$SEQUENCE_ID), sep=", "))},
           if(dispTotalSeq==T){do.call(paste, c(as.list(unique(changeo.new$SEQUENCE_INPUT)), sep=", "))}))
        }
      }else if(length(which(a.dist[,a]<=tr))==0 && length(grep(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep=""), changeo.sub$JUNCTION_AA))>1){
        temp<-c(temp,uniqueCDR3.sub[a])

        changeo.new<-changeo.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",changeo.sub$JUNCTION_AA),perl=T),]

        Jgene.uni<-unlist(lapply(changeo.new$J_CALL,function(x){strsplit(x, split=" |[*]|,|S")[[1]]}))
        Jgene.uni<-unique(Jgene.uni[grep("J",Jgene.uni)])

        tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[a]), sep=", ")), # shared CDR3 seq.
         l, # CDR3 length
         length(uniqueCDR3.sub[a]), # number_shared_CDR3
         nrow(changeo.new), # number all sequences, belonging to clone
         do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",changeo.new$JUNCTION_AA),perl=T))})), sep=", ")), # sequence count
         V[i], # V_gene
         do.call(paste, c(as.list(unique(changeo.new$V_CALL)), sep=", ")), # V_gene & allele
         if(length(Jgene.uni)>0){do.call(paste, c(as.list(Jgene.uni), sep=", "))}else{"No J"}, # J_gene
         do.call(paste, c(as.list(unique(changeo.new$J_CALL)), sep=", ")), # J_gene & allele

         if(dispD==T){do.call(paste, c(as.list(unique(changeo.new$D_CALL)), sep=", "))},
         if(dispCDR3aa==T){do.call(paste, c(as.list(changeo.new$JUNCTION_AA), sep=", "))},
         if(dispCDR3nt==T){do.call(paste, c(as.list(changeo.new$JUNCTION), sep=", "))},
         if(dispFunctionality.list==T){do.call(paste, c(as.list(changeo.new$FUNCTIONAL), sep=", "))},
         if(dispFunctionality.ratio==T){length(grep("^T",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
         if(dispFunctionality.ratio==T){length(grep("^F",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
         if(dispFunctionality.ratio==T){length(grep("^T|^F",changeo.new$FUNCTIONAL,invert=T))/nrow(changeo.new)},
         if(dispJunctionFr.list==T){do.call(paste, c(as.list(changeo.new$IN_FRAME), sep=", "))},
         if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="T"))/nrow(changeo.new)},
         if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="F"))/nrow(changeo.new)},
         if(dispJunctionFr.ratio==T){length(grep("^T|^F",changeo.new$IN_FRAME,invert=T))/nrow(changeo.new)},
         if(dispSeqID==T){do.call(paste, c(as.list(changeo.new$SEQUENCE_ID), sep=", "))},
         if(dispTotalSeq==T){do.call(paste, c(as.list(unique(changeo.new$SEQUENCE_INPUT)), sep=", "))}))
      }
    }
  }else{ ### several CDR3 length
  for(l in uniqueCDR3length){
    temp<-vector()
    uniqueCDR3.sub<-unique(changeo.sub$JUNCTION_AA[which(CDR3length==l)])
    a.dist<-adist(unique(uniqueCDR3.sub))
    a.dist[upper.tri(a.dist,diag = T)]<-NA
    tr<-floor(l*(1-identity))
    for(a in nrow(a.dist):1){
      if(length(which(a.dist[a,]<=tr))>0){
        if(a==nrow(a.dist) || sum(apply(data.frame(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),1,function(x){length(grep(gsub("[*]","-",x),gsub("[*]","-",temp),perl=T))}))<length(intersect(c(which(a.dist[a,]<=tr),a),which(!is.na(a.dist[a,]))))){
          temp<-c(temp,do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")))

          changeo.new<-changeo.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")), sep="|")),gsub("[*]","-",changeo.sub$JUNCTION_AA),perl=T),]

          Jgene.uni<-unlist(lapply(changeo.new$J_CALL,function(x){strsplit(x, split=" |[*]|,|S")[[1]]}))
          Jgene.uni<-unique(Jgene.uni[grep("J",Jgene.uni)])

          tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), sep=", ")), # shared CDR3 seq.
           l, # CDR3 length
           length(uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]), # number_shared_CDR3
           nrow(changeo.new), # number all sequences, belonging to clone
           do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[c(which(a.dist[a,]<=tr),a)]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",changeo.new$JUNCTION_AA),perl=T))})), sep=", ")), # sequence count
           V[i], # V_gene
           do.call(paste, c(as.list(unique(changeo.new$V_CALL)), sep=", ")), # V_gene & allele
           if(length(Jgene.uni)>0){do.call(paste, c(as.list(Jgene.uni), sep=", "))}else{"No J"}, # J_gene
           do.call(paste, c(as.list(unique(changeo.new$J_CALL)), sep=", ")), # J_gene & allele

           if(dispD==T){do.call(paste, c(as.list(unique(changeo.new$D_CALL)), sep=", "))},
           if(dispCDR3aa==T){do.call(paste, c(as.list(changeo.new$JUNCTION_AA), sep=", "))},
           if(dispCDR3nt==T){do.call(paste, c(as.list(changeo.new$JUNCTION), sep=", "))},
           if(dispFunctionality.list==T){do.call(paste, c(as.list(changeo.new$FUNCTIONAL), sep=", "))},
           if(dispFunctionality.ratio==T){length(grep("^T",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
           if(dispFunctionality.ratio==T){length(grep("^F",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
           if(dispFunctionality.ratio==T){length(grep("^T|^F",changeo.new$FUNCTIONAL,invert=T))/nrow(changeo.new)},
           if(dispJunctionFr.list==T){do.call(paste, c(as.list(changeo.new$IN_FRAME), sep=", "))},
           if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="T"))/nrow(changeo.new)},
           if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="^F"))/nrow(changeo.new)},
           if(dispJunctionFr.ratio==T){length(grep("frame",changeo.new$IN_FRAME,invert=T))/nrow(changeo.new)},
           if(dispSeqID==T){do.call(paste, c(as.list(changeo.new$SEQUENCE_ID), sep=", "))},
           if(dispTotalSeq==T){do.call(paste, c(as.list(unique(changeo.new$SEQUENCE_INPUT)), sep=", "))}))
        }
      }else if(length(which(a.dist[,a]<=tr))==0 && length(grep(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep=""), changeo.sub$JUNCTION_AA))>1){
        temp<-c(temp,uniqueCDR3.sub[a])

        changeo.new<-changeo.sub[grep(do.call(paste, c(as.list(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")), sep="|")),gsub("[*]","-",changeo.sub$JUNCTION_AA),perl=T),]

        Jgene.uni<-unlist(lapply(changeo.new$J_CALL,function(x){strsplit(x, split=" |[*]|,|S")[[1]]}))
        Jgene.uni<-unique(Jgene.uni[grep("J",Jgene.uni)])

        tempout<-rbind(tempout,c(do.call(paste, c(as.list(uniqueCDR3.sub[a]), sep=", ")), # shared CDR3 seq.
         l, # CDR3 length
         length(uniqueCDR3.sub[a]), # number_shared_CDR3
         nrow(changeo.new), # number all sequences, belonging to clone
         do.call(paste, c(as.list(apply(data.frame(paste("^",gsub("[*]","-",uniqueCDR3.sub[a]),"$",sep="")),1,function(x){length(grep(x,gsub("[*]","-",changeo.new$JUNCTION_AA),perl=T))})), sep=", ")), # sequence count
         V[i], # V_gene
         do.call(paste, c(as.list(unique(changeo.new$V_CALL)), sep=", ")), # V_gene & allele
         if(length(Jgene.uni)>0){do.call(paste, c(as.list(Jgene.uni), sep=", "))}else{"No J"}, # J_gene
         do.call(paste, c(as.list(unique(changeo.new$J_CALL)), sep=", ")), # J_gene & allele

         if(dispD==T){do.call(paste, c(as.list(unique(changeo.new$D_CALL)), sep=", "))},
         if(dispCDR3aa==T){do.call(paste, c(as.list(changeo.new$JUNCTION_AA), sep=", "))},
         if(dispCDR3nt==T){do.call(paste, c(as.list(changeo.new$JUNCTION), sep=", "))},
         if(dispFunctionality.list==T){do.call(paste, c(as.list(changeo.new$FUNCTIONAL), sep=", "))},
         if(dispFunctionality.ratio==T){length(grep("^T",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
         if(dispFunctionality.ratio==T){length(grep("^F",changeo.new$FUNCTIONAL))/nrow(changeo.new)},
         if(dispFunctionality.ratio==T){length(grep("^T|^F",changeo.new$FUNCTIONAL,invert=T))/nrow(changeo.new)},
         if(dispJunctionFr.list==T){do.call(paste, c(as.list(changeo.new$IN_FRAME), sep=", "))},
         if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="T"))/nrow(changeo.new)},
         if(dispJunctionFr.ratio==T){length(which(changeo.new$IN_FRAME=="F"))/nrow(changeo.new)},
         if(dispJunctionFr.ratio==T){length(grep("^T|^F",changeo.new$IN_FRAME,invert=T))/nrow(changeo.new)},
         if(dispSeqID==T){do.call(paste, c(as.list(changeo.new$SEQUENCE_ID), sep=", "))},
         if(dispTotalSeq==T){do.call(paste, c(as.list(unique(changeo.new$SEQUENCE_INPUT)), sep=", "))}))
      }
    }
  }
}
}
return(data.frame(tempout,stringsAsFactors = F))
}
clonRel<-do.call(rbind.data.frame, clonelist)
stopCluster(cl)

if(is.data.frame(clonRel) && nrow(clonRel)>0){
  colnames(clonRel)<-c("unique_CDR3_sequences_AA",
   "CDR3_length_AA",
   "number_of_unique_sequences",
   "total_number_of_sequences",
   "sequence_count_per_CDR3",
   "V_gene",
   "V_gene_and_allele",
   "J_gene",
   "J_gene_and_allele",
   if(dispD==T){"D_gene"},
   if(dispCDR3aa==T){"all_CDR3_sequences_AA"},
   if(dispCDR3nt==T){"all_CDR3_sequences_nt"},
   if(dispFunctionality.list==T){"Functionality_all_sequences"},
   if(dispFunctionality.ratio==T){c("Func_productive_sequences","Func_unproductive_sequences", "Func_unknown")},
   if(dispJunctionFr.list==T){"Junction_frame_all_sequences"},
   if(dispJunctionFr.ratio==T){c("JF_in_frame","JF_out_of_frame","JF_unknown")},
   if(dispSeqID==T){"Sequence_IDs"},
   if(dispTotalSeq==T){"Total_sequences_nt"})

  clonRel<-clonRel[duplicated(clonRel[,"unique_CDR3_sequences_AA"])==F,]
  if(length(intersect(which(clonRel[,"total_number_of_sequences"]==1),
    which(clonRel[,"sequence_count_per_CDR3"]==1)))>0){
    clonRel<-clonRel[-intersect(which(clonRel[,"total_number_of_sequences"]==1),
      which(clonRel[,"sequence_count_per_CDR3"]==1)),]
  }
} else {
  clonRel<-NULL
  cat("... no clones\n")
}
return(data.frame(clonRel,check.names=F,stringsAsFactors = F))
}

clonesID <- function(clones) {
  seqID = strsplit(clones$Sequence_IDs, ", ")
  names(seqID) <- seq(1, length(seqID))
  clone_df <- list()
  for (i in 1:length(seqID)) {
    clone_df[[i]] <- data.frame(SEQUENCE_ID = seqID[[i]], CLONE = names(seqID)[i])
  }

  clone_df <- do.call(rbind, clone_df)
  clone_df <- apply(clone_df, 2, as.character)

  if (length(which(duplicated(clone_df[, 1]))) > 0) {
    add.out <- vector()
    dupli.set <- clone_df[which(duplicated(clone_df[, 1]) == TRUE), 1]
    for (i in 1:length(dupli.set)) {
      add.out <- rbind(add.out, c(dupli.set[i], do.call(paste, c(as.list(unique(clone_df[which(clone_df[, 1] == dupli.set[i]), 2])), sep = "_"))))
    }
    clone_df <- clone_df[-which(clone_df[, 1] %in% dupli.set), ]
    clone_df <- rbind(clone_df, add.out)
  }

  clone_df <- as.data.frame(clone_df, stringsAsFactors = FALSE)
  return(clone_df)
}

changeo <- readr::read_tsv(opt$changeo, col_types = readr::cols())
changeo <- changeo %>% filter(FUNCTIONAL == 'TRUE')
changeo$BARCODE <- gsub('_c.*', '', changeo$SEQUENCE_ID)
changeo <- changeo %>% select(BARCODE, everything())

clones <- findClones(changeo)
clone_df <- clonesID(clones)

changeo_new <- left_join(changeo, clone_df, by = 'SEQUENCE_ID')
unassigned_IGH <- which(changeo_new$LOCUS == 'IGH' & changeo_new$FUNCTIONAL == TRUE & !is.na(changeo_new$JUNCTION) & is.na(changeo_new$CLONE))
suppressWarnings(start_unassigned <- max(as.numeric(changeo_new$CLONE), na.rm = TRUE)+1)
library(magrittr)
changeo_new[unassigned_IGH, ] %<>% mutate(CLONE = as.character(seq(start_unassigned, start_unassigned+length(unassigned_IGH)-1)))
changeo_new <- split(changeo_new, changeo_new$BARCODE)
changeo_new <- parallel::mclapply(changeo_new, function(x) {
  if(any(!is.na(x$CLONE))){
    if(length(which(is.na(x$CLONE))) < 1){
      clone = as.character(x$CLONE)             
    } else {
      clone = paste(as.character(x$CLONE[-which(is.na(x$CLONE))]), collapse ='_')
    }
    x$CLONE = clone
  } else {
    x$CLONE = NA
  }
  return(x)
}, mc.cores = parallel::detectCores())
names(changeo_new) <- NULL
changeo_new <- do.call(rbind, changeo_new)
write.table(changeo_new, opt$output, col.names = TRUE, row.names = FALSE, quote = FALSE, sep ='\t')