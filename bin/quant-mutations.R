#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("shazam"))

NPROC <- parallel::detectCores()

opt_list <- list(make_option(c("-d", "--db"), dest="DB", type = 'character', help="input tsv file"),
	make_option(c("-r", "--reg"), dest="REG", type = 'character', help = 'region definition', ),
	make_option(c("-m", "--mut"), dest="MUT", type = 'character', help = 'mutation definition'),
	make_option(c("-f", "--freq"), dest="FREQ", type = 'logical', help = 'frequency'),
	make_option(c("-c", "--comb"), dest="COMB", type = 'logical', help = 'combine'),
	make_option(c("-o", "--out"), dest="OUT", type = 'character', help = 'out file name'))

# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))
db <- airr::read_rearrangement(opt$DB)
if (!'germline_alignment_d_mask' %in% colnames(db)){
	stop('Please run CreateGermlines.py first.')
}

if (opt$REG == 'NULL'){
	reg = NULL	
} else{
	reg = reg
}

if (opt$MUT == 'NULL'){
	mut = NULL	
} else{
	mut = mut
}

res = observedMutations(db,
	sequenceColumn = "sequence_alignment",
    germlineColumn = "germline_alignment_d_mask",
	regionDefinition=reg,
	mutationDefinition=mut,
	frequency=opt$FREQ, 
	combine=opt$COMB,
	nproc=NPROC)
# Write genotyped data
alakazam::writeChangeoDb(db, opt$OUT)