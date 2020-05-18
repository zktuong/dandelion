#!/usr/bin/env Rscript
# @Author: kt16
# @Date:   2019-11-25 22:59:25
# @Last Modified by:   kt16
# @Last Modified time: 2020-02-16 20:49:02
# @Notes:	plot the resulting network
options(warn=-1)
suppressMessages(library(optparse))
option_list = list(
        make_option(c("-n", "--name"), type = "character", default = NULL,
              help = "name of sample", metavar = "character"),
        make_option(c("-v", "--vertex"), type = "character", default = NULL,
              help = "name of sample", metavar = "character"),
        make_option(c("-e", "--edge"), type = "character", default = NULL,
              help = "name of sample", metavar = "character"),        
        make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "name of sample", metavar = "character"),
        make_option(c("-l", "--layout"), type = "character", default = NULL,
              help = "name of sample", metavar = "character")
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$name)){
  print_help(opt_parser)
  stop("Argument must be supplied (name of sample).n", call. = TRUE)}
if (is.null(opt$vertex)){
  print_help(opt_parser)
  stop("Argument must be supplied (vertex file).n", call. = TRUE)}
if (is.null(opt$edge)){
  print_help(opt_parser)
  stop("Argument must be supplied (edge file).n", call. = TRUE)}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Argument must be supplied (out file).n", call. = TRUE)}
if (is.null(opt$layout)){
  print_help(opt_parser)
  stop("Argument must be supplied (layout file).n", call. = TRUE)}

## program ##
nodes <- as.data.frame(readr::read_tsv(opt$vertex, col_names = FALSE, col_types = readr::cols()))
edge <- as.data.frame(readr::read_tsv(opt$edge, col_names = FALSE, col_types = readr::cols()))

suppressMessages(library(igraph))
g <- graph.empty(n = 0, directed = FALSE)
freq <- as.numeric(nodes[, 2])
# max_freq <- (sum(freq))
# frequency <- freq/max_freq
g <- igraph::add.vertices(g, length(nodes[, 1]), name = as.character(nodes[, 1]))	
names <- V(g)$name
ids <- 1:length(names)
names(ids) <- names
from <- as.character(edge[, 1])
to <- as.character(edge[, 2])
edges <- matrix(c(ids[from], ids[to]), nc = 2)
g <- add.edges(g, t(edges), weight = 1)
V(g)$label <- V(g)$name
# V(g)$size <- frequency

V(g)$size <- 0.05
V(g)$label.cex <- 1e-100
E(g)$size <- (0.8/(edge[, 3]+1))
mat <- match(from, as.character(nodes[, 1]))

if (!file.exists(opt$layout)) {
  cat('writing new layout file\n')
  # layout <- layout_with_graphopt(g, niter = 800, charge = 0.00001, spring.constant = 2)
  layout <- layout_with_fr(g)
  write.table(layout, file = opt$layout, sep = "\t", row.names = FALSE, col.names = FALSE)
} else {
  layout <- as.matrix(readr::read_tsv(opt$layout, col_names = FALSE, col_types = readr::cols()))
} 

colour_i <- as.numeric(factor(nodes[, 3], levels = c('IGHA1', 'IGHA2', 'IGHD', 'IGHE', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHM', 'NA', NA, 'IGKC', 'IGLC1', 'IGLC2', 'IGLC3', 'IGLC4', 'IGLC5', 'IGLC6', 'IGLC7')))
colour_i1 <- c('#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f', '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#f2f2f2', '#f2f2f2', '#f2f2f2', '#f2f2f2', '#f2f2f2', '#f2f2f2', '#f2f2f2', '#f2f2f2', '#f2f2f2', '#f2f2f2')
V(g)$color <- colour_i1[colour_i]
V(g)$frame.color <- colour_i1[colour_i]
edge_col_i <- colour_i1[colour_i[mat]]
E(g)$color <- edge_col_i 
g1 <- simplify(g, edge.attr.comb = 'first')

pdf(paste0(opt$out, '/isotype_network_R.pdf'), h = 5, w = 5)
plot(g1, layout = layout, main = opt$name, edge.width = E(g1)$size)
legend('topleft',legend=c('IGHA1', 'IGHA2', 'IGHD', 'IGHE', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHM'), pt.cex=1, col=colour_i1, pch=16, bty = "n", cex = 0.4)
dev.off()

colour_d <- as.numeric(factor(nodes[, 4], levels = c('A29', 'A31', 'A35', 'A36', 'A37')))
colour_d1 <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd')
V(g)$color <- colour_d1[colour_d]
V(g)$frame.color <- colour_d1[colour_d]
edge_col_d <- colour_d1[colour_d[mat]]
E(g)$color <- edge_col_d  
g1 <- simplify(g, edge.attr.comb = 'first')

pdf(paste0(opt$out, '/sample_network_R.pdf'), h = 5, w = 5)
plot(g1, layout = layout, main = opt$name, edge.width = E(g1)$size)
legend('topleft',legend=c('A29', 'A31', 'A35', 'A36', 'A37'), pt.cex=1, col=colour_d1, pch=16, bty = "n", cex = 0.4)
dev.off()

colour_t <- as.numeric(factor(nodes[, 5], levels = c("CAE", "TCL", "SCL", "TLN", "SPL", "LNG", "MLN", "BMA", "LIV", "DUO", "ILE", "KID", "OES", "OME", "DKM", "THY", "TIL", "SKM")))
colour_t1 <- c('#4E79A7','#A0CBE8','#F28E2B','#FFBE7D','#59A14F','#8CD17D','#B6992D','#F1CE63','#499894','#86BCB6','#E15759','#FF9D9A','#79706E','#BAB0AC','#D37295','#FABFD2','#B07AA1','#D4A6C8')
V(g)$color <- colour_t1[colour_t]
V(g)$frame.color <- colour_t1[colour_t]
edge_col_t <- colour_t1[colour_t[mat]]
E(g)$color <- edge_col_t  
g1 <- simplify(g, edge.attr.comb = 'first')

pdf(paste0(opt$out, '/tissue_network_R.pdf'), h = 5, w = 5)
plot(g1, layout = layout, main = opt$name, edge.width = E(g1)$size)
legend('topleft',legend=c("CAE", "TCL", "SCL", "TLN", "SPL", "LNG", "MLN", "BMA", "LIV", "DUO", "ILE", "KID", "OES", "OME", "DKM", "THY", "TIL", "SKM"), pt.cex=1,col=colour_t1, pch=16, bty = "n", cex = 0.4)
dev.off()

colour_c <- as.numeric(factor(nodes[, 6], levels = c("IGH", "IGK", "IGL")))
colour_c1 <- c('#c7c7c7', '#1F77B4', '#FF7F0E')
V(g)$color <- colour_c1[colour_c]
V(g)$frame.color <- colour_c1[colour_c]
edge_col_c <- colour_c1[colour_c[mat]]
E(g)$color <- edge_col_c
g1 <- simplify(g, edge.attr.comb = 'first')

pdf(paste0(opt$out, '/locus_network_R.pdf'), h = 5, w = 5)
plot(g1, layout = layout, main = opt$name, edge.width = E(g1)$size)
legend('topleft',legend=c("IGH", "IGK", "IGL"), pt.cex=1,col=colour_c1, pch=16, bty = "n", cex = 0.4)
dev.off()

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

colour_cl <- as.numeric(factor(nodes[, 7]))
colour_cl1 <- gg_color_hue(length(unique(colour_cl)))
V(g)$color <- colour_cl1[colour_cl]
V(g)$frame.color <- colour_cl1[colour_cl]
edge_col_cl <- colour_cl1[colour_cl[mat]]
E(g)$color <- edge_col_cl
g1 <- simplify(g, edge.attr.comb = 'first')
pdf(paste0(opt$out, '/clone_network_R.pdf'), h = 5, w = 5)
plot(g1, layout = layout, main = opt$name, edge.width = E(g1)$size)
dev.off()
