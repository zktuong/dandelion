if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings","GenomicAlignments", "IRanges"))
install.packages(c('shazam', 'alakazam', 'tigger', 'airr', 'optparse'))