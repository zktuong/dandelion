.libPaths()
# install.packages(c("BiocManager"), repos = "https://cloud.r-project.org")
# BiocManager::install(c("edgeR", "Biostrings", "GenomicAlignments", "IRanges"))
# install.packages(c("shazam", "alakazam", "tigger", "scoper"), repos = "https://cloud.r-project.org")

pkgs <- c("shazam", "alakazam", "tigger", "airr", "optparse", "edgeR", "scoper")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Optional R package '%s' not installed, skipping...", pkg))
  }
}
quit(status = 0)
