#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
repo_dir <- if (length(args) >= 1) args[[1]] else "/home/data/haoyu/software/gRodon2"

if (!dir.exists(repo_dir)) {
  stop("gRodon2 source directory does not exist: ", repo_dir)
}

r_bin <- file.path(R.home("bin"), "R")
cmd_args <- c("CMD", "INSTALL", repo_dir)
status <- system2(r_bin, cmd_args)
if (!identical(status, 0L)) {
  stop("R CMD INSTALL failed for: ", repo_dir)
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(gRodon)
})

path_to_genome <- system.file(
  "extdata",
  "GCF_000349925.2_ASM34992v2_cds_from_genomic.fna.gz",
  package = "gRodon"
)
genes <- readDNAStringSet(path_to_genome)
highly_expressed <- grepl(
  "^(?!.*(methyl|hydroxy)).*0S ribosomal protein",
  names(genes),
  ignore.case = TRUE,
  perl = TRUE
)
print(predictGrowth(genes, highly_expressed))
