#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(gRodon)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: 08_run_grodon_one.R <genome_id> <cds_fasta> <he_ids_file>")
}

genome_id <- args[[1]]
cds_fasta <- args[[2]]
he_file <- args[[3]]

genes <- readDNAStringSet(cds_fasta)
he_ids <- readLines(he_file, warn = FALSE)
highly_expressed <- names(genes) %in% he_ids
n_he <- sum(highly_expressed)

if (length(genes) == 0) {
  stop("No CDS found")
}
if (n_he < 10) {
  stop(paste("Too few highly expressed genes:", n_he))
}

res <- suppressWarnings(predictGrowth(genes, highly_expressed))
fields <- c(
  "RESULT",
  genome_id,
  length(genes),
  n_he,
  res$CUBHE,
  res$GC,
  res$GCdiv,
  res$ConsistencyHE,
  res$CUB,
  res$CPB,
  res$FilteredSequences,
  res$d,
  res$LowerCI,
  res$UpperCI
)
cat(paste(fields, collapse = "\t"), "\n")
