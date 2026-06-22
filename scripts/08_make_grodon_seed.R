#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(gRodon)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "data/external/grodon"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

path_to_genome <- system.file(
  "extdata",
  "GCF_000349925.2_ASM34992v2_cds_from_genomic.fna.gz",
  package = "gRodon"
)
genes <- readDNAStringSet(path_to_genome)
is_ribo <- grepl(
  "^(?!.*(methyl|hydroxy)).*0S ribosomal protein",
  names(genes),
  ignore.case = TRUE,
  perl = TRUE
)
ribosomal_cds <- genes[is_ribo]
names(ribosomal_cds) <- paste0("gRodon_seed_ribo_", seq_along(ribosomal_cds))
ribosomal_proteins <- translate(ribosomal_cds, if.fuzzy.codon = "X")

writeXStringSet(ribosomal_cds, file.path(out_dir, "grodon_ribosomal_seed_cds.fna"))
writeXStringSet(ribosomal_proteins, file.path(out_dir, "grodon_ribosomal_seed.faa"))
cat("Wrote", length(ribosomal_cds), "ribosomal protein seeds to", out_dir, "\n")
