#!/usr/bin/env Rscript

library(gRodon)
library(Biostrings)

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
