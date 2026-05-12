#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
})

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript CFA.R <pheno1_id> <pheno2_id> <pheno3_id> <out_dir>")
}

pheno1_id <- args[1]
pheno2_id <- args[2]
pheno3_id <- args[3]
out_dir <- args[4]

trait_names <- c(
  pheno1_id, # AD
  pheno2_id, # SCZ
  pheno3_id  # MDD
)

hm3 <- "ref/ldsc/eur_w_ld_chr/w_hm3.snplist"
ld <- "ref/ldsc/eur_w_ld_chr/"
wld <- "ref/ldsc/weights_hm3_no_hla/"

# LOAD LDSC output
ldsc <- file.path(out_dir, "ldsc", "LDSC_INT.rds")
LDSC_INT <- readRDS(ldsc)

# INT.model <- "
# Psych =~ SCZ + MDD
# Psych ~~ AD  # Correlation between the factor and AD
# "

INT.model <- paste0(
  "F1 =~ L1*", pheno2_id, " + L1*", pheno3_id, "\n",
  "F1 ~~ ", pheno1_id, "\n",
  pheno1_id, " ~~ ", pheno1_id
)

sem_dir <- file.path(out_dir, "sem")
dir.create(sem_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(sem_dir, "Psych_AD_TFA_log.txt")
sink(log_file)

IntResults <- usermodel(
  covstruc = LDSC_INT, 
  model = INT.model, 
  std.lv = TRUE
)

sink()

if (!is.null(IntResults)) {
  saveRDS(IntResults, file = file.path(sem_dir, "Psych_AD_results.rds"))
    write.table(
    IntResults$results, 
    file = file.path(sem_dir, "Psych_AD_results.tsv"), 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE
  )
  
  cat("\n--- Model Fit Indices ---\n")
  print(IntResults$modelfit)
} else {
  cat("\nError: Model failed to produce results. Check the log file.\n")
}