#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
})

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript CFA.R <pheno1_id> <pheno2_id> <pheno3_id> <pheno4_id> <out_dir>")
}

pheno1_id <- args[1]
pheno2_id <- args[2]
pheno3_id <- args[3]
pheno4_id <- args[4]
out_dir <- args[5]

trait_names <- c(
  pheno1_id,
  pheno2_id,
  pheno3_id,
  pheno4_id
)

hm3 <- "ref/ldsc/eur_w_ld_chr/w_hm3.snplist"
ld <- "ref/ldsc/eur_w_ld_chr/"
wld <- "ref/ldsc/weights_hm3_no_hla/"

# LOAD LDSC output
ldsc <- file.path(out_dir, "ldsc", "LDSC_INT.rds")
LDSC_INT <- readRDS(ldsc)


# exploratory factor analysis (EFA) # model 3 
# Load AD, MDD and SCZ into a single latent factor 
# F1 =~ AD + SCZ + MDD

# contruct multiv LDSC matrix
# CFA

# IntResults <- usermodel(covstruc=covstruc, model=INT.model)
covstruc <- LDSC_INT
INT.model <- paste0("F1 =~ ", paste(trait_names, collapse=" + "))

# IntResults <- usermodel(covstruc=covstruc, model=INT.model)
# IntResults <- usermodel(
#  covstruc=covstruc,
#  model=INT.model,
#  std.lv=TRUE
# )


sem_dir <- file.path(out_dir, "sem")
dir.create(sem_dir, recursive=TRUE, showWarnings=FALSE)

sink(file.path(sem_dir, "CFA_console_output.txt"))

IntResults <- usermodel(
  covstruc=covstruc,
  model=INT.model,
  std.lv=TRUE
)

sink()

if (!is.null(IntResults)) {
  saveRDS(IntResults, file=file.path(sem_dir, "CFA_results.rds"))
  write.table(IntResults, file=file.path(sem_dir, "CFA_results.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}

print("Model fit")
IntResults$modelfit