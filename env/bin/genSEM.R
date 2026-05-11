#!/usr/bin/env Rscript

# remotes::install_github("GenomicSEM/GenomicSEM")

suppressPackageStartupMessages({
  library(GenomicSEM)
})


args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 16) stop("Usage: Rscript genSEM.R <pheno1_id> <pheno1_sumstats> <pheno2_id> <pheno2_sumstats> <pheno3_id> <pheno3_sumstats> <t1_cases> <t1_controls> <t2_cases> <t2_controls> <t3_cases> <t3_controls> <t1_pop_prev> <t2_pop_prev> <t3_pop_prev> <out_dir>")

pheno1_id       <- args[1]
pheno1_sumstats <- args[2]
pheno2_id       <- args[3]
pheno2_sumstats <- args[4]
pheno3_id       <- args[5]
pheno3_sumstats <- args[6]
t1_cases        <- as.numeric(args[7])
t1_controls     <- as.numeric(args[8])
t2_cases        <- as.numeric(args[9])
t2_controls     <- as.numeric(args[10])
t3_cases        <- as.numeric(args[11])
t3_controls     <- as.numeric(args[12])
t1_pop_prev     <- as.numeric(args[13])
t2_pop_prev     <- as.numeric(args[14])
t3_pop_prev     <- as.numeric(args[15])
out_dir         <- args[16]

# analytics pipeline with paths
# *2 -> assemble as data pipeline

files <- c(
  pheno1_sumstats,
  pheno2_sumstats,
  pheno3_sumstats
)

trait_names <- c(
  pheno1_id,
  pheno2_id,
  pheno3_id
)

hm3 <- "ref/ldsc/eur_w_ld_chr/w_hm3.snplist"
ld <- "ref/ldsc/eur_w_ld_chr/"
wld <- "ref/ldsc/weights_hm3_no_hla/"

munge_dir <- file.path(out_dir, "munge")
ldsc_dir <- file.path(out_dir, "ldsc")

dir.create(munge_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(ldsc_dir, recursive=TRUE, showWarnings=FALSE)

N <- c(
  4 * t1_cases * t1_controls / (t1_cases + t1_controls),
  4 * t2_cases * t2_controls / (t2_cases + t2_controls),
  4 * t3_cases * t3_controls / (t3_cases + t3_controls)
)

old_wd <- getwd()
setwd(munge_dir)

munge(
  files=normalizePath(file.path(old_wd, files)),
  hm3=normalizePath(file.path(old_wd, hm3)),
  trait.names=trait_names,
  N=N
)

setwd(old_wd)

traits <- file.path(
  munge_dir,
  paste0(trait_names, ".sumstats.gz")
)

sample.prev <- c(
  t1_cases / (t1_cases + t1_controls),
  t2_cases / (t2_cases + t2_controls),
  t3_cases / (t3_cases + t3_controls)
)

population.prev <- c(
  t1_pop_prev,
  t2_pop_prev,
  t3_pop_prev
)

LDSC_INT <- ldsc(
  traits=traits,
  sample.prev=sample.prev,
  population.prev=population.prev,
  ld=ld,
  wld=wld,
  trait.names=trait_names
)

saveRDS(
  LDSC_INT,
  file=file.path(ldsc_dir, "LDSC_INT.rds")
)

write.table(
  LDSC_INT$S,
  file=file.path(ldsc_dir, "S_matrix.tsv"),
  sep="\t",
  quote=FALSE,
  col.names=NA
)

write.table(
  LDSC_INT$V,
  file=file.path(ldsc_dir, "V_matrix.tsv"),
  sep="\t",
  quote=FALSE,
  col.names=NA
)

write.table(
  LDSC_INT$I,
  file=file.path(ldsc_dir, "I_matrix.tsv"),
  sep="\t",
  quote=FALSE,
  col.names=NA
)

write.table(
  LDSC_INT$N,
  file=file.path(ldsc_dir, "N_matrix.tsv"),
  sep="\t",
  quote=FALSE,
  col.names=NA
)


capture.output(
  print(LDSC_INT),
  file=file.path(ldsc_dir, "LDSC_results.txt")
)


# exploratory factor analysis (EFA)




