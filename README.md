# GenomicSEM pipeline 

---

### Run QC on SLURM cluster

````bash
ssh <user>@falconlogin.cf.ac.uk \
'cd /shared/home1/<user>/genSEM && \
apptainer exec env/genSEM.sif python - \
--pheno-id AD \
--sumstats GCST90027158_buildGRCh38.tsv \
--out-dir /shared/home1/<user>/genSEM/out \
--maf 0.01 \
--snp-col variant_id \
--a1-col effect_allele \
--a2-col other_allele \
--beta-col beta \
--se-col standard_error \
--af_col effect_allele_frequency \
--p-col p_value \
--pos-col base_pair_location \
--chr-col chromosome \
--genome_build GRCh38 \
--remove_mhc \
--n_cases 111326 \
--n_controls 677663 \
--falcon-user <user>' < qc_gwas.py
```

---

# Run LDSC (genSEM)

´´´´bash
ssh <user>@falconlogin.cf.ac.uk \
'cd /shared/home1/<user>/genSEM && \
apptainer exec env/genSEM.sif Rscript - \
AD out/QC/AD/AD.tsv \
SCZ out/QC/SCZ/SCZ.tsv \
MDD out/QC/MDD/MDD.tsv \
111326 677663 \
67390 94015 \
135458 344901 \
0.07 0.01 0.15 \
out' < genSEM.R
´´´ 


---

## Run CFA.R

ssh <user>@falconlogin.cf.ac.uk \
'cd /shared/home1/<user>/genSEM && \
apptainer exec env/genSEM.sif Rscript - \
AD SCZ MDD out' < CFA.R

---

## Run TFA on AD loading with MDD+SCZ latent factor

ssh <user>@falconlogin.cf.ac.uk \
'cd /shared/home1/<user>/genSEM && \
apptainer exec env/genSEM.sif Rscript - \
AD SCZ MDD out' < TFA.R

