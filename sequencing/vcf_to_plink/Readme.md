# Convert VCF file to PLINK file (DNAnexus Platform App)

## What does this app do?

This app uses PLINK 1.9 to convert VCF files (`*.vcf.gz`) to PLINK bed/bim/fam files (`*.bed`, `*.bed`, `*.fam`).

## What data are required for this app to run?

This app requires an equal number of VCF files (`*.vcf.gz`) and their associated index files (`*.vcf.gz.tbi`).

## What does this app output?

This app outputs an equal number of binary genotype files (`*.bed`), extend mapping files (`*.bim`), and FAM files (`*.fam`) as that of VCF files. PLINK is a standard format for whole-genome association analysis of genotype/phenotype data.