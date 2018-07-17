# Calculate PCs from VCF (DNAnexus Platform App)

## What does this app do?

This app uses PLINK 1.9 and eigenstrat 6.0 to calculate principal components from VCF (`*.vcf.gz`) or PLINK files (`*.bed`, `*.bim`, `*.fam`). The app can process both types of data at the same time.

## What data are required for this app to run?

This app takes the following inputs:

For VCF data, this app requires an equal number of VCF (variant calling format) files (`*.vcf.gz`) and their associated index files (`*.vcf.gz.tbi`)

For PLINK, this app requires an equal number of binary genotype files (`*.bed`), extend mapping files (`*.bim`), and FAM files (`*.fam`) -- first 6 columns of the pedigree file. PLINK is a standard format for whole-genome association analysis of genotype/phenotype data.


## What does this app output?

This app outputs the eigenvectors, eigenvalues, TW Statistics (Tracy-Widom statistics for statistically significant principal components), and samples excluded from analysis.

