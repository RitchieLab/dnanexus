# VCF QC (DNAnexus Platform App)

## What does this applet do?

This app uses [GATK ApplyRecalibration](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php) and [GATK VariantFiltration](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php) to apply filters to VCF files.

## What data are required for this app to run?

This app requires:
1. The array of VCF (variant calling format) files (`*.vcf.gz`) and their associated index files (`*.vcf.gz.tbi`). VCF is the standard format for variant calling.
2. The GATK 3 tar ball
3. The reference genome in FASTA format (`*.fasta` or `*.fasta`), the Picard FASTA index (`*.dict`), and the satools FASTA index (`*.fasta.fai`)

The user can also provide the BED interval of exon target and other options for variant filtering.

## What does this applet output?

This app outputs filtered VCF files (`*.vcf.gz`) and their associated index files (`*.vcf.gz.tbi`).

