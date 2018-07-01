# VCF QC (DNAnexus Platform App)

## What does this applet do?

This app uses [GATK ApplyRecalibration](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php) and [GATK VariantFiltration](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php) to apply filters to VCF files.

## What data are required for this app to run?

This app requires the array of VCF (variant calling format) files (`*.vcf.gz`) and their associated index files (`*.vcf.gz.tbi`). VCF is the standard format for variant calling. The user can also provide the BED interval of exon target and other options for variant filtering.


## What does this applet output?

This app outputs filtered VCF files (`*.vcf.gz`) and their associated index files (`*.vcf.gz.tbi`).

