## List of apps and descriptions of them in this directory. 

* bcftools_view: BCFTools View
  - Calls "bcftools view".  Still in experimental stages.
* calc_ibd: Calculate IBD from VCF or PLINK file
  * Calculates a pairwise IBD estimate from either VCF or PLINK files using PLINK 1.9.
* call_bqsr: Base Quality Score Recalibration
  * Call [GATK BaseRecalibrator](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php) and return the tables for use in HaplotypeCaller
* call_genotypes: 
  * Obsolete, do not use; use geno_p instead.  Calls GATK GenotypeGVCFs.
* call_hc: 
  * Call [GATK HaplotypeCaller](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) and return gVCF files
* call_vqsr:
  * Calls [GATK VariantRecalibrator](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php) and returns the files needed to apply the recalibration
* cat_variants: combine_variants
  * Combines non-overlapping VCF files with the same subjects.  A reimplementation of GATK CatVariants (GATK CatVariants available upon request)
* combine_variants: combine_variants
  * Calls [GATK CombineVariants](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineVariants.php) to merge VCF files
* gen_ancestry: 
  * Determine Ancestry from PCA.  Uses an eigenvector file and training dataset listing known ancestries.  Runs QDA to determine posterior ancestries for all samples, even those in the training set.
* gen_related_todrop: 
  * Uses a PLINK IBD file to determine the minimal set of samples to drop in order to generate an unrelated sample set.  Uses a minimum vertex cut algorithm of the related samples to get 
* geno_p: 
  - Calls [GATK GenotypeGVCFs](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php) in parallel by chromosome
* merge_gvcfs:
  * Calls [GATK CombineGVCFs](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineGVCFs.php)
* plink_merge: 
  * Merge PLINK bed/bim/fam files using PLINK 1.9
* select_variants: VCF QC
  * Calls [GATK SelectVariants](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php)
* variant_annotator: VCF QC
  * Calls [GATK VariantAnnotator](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_VariantAnnotator.php)
* vcf_annotate: Annotate VCF File
  * Use a variety of tools to annotate a sites-only VCF.  
* vcf_concordance: VCF Concordance
  * Generate concordance metrics from VCF file(s) using [GATK GenotypeConcordance](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeConcordance.php).  Not recommended for large files.
* vcf_gen_lof: 
  * Subset a VCF from vcf_annotate based on the given annotations to get a sites-only VCF of loss-of-function variants.
* vcf_pca: 
  * Uses PLINK 1.9 and eigenstrat 6.0 to calculate principal components from VCF or PLINK bed/bim/fam files.
* vcf_qc:
  * Calls [GATK ApplyRecalibration](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php) and [GATK VariantFiltration](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php) to apply filters to VCF files.
* vcf_query:
  * Calls "bcftools query" to extract annotations from the VCF file
* vcf_sitesonly: VCF QC
  * Generates a sites-only file from full VCF files.
* vcf_slice: Slice VCF File(s)
  * Return a small section of a VCF file (similar to tabix).  For large output, many small regions, or subsetting samples, use subset_vcf instead.
* vcf_summary: VCF Summary Statistics
  * Generate summary statistics for a VCF file (by sample and by variant)
* vcf_to_plink: 
  * Uses PLINK 1.9 to convert VCF files to PLINK bed/bim/fam files
