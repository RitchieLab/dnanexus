# dnanexus
Dnanexus Apps and Scripts

## applets 
* binning_step0: BioBin Pipeline
 * biobin_pipeline
* binning_step1: BioBin Pipeline
 * biobin_pipeline
* binning_step2: BioBin Pipeline
 * biobin_pipeline
* binning_step3: BioBin Pipeline
 * biobin_pipeline
* impute2_group_join: Impute2_group_join
 * This app can be used to merge multiple imputed impute2 files
* plato_biobin: PLATO BioBin Regression Analysis
 * PLATO_BioBin
* vcf_batch: VCF Batch effect tester 
 * vcf_batch

## apps
* association_result_annotation: Annotate GWAS, PheWAS Assocaitions
 * association_result_annotation
* biobin: 
 * This app runs the latest development build of the rare variant binning tool BioBin.
* generate_phenotype_matrix: Generate Phenotype Matrix
 * generate_phenotype_matrix
* genotype_case_control: Generate Case/Control by Genotype
 * App provides case and control number by each genotype
* impute2: imputation
 * This will perfrom imputation using Impute2
* impute2_to_plink: Impute2 To PLINK
 * Convert Impute2 file to PLINK files
* plato_single_variant: PLATO - Single Variant Analysis 
 * Apps allows you to run single variant association testing against single phenotype (GWAS) or multiple phenotype (PheWAS) test
* rl_sleeper_app: sleeper
 * This App provides some useful tools when working with data in DNANexus. This App is designed to be run on the command line with "dx run --ssh RL_Sleeper_App" in the project that you have data that you want to explore (use "dx select" to switch projects as needed).
* shapeit2: SHAPEIT2
 * This app do phasing using SHAPEIT2
* strand_align: Strand Align
 * Strand Align prior to phasing
* vcf_annotation_formatter: 
 * Extracts and reformats VCF annotations (CLINVAR, dbNSFP, SIFT, SNPEff)
* **QC_apps subfolder:**
  * drop_marker_sample: Drop Markers and/or Samples (PLINK)
    * drop_marker_sample
 * drop_relateds: Relatedness Filter (IBD)
   * drop_relateds
 * extract_marker_sample: Drop Markers and/or Samples (PLINK)"
   * extract_marker_sample
 * maf_filter: Marker MAF Rate Filter (PLINK)
   * maf_filter
 * marker_call_filter: Marker Call Rate Filter (PLINK)
   * marker_call_filter
 * missing_summary: Missingness Summary (PLINK)
   * Returns missingness rate by sample
 * pca: Principal Component Analysis using SMARTPCA 
    * pca
 * sample_call_filter: Sample Call Rate Filter (PLINK)
   * sample_call_filter

## scripts
* cat_vcf.py
  * 
* download_intervals.py
  * 
* download_part.py
  * 
* estimate_size.py
  * 
* interval_pad.py
  * This reads a bed file from standard input, pads the intervals, sorts and then outputs the intervals guranteed to be non-overlapping
* update_applet.sh
  * 

## sequencing
* bcftools_view: 
  * Calls "bcftools view".  Still in experimental stages.
* calc_ibd: 
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
  * Calls [GATK GenotypeGVCFs](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php) in parallel by chromosome
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
  * Calls "bcftools query" to extract annotations from the VCF file.  Used in the stripping of files for MEGAbase
* vcf_sitesonly: VCF QC
  * Generates a sites-only file from full VCF files.
* vcf_slice: Slice VCF File(s)
  * Return a small section of a VCF file (similar to tabix).  For large output, many small regions, or subsetting samples, use subset_vcf instead.
* vcf_summary: VCF Summary Statistics
  * Generate summary statistics for a VCF file (by sample and by variant)
* vcf_to_plink: 
  * Uses PLINK 1.9 to convert VCF files to PLINK bed/bim/fam files
