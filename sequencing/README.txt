Last revised: 02/10/2016

Please note: Due to a problem that arose during merging of gvcf files, 6 GHS samples (out of 51,289 in the F50K) did not get joint called with the rest of the samples and thus are not present in the F50K jVCF files as well as downstream files (e.g. plink files, PCA etc). The 6 missing samples are listed below and will be included in joint calling during the next data freeze.
GHS_PT116761_180731223
GHS_PT105584_180731343
GHS_PT735774_180731353
GHS_PT597988_180731397
GHS_PT1154675_180731407
GHS_PT997609_180731424

Data Structures
===============

The following describes the folders contained within the recalling project.
Any folder not listed in this directory structure will contain either temporary
data or old data from a failed step in the process and can be safely ignored.

+-bqsr
|   Contains the tables needed for Base Quality Score Recalibration.
|   Not intended for use in analysis.
+-gvcf
|   Contains the individual genomic VCF files as output by GATK's HaplotypeCaller
|   Also contains the merged gvcf files (both tiers)
|   Not intended for use in analysis.
+-IBD
|   Contains the relatedness estimate obtained from PLINK.  See the "IBD and PCA"
|   section of this document for details on how this estimate was obtained.
+-LOF
|   Contains the sites-only LOF rollups at various levels.  The LOF rollups are
|   a subset of variants that meet a given criteria.  See the "VCF Annotation and
|   LOF Rollup" section for details on how these files are created and how to use
|   them effectively.
+-PCA
| | Contains the principal component eigenvectors and eigenvalues as returned 
| | by eigenstrat, as well as the results of genetically informed ancestry determination.
| | See the "IBD and PCA" section of this document for details on how the PCs were
| | generated, as well as "Genetically informed ancestry determination" section for more 
| | details on how the ancestry file was generated.
| +-GHS_only
| |   Same as above but for GHS samples only (8 positive control samples were excluded) 
+-PLINK
| | Contains PLINK files obtained from VCF to PLINK conversion as output by 
| | PLINK 1.9's "--make-bed" argument.  For all PLINK files, the single ID
| | given in the VCF file is both the FID and IID in the PLINK files.  For each 
| | subdirectory mentioned, PLINK files are provided both by chromosome in a
| | "by_chrom" subfolder as well as concatenated in the top level folder.
| | The PLINK arguments in common to all of the following subfolders is:
| | --double-id --id-delim "' '" --set-missing-var-ids @:#:\$1 --vcf-filter -allow-no-sex
| +-all
| |   Contains all (passed) variants in the PLINK file.  Each variant site
| |   includes the referent allele and the most common allele, with less frequent
| |   alleles at multiallelic sites set to missing.
| +-biallelic
| |   As above, but only includes variants that are truly biallelic (after 
| |   applying any genotype level filtration).  The additional argument to  PLINK is:
| |   --biallelic-only
| +-snps
| |   As above, but only containing variants that are biallelic SNPs (after 
| |   filtration).  The additional PLINK arguments were:
| |   --biallelic-only --snps-only
| +-GHS_only
| | +-all
| | +-biallelic
| | +-snps
| |   Same as above but for GHS samples only (8 positive control samples were excluded) 
+-vcf
  | Contains the VCF files as output by various GATK utilities.  Within these
  | folders are the VCFs intended for use in analyses.  Note that most folders
  | have a "by_chrom" subfolder which contains the chromosomal-level VCFs for
  | use in analyses and apps that can be parallelized, as well as "GHS_only" 
  | subfolder which contains GHS samples only (with 8 positive control samples 
  | excluded).  All VCF files within these directories should contain exactly 
  | the same number of variants; the differences lie in the annotations and 
  | filtration applied to each variant.
  +-filtered_vcf
  |   Contains the VCF after applying both the Variant Quality Score Recalibration
  |   filter as well as a genotype-level quality filter.  See the "VCF QC" section
  |   of this document for details on Quality Control of VCF files.  All data from
  |   the raw VCF is contained; only annotations of substandard quality have been
  |   added.  This is the VCF recommended for most analyses.  Note that the 
  |   AC, AF, and AN annotations in the INFO field do not take the genotype-level
  |   filters into account.  For those annotations calculated correctly, please 
  |   use the "masked" versions of the VCF files.
  +-masked_vcf
  |   Contains the VCF after applying both VQSR and genotype-level quality filters
  |   as above.  Variant calls failing the genotype-level filters have been
  |   manually set to no-call to accommodate tools that cannot read the "FT" 
  |   genotype annotation.  Additionally, the AC, AF, and AN annotations in the
  |   INFO field are updated to accommodate the genotype-level filters used.
  |   Variant level filters are not removed from this file, but the annotations
  |   still exist in the FILTER column of the VCF.
  +-raw_vcf
  |   Contains the VCF as output by GATK's GenotypeGVCFs tool.  This is not 
  |   intended for routine analysis, but is provided as a means for a user to 
  |   perform his/her own QC of the VCF file if the QC procedures defined here 
  |   are deemed inadequate.  Note that only chromosomal-level VCFs are provided
  |   for the raw VCF files.
  +-vcf_headers
  | | Contains sites-only VCF files (first 8 columns) corresponding to the 
  | | various VCF files mentioned above.  Sites-only files can be useful for 
  | | obtaining variant-level information about VCF files without needing to
  | | download the large VCF files that contain all genetic information.
  | +-annotated
  | |   Contains the masked sites-only VCF file with added annotations used 
  | |   in determining LOF rollups.  See the "VCF Annotation and LOF Rollup" 
  | |   section of this document for details on the annotations available and
  | |   how they are used to generate an LOF rollup file.
  | +-filtered
  | |   The sites-only VCF corresponding to the "filtered_vcf" folder above.
  | +-masked
  | |   The sites-only VCF corresponding to the "masked_vcf" folder above.
  | +-raw
  |     The sites-only VCF corresponding to the "raw_vcf" folder above.
  +-vqsr
     Contains the supporting files necessary to run Variant Quality Score
     Recalibration, both SNP and INDEL.  Not intended for general analysis
     use, but can be informative in assessing the effectiveness of the 
     VQSR filtration step.
     
VCF QC
======

VCF Files as produces by GATK's HaplotypeCaller (and by extension GenotypeGVCFs)
are very sensitive to perceived variation and will call many sequencing errors
as true variation, yielding a high false positive rate without further QC.
Using the GATK best practices as a guide, we apply two main methods of filtration:
 - Variant Quality Score Recalibration
 - Genotype Quality

For VQSR, we followed https://www.broadinstitute.org/gatk/guide/article?id=1259
as a base, using the describe resources and the following annotations to inform
the VQSR model:
 - QD: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_QualByDepth.php
 - FS: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_FisherStrand.php
 - SOR: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
 - MQ (SNPs only): https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_RMSMappingQuality.php
 - MQRankSum: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityRankSumTest.php
 - ReadPosRankSum: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_ReadPosRankSumTest.php
 - InbreedingCoeff: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_InbreedingCoeff.php

Additionally, we use a maximum of 4 Gaussians when fitting the INDEL model.
 
For both SNP and INDEL VQSR, we provide the following tranche levels:
 - 99.9%
 - 99.5%
 - 99.0%
 - 95.0%
 - 90.0%
 
When filtering VQSR variants, we used the following truth sensitivity thresholds:
 -SNPS:   99.5%
 -INDELS: 99.0%
 
Following the guidance in http://www.biomedcentral.com/1471-2105/15/125, we also
use a genotype-level filter requiring the Genomic Quality (GQ) to be at least 20
for all high-quality genotype calls.  This corresponds to a 99% confidence in the
genotype call for every individual call.

Additionally, to asses the quality of the QC steps above, the VCF files are 
jointly called with the positive control samples NA12877, NA12878, NA12880,
NA12882, NA12889, NA12890, NA12891, and NA12892. These samples were excluded in 
certain downstream files (GHS_only subfolders).

IBD and PCA
===========

For convenience, we provide IBD estimation from PLINK and a principal component
analysis from eigenstrat.  For both analyses, we used only the biallelic SNPs on
the autosomes.

The IBD analysis results are provided as a gzipped IBS file as produced by 
PLINK 1.9's "--genome gz" command (see https://www.cog-genomics.org/plink2/ibd).
For IBD, we used only SNPs with a MAF >5% and Hary-Wienberg equillibrium p-value 
less than 0.000001.  We then further LD-pruned the dataset using the 
"--indep-pairwise 50 5 0.5" argument to PLINK.  This corresponds to eliminating 
one marker in every pair that has an LD R^2 value of 0.8, examining markers 
within a 50-marker window, sliding by 5 markers at a time.  This yielded a final 
set of 18,992 markers used in calculating IBD.  Note that all pairwise 
individuals are provided in the genome report; it is up to the user to specify 
an appropriate relatedness threshold for his/her analysis.

The PCA analysis results are provided as an eigenstrat (6.0) eigenvector file
using fast eigenvector approximation projected onto the 1000 genomes data 
(5/8/2013 release).  For PCA, we used only SNPs with a MAF >1%
and a Hardy-Weinberg equilibrium p-value less than 0.000001.  We then performed 
LD pruning exactly as above for IBD.  This yielded a final set of 32,618 markers 
used in calculating the principal components.  We did not perform any outlier 
removal using eigenstrat, so PCs are given for all samples.  For a more thorough
PCA calculation, a user may wish to recalculate the principal components based 
only on the samples included in a given analysis. We also performed PCA excluding
the 8 positive control samples. The same analysis yielded a final set of 32,613
markers (GHS_only subfolder).

Genetically informed ancestry determination
===========================================

To calculate the ancestry for each person in the population, we used the top 2 eigenvectors (configurable) for each sample included in the Freeze 40K race file to generate a training set for a quadratic discriminant analysis (excluding any samples which were labeled as "Other").  We then calculated the posteriori probabilities for every sample in the joint-called Freeze 50K data using the top 2 eigenvectors and assigned the race with the highest probability, provided that the probability was over 80% (configurable).  Samples not passing the given probability threshold were labeled as "Other".  For convenience, we provide both the race classification and the posteriori probability of the most common class in columns 2 and 3, respectively. Please note, the ancestry file contains the "1000 genomes" samples in addition to the 51,283 GHS samples.

VCF Annotation and LOF Rollup
=============================

The masked sites-only VCF was annotated in order to determine if each variant
is potentially detrimental and should receive increased scrutiny from researchers.
We have designed the loss of function rollups to be comparable to the LOF rollups
provided alongside the pVCF files released by Regeneron.  From the RGN readme 
files, there are 5 levels of LOF rollups:

- M1: LOF variants only
- M2: LOF and NS variants (MAF<1%)
- M3: LOF and predicted deleterious NS variants (5 algorithms; no MAF cutoff)
- M4: LOF and predicted deleterious NS variants (1 algorithm; MAF<1%)
- M5: LOF variants only (MAF<1%)

Within the annotated VCF, we provide the following variant-level annotations
that are then used to determine if a particular variant is loss of function, 
nonsynonymous, or deleterious.  The annotations provided are:
 - SNPEff 4.1l using GRCh 37 and Ensembl 75 databases.
 - dbNSFP 2.9
 - SIFT 5.2.2
 - CLINVAR, 10/2015 release

For each annotation, each variant is annotated in a style similar to the SNPEff
"ANN" field, see http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf for
details.  Specifically, each alternate allele is annotated and those annotations
are combined for each line.  The annotated VCF file will match all of the other
VCF files line for line.

We used the following definitions of the terms used in the RGN readme above:
- Loss Of Function Variant
  A variant that has one of the following SNPEff roles:
    - chromosome_number_variation
    - exon_loss_variant
    - frameshift_variant
    - stop_gained
    - stop_lost
    - start_lost
    - splice_acceptor_variant
    - splice_donor_variant
    - rare_amino_acid_variant
    - transcript_ablation
    - disruptive_inframe_insertion
    - disruptive_inframe_deletion
  Note that this should consist of all SNPEff roles with a HIGH impact modifier,
  plus the addition of the disruptive insertion/deletion.
- Nonsynonymous Variant
  A variant that is identified as a Loss of Function variant above, or has one 
  of the following SNPEff roles:
    - missense_variant
    - inframe_insertion
    - inframe_deletion
    - 5_prime_UTR_truncation
    - 3_prime_UTR_truncation
    - splice_region_variant
    - splice_branch_variant
    - coding_sequence_variant
    - regulatory_region_ablation
    - TFBS_ablation
    - 5_prime_UTR_premature_start_codon_gain_variant
    - non_canonical_start_codon
- Predicted Deleterious
  A variant is predicted deleterious if a consensus of M (user-specified) 
  algorithms determine that a variant is deleterious.  The algorithms used to
  determine deleteriousness are provided by dbNSFP; this means that predicted
  deleteriousness is limited to SNPs only.  The algorithms available for 
  prediction of effect are as follows (an asterisk indicates that this algorithm
  is used in generating the M1-M5 VCF files released):
    * SIFT
    * PolyPhen2 (HDIV training set)
    * PolyPhen2 (HVAR training set)
    * LRT
    * MutationTaster
    - MutationAssessor
    - FATHMM
    - MetaSVM
    - MetaLR
    - PROVEAN

If multiple annotations exist for a given variant (as would be in the case of 
multiallelic variants or variant annotations specific to a transcript), a variant
is considered to be in a given class if ANY annotation meets the specification
above.  As an example, if a variant had possible alleles of C/T/G, and was 
annotated in SNPEff as "C|exon_loss_variant,G|synonymous_variant", this variant
would be considered to meet the definition of loss of function, even though only
one possible allele confers lost function.

Given the above definition, a variant is in each of the following classes if the
following conditions are met:

M1: Loss of Function Variant
M2: Loss of Function Variant OR Nonsynonymous Variant with MAF < 1%
M3: Loss of Function Variant OR Predicted Deleterious with all 5 algorithms
M4: Loss of Function Variant OR Predicted Deleterious with 1 algorithm and MAF <1%
M5: Loss of Function Variant with MAF<1%

Note that above, M5 is a subset of M1, and M1 is a subset of M2, M3, and M4.  
M2, M3, and M4 are neither subsets nor supersets of each other.

The LOF rollups are provided only as sites-only VCFs.  If a user wishes to 
generate a more specific loss of function set that is not provided, the app 
"Generate LOF Rollup" will take the annotated sites-only VCF and filter the 
results according to the SNPEff, dbNSFP, or CLINVAR annotations.  In order to 
obtain the full VCF corresponding to the appropriate LOF mask, a user can use
the "bcftools isec" command, or the GATK SelectVariants tool.  If using the 
latter, the app "Subset VCF" (named "select_variants") provides a mechanism
to generate the full VCF entirely in DNANexus.
