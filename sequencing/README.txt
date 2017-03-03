Last revised: 03/02/2017

Please note: Due to a problem that arose during merging of gvcf files, 1 GHS samples (out of 61,019 in the F60K) did not get joint called with the rest of the samples and thus are not present in the F50K jVCF files as well as downstream files (e.g. plink files, PCA etc). The 6 missing samples are listed below and will be included in joint calling during the next data freeze.
GHS_PT984567_182976186

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
|   Contains the VCF LOF rollups at various levels.  The LOF rollups are
|   a subset of variants that meet a given criteria.  See the "VCF Annotation and
|   LOF Rollup" section for details on how these files are created and how to use
|   them effectively.
+-PCA
| | Contains the principal component eigenvectors and eigenvalues as returned 
| | by eigenstrat, as well as the results of genetically informed ancestry determination.
| | See the "IBD and PCA" section of this document for details on how the PCs were
| | generated, as well as "Genetically informed ancestry determination" section for more 
| | details on how the ancestry file was generated.
| +-GHS_only/onTarget/masked/projected
| |   Same as above but pr
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
| | Contains the VCF files as output by various GATK utilities.  Within these
| | folders are the VCFs intended for use in analyses.  Note that most folders
| | have a "by_chrom" subfolder which contains the chromosomal-level VCFs for
| | use in analyses and apps that can be parallelized, as well as "GHS_only" 
| | subfolder which contains GHS samples only (with 8 positive control samples 
| | excluded).  All VCF files within these directories should contain exactly 
| | the same number of variants; the differences lie in the annotations and 
| | filtration applied to each variant.
| +-filtered_vcf
| |   Contains the VCF after applying both the Variant Quality Score Recalibration
| |   filter as well as a genotype-level quality filter.  See the "VCF QC" section
| |   of this document for details on Quality Control of VCF files.  All data from
| |   the raw VCF is contained; only annotations of substandard quality have been
| |   added.  This is the VCF recommended for most analyses.  Note that the 
| |   AC, AF, and AN annotations in the INFO field do not take the genotype-level
| |   filters into account.  For those annotations calculated correctly, please 
| |   use the "masked" versions of the VCF files.
| +-masked_vcf
| |   Contains the VCF after applying both VQSR and genotype-level quality filters
| |   as above.  Variant calls failing the genotype-level filters have been
| |   manually set to no-call to accommodate tools that cannot read the "FT" 
| |   genotype annotation.  Additionally, the AC, AF, and AN annotations in the
| |   INFO field are updated to accommodate the genotype-level filters used.
| |   Variant level filters are not removed from this file, but the annotations
| |   still exist in the FILTER column of the VCF.
| +-raw_vcf
| |   Contains the VCF as output by GATK's GenotypeGVCFs tool.  This is not 
| |   intended for routine analysis, but is provided as a means for a user to 
| |   perform his/her own QC of the VCF file if the QC procedures defined here 
| |   are deemed inadequate.  Note that only chromosomal-level VCFs are provided
| |   for the raw VCF files.
| +-vcf_headers
| | | Contains sites-only VCF files (first 8 columns) corresponding to the 
| | | various VCF files mentioned above.  Sites-only files can be useful for 
| | | obtaining variant-level information about VCF files without needing to
| | | download the large VCF files that contain all genetic information.
| | +-annotated
| | |   Contains the masked sites-only VCF file with added annotations used 
| | |   in determining LOF rollups.  See the "VCF Annotation and LOF Rollup" 
| | |   section of this document for details on the annotations available and
| | |   how they are used to generate an LOF rollup file.
| | +-filtered
| | |   The sites-only VCF corresponding to the "filtered_vcf" folder above.
| | +-masked
| | |   The sites-only VCF corresponding to the "masked_vcf" folder above.
| | +-raw
| |     The sites-only VCF corresponding to the "raw_vcf" folder above.
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

The PCA analysis results are provided as an eigenstrat (6.1.4) eigenvector file
using fast eigenvector approximation projected onto the 1000 genomes data 
(lifted over b37).  For PCA, we used only SNPs with a MAF >1%
and a Hardy-Weinberg equilibrium p-value less than 0.000001.  We then performed 
LD pruning exactly as above for IBD.  This yielded a final set of 32,618 markers 
used in calculating the principal components.  We did not perform any outlier 
removal using eigenstrat, so PCs are given for all samples.  For a more thorough
PCA calculation, a user may wish to recalculate the principal components based 
only on the samples included in a given analysis. We also performed PCA excluding
the 8 positive control samples. The same analysis yielded a final set of 32,613
markers (GHS_only subfolder).

VCF Annotation and LOF Rollup
=============================

The masked VCFs were annotated in order to determine if each variant
is potentially detrimental and should receive increased scrutiny from researchers.
The LOF roll ups are based off of what types of requests were made on the 50k release 
by the Return of Results project and request by other researchers. As Regeneron proved 5
levels of rollups we have as well.


- M1: VEP: Moderate Impact (Highest annotation from Ensembl or RefSeq), MAF<5%, present in CLINVAR, PRESENT in HGMD
- M2: VEP: HIGH, (Highest annotation from Ensembl or RefSeq) MAF<1%, CLINVAR 1*+ Not Begnin, HGMD (DFP, DP, DM? or DM)
- M3: VEP: HIGH, (Cannonical RefSeq transcripts) MAF<1%, CLINVAR 2*+ (Pathogenic, Likely Pathogenic), HGMD (DM: Disease Causing)
- M4: All Varaints MAF<1%
- M5: All Varaints MAF<.1%

Within the annotated VCF, we provide the following variant-level annotations
that are then used to determine if a particular variant is loss of function, 
nonsynonymous, or deleterious.  The annotations provided are:
 - VEP annotation using Ensembl Tools Release-87, for GRCh 38, anntoated with 
	Ensembl and Refseq annotations, utilizing the merged data set and the everything flag
 - dbNSFP v3.3
 - CLINVAR, Jan 2017 release
 - HGMD 2016.4 b38 release

If multiple annotations exist for a given variant (as would be in the case of 
multiallelic variants or variant annotations specific to a transcript), a variant
is considered to be in a given class if ANY annotation meets the specification
above.  As an example, if a variant had possible alleles of C/T/G, and was 
annotated in SNPEff as "C|exon_loss_variant,G|synonymous_variant", this variant
would be considered to meet the definition of loss of function, even though only
one possible allele confers lost function.

Note that above, M5 is a subset of M4, and M2 is a subset of M1 and M4, M3 is subset of M2.

