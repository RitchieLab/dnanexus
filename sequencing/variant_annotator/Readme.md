# GATK3 VariantAnnotator 


## What does this app do?

This app runs the [GATK VariantAnnotator](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_VariantAnnotator.php).
This tool is meant to annotate the variant calls in a given `VCF` file based on their context.  
Given an array of `VCF` files, this app will launch a separate subjob for each `VCF` and run GATK's VariantAnnotator.  Users can request that a sites-only `VCF` (no sample level information), a full `VCF`, or both be produced.


## What data are required for this app to run?

This app requires:

- An array of `VCF` files to be annotated (`*.vcf.gz`).
- A matching array of index files for those `VCFs` (`*.vcf.gz.tbi`).
- The GATK3 `*.jar` file.
- The reference genome the variants were called against in `FASTA` format.
- A boolean value indicating whether the app should produce a sites-only `VCF` (no sample level information).
- A boolean value indicating whether the app should produce a full `VCF` (with sample level information).

## What does this app output?

This app has two sets of potential outputs:

- An annotated sites-only `VCF` and the corresponding `TBI` index.
- An annotated `VCF` complete with sample information and the corresponding `TBI` index.
