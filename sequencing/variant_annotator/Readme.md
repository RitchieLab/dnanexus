# GATK3's VariantAnnotator 


## What does this app do?

This app run's the [GATK VariantAnnotator](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_VariantAnnotator.php).
This tool is meant to annotate the variant calls in a given `VCF` file based on their context.  


## What data are required for this app to run?

This app requires:

- An array of `vcf` files to be annotated.
- A matching array of `tbi` index files for those `vcfs`.
- The GATK3 `jar` file.
- The reference genome the variants were called against in `fasta` format.
- A boolean value indicating whether the app should produce a sites-only `vcf` (no sample level information).
- A boolean value indicating whether the app should produce a full `vcf` (with sample level information).

## What does this app output?

This app has two sets of potential outputs:

- An annotated sites-only `vcf` and the corresponding `tbi` index.
- An annotated `vcf` complete with sample information and the corresponding `tbi` index.