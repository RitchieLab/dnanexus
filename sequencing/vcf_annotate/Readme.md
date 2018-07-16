<!-- dx-header -->
# Annotate VCF File (DNAnexus Platform App)

Use a variety of tools to annotate a sites-only VCF.
<!-- /dx-header -->

<!--
Annotation is the process of adding additional information about each variant in a VCF file, derived from public databases.

This script annotates one or more VCF files with one or more of the following databases:

* [Variant Effect Predictor (VEP)](https://useast.ensembl.org/info/docs/tools/vep/index.html): a tool that aggregates multiple annotation databases. The annotations for each variant are all concatenated together in a single INFO field (CSQ by default).
* [Human Gene Mutation Database (HGMD)](http://www.hgmd.cf.ac.uk/ac/index.php): Database of known human coding variants.
* [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/): Database of known clinically relevant variants.
* [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP): Database of functional prediction scores for non-synonymous variants.

# Input

* List of one or more VCF files.
* Each VCF file must be indexed using tabix, which must be specified in a list the same length of, and in the same order as, the VCF list.
* VCF and TBI files for each of the annotation databases listed above that you want to apply.

## VEP

VEP annotations can be applied by installing and running the VEP executable instead of providing the VEP VCF. If you wish to use this option, download the VEP [source](https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html) and [caches](ftp://ftp.ensembl.org/pub/release-87/variation/VEP/) files.

# Output

* VCF and index lists the same length, and in the same order, as the input lists.
-->
