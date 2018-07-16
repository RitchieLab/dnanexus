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
* By default, all annotations are applied. To disable any of the annotations, set its corresponding environment variable to 'false'.
* $annotate_header is an environment variable that, when $VEP=true, determines whether to annotate using the VEP executable and database (true), or to use a pre-computed VCF file with all the annotations for each position stored in the INFO/CSQ field.

# Output

* VCF and index lists the same length, and in the same order, as the input lists.
-->
