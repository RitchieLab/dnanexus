# VCF File Annotation (DNAnexus Platform App)

# What does this app do?
This app annotates a VCF File with the most recent ClinVar TSV file (downloaded from `ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz`) and with a user-provided HGMD PRO vcf file. Because of VCF Normalization differences, some annotations will be skipped. The app does not function with UCSC chromosome labeling. As a result, the input VCF file(s) cannot include `chr` prefixes.

# What are the inputs to this app?
This app takes as input an array of VCF (`*.vcf.gz`) file(s), an optional array of VCF file indices (`*.vcf.gz.tbi`), and a TSV (`*.tsv`) file containing the clinical variants mentioned above. Optionally, a VCF (`*.vcf.gz`) file from the Human Gene Mutation Database (HGMD) may also be provided.
# What are the outputs to this app?
This app outputs an array of annotated VCF file(s) and their associated index file(s). 

