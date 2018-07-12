<!-- dx-header -->
# VCF File Annotation (DNAnexus Platform App)

# What does this app do?
It annotates a VCF File with most recent ClinVar TSV file ( downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz ) and with a user provided HGMD PRO vcf file. Because of VCF Normalization difference some annotations will be skipped. It does not function with UCSC chromosome labeling. As a result, the input VCF file(s) will need to have their chr prefix removed.

# What are the inputs to this app?
This app takes as input an array of VCF(`*.vcf.gz`) file(s), an optional array of indexed VCF file(s) (`*.vcf.gz.tbi`), and a TSV (`*.tsv`) containing clinical variants mentioned aboved. Additionally, a VCF file from the Human Gene Mutation Database (HGMD) may be optionally provided. 

# What are the outputs to this app?
This app outputs an array of annotated VCF file(s) and their associated index file(s). 

