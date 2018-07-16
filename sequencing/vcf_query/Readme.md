<!-- dx-header -->
# vcf_query (DNAnexus Platform App)

## What does this app do?
This app calls bcftools query to extract annotations from VCF files (`*.vcf.gz`).

## What data are required for this app to run?

This app requires VCF file(s) (`*.vcf.gz`) along with its respective index file(s) (`*.vcf.gz.tbi`) as input. Also, as a required input, the app takes in a string input consisting of the query to be performed on the VCF file. This can contain fields of the VCF file in a tab-delimited format. For example, the default query the app will perform is %CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%INFO/QD\\t%INFO/AN\\t%INFO/AC[\\t%GT]\\n.

Additonally, the app has an option to concatenate the results into a single file. 

## What does this app output?

The applet outputs the VCF file annotations as a (`*.query.gz`) file.