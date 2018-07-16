<!-- dx-header -->
# Get Sites-Only VCF

## What does this app do?

This app generates a 'sites-only' file from VCF file(s) (`*.vcf.gz`).

## What data are required for this app to run?

This app requires VCF file(s) (`*.vcf.gz`) as an input. 

## What does this app output?

This app outputs a site-only VCF file(s) (`*.vcf.gz`) along with its index file. The output VCF file (`*.vcf.gz`) will contain the first eight columns of the original VCF file.

## How does this app work?

This app uses the cut command to get a sites-only subset (the first 8 columns) of the VCF file (`*.vcf.gz`). 