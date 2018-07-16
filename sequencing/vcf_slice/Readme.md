# Slice VCF File(s) (DNAnexus Platform App)

Return a small section of a VCF file (similar to tabix).  For large output, 
many small regions, or subsetting samples, use the app `select_variants` instead.

## App Inputs

This app requires 1 or more variant files in VCF format (`*.vcf.gz`) to be
filtered to the target region and a corresponding Tabix index file 
(`*.vcf.gz.tbi`) as well as a region file in UCSC list format (`*.list`)
which specifies the target regions to filter VCF files to.

The region file should appear in the following format: `chrom:start-end`

Example:
```
chr1:10273-104000
chr1:105000-109000
chr1:109000-112000
```

Or:
```
1:10273-104000
1:105000-109000
1:109000-112000
```

The chromosome names should correspond to the names found in the VCF files.

## App Output

The app returns VCF files in `*.vcf.gz` format, with corresponding tabix index 
files in `*.vcf.gz.tbi* format.
