 #!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

main() {

	echo "Value of variants_vcfgzs: '$variants_vcfgzs'"

	VCF_LIST=$(mktemp)
	# split into sub jobs
	for i in "${!variants_vcfgzs[@]}"; do	
		# input : single VCF.gz file 
		SUBJOB=$(dx-jobutil-new-job get_so -ivariants_vcfgz="${variants_vcfgzs[$i]}")
		# pass output 
			dx-jobutil-add-output vcf_out --array "$SUBJOB:vcf_out" --class=jobref
		dx-jobutil-add-output vcfidx_out --array "$SUBJOB:vcfidx_out" --class=jobref	
	done
}

# get sites only 
get_so() {

	echo "Value of variants_vcfgz: '$variants_vcfgz'"

	OUTDIR=$(mktemp -d)

	PREFIX=$(dx describe --name "$variants_vcfgz" | sed 's/\.vcf.\(gz\)*$//')

	# get the first 8 columns of the VCF file 
	dx cat "$variants_vcfgz" | pigz -dc | cut -f 1-8 | bgzip -c > $OUTDIR/header.$PREFIX.vcf.gz	
	tabix -p vcf $OUTDIR/header.$PREFIX.vcf.gz

	# upload the output VCF along with this index
	vcf_hdr_out=$(dx upload "$OUTDIR/header.$PREFIX.vcf.gz" --brief)
	vcfidx_hdr_out=$(dx upload "$OUTDIR/header.$PREFIX.vcf.gz.tbi" --brief)
	dx-jobutil-add-output vcf_out "$vcf_hdr_out" --class=file
	dx-jobutil-add-output vcfidx_out "$vcfidx_hdr_out" --class=file
	
}
