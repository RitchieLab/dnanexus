#!/bin/bash
# vcf_qc 0.0.1

set -e -x -o pipefail


main() {

	echo "Value of query_str: '$query_str'"
    echo "Value of variants_vcfgzs: '$variants_vcfgzs'"
    echo "Value of variants_vcfgztbis: '$variants_vcfgztbis'"
    echo "value of DX_ASSETS_ID: $DX_ASSETS_ID"

	# Sanity check - make sure vcf + vcfidx have same # of elements
	if test "${#variants_vcfgztbis[@]}" -ne "${#variants_vcfgzs[@]}"; then
		dx-jobutil-report-error "ERROR: Number of VCFs and VCF indexes do NOT match!"
	fi

	# first, we need to match up the VCF and tabix index files
	# To do that, we'll get files of filename -> dxfile ID
	VCF_LIST=$(mktemp)
	for i in "${!variants_vcfgzs[@]}"; do	
		dx describe --json "${variants_vcfgzs[$i]}" | jq -r ".name,.id" | tr '\n' '\t' | sed 's/\t$/\n/' >> $VCF_LIST
	done
	
	VCFIDX_LIST=$(mktemp)
	for i in "${!variants_vcfgztbis[@]}"; do	
		dx describe --json "${variants_vcfgztbis[$i]}" | jq -r ".name,.id" | tr '\n' '\t' | sed -e 's/\t$/\n/' -e 's/\.tbi\t/\t/' >> $VCFIDX_LIST
	done
	
	# Now, get the prefix (strip off any .tbi) and join them
	JOINT_LIST=$(mktemp)
	join -t$'\t' -j1 <(sort -k1,1 $VCF_LIST) <(sort -k1,1 $VCFIDX_LIST) > $JOINT_LIST
		
	# Ensure that we still have the same number of files; throw an error if not
	if test $(cat $JOINT_LIST | wc -l) -ne "${#variants_vcfgzs[@]}"; then
		dx-jobutil-report-error "ERROR: VCF files and indexes do not match!"
	fi
	
	# and loop through the file, submitting sub-jobs
	CAT_ARGS=""
	while read VCF_LINE; do
		VCF_DXFN=$(echo "$VCF_LINE" | cut -f2)
		VCFIDX_DXFN=$(echo "$VCF_LINE" | cut -f3)		
	
		SUBJOB=$(dx-jobutil-new-job run_query -ivariants_vcfgzs:file="$VCF_DXFN" -ivariants_vcfgztbis:file="$VCFIDX_DXFN" -iquery_str:string=$query_str)
		CAT_ARGS="$CAT_ARGS -iquery_in:array:file=$SUBJOB:query_gzs"
		
		if test "$cat_results" = "false"; then
			# for each subjob, add the output to our array
    		dx-jobutil-add-output query_gzs --array "$SUBJOB:query_gzs" --class=jobref
    	fi
		
	done < $JOINT_LIST
	
	if test "$cat_results" = "true"; then
		CATJOB=$(dx-jobutil-new-job cat_results $CAT_ARGS -iprefix:string=$prefix)
		dx-jobutil-add-output query_gzs --array "$CATJOB:query_gzs" --class=jobref
	fi
	
}

run_query() {

    echo "Value of variants_vcfgzs: '$variants_vcfgzs'"
    echo "Value of variants_vcfgztbis: '$variants_vcfgztbis'"
	echo "Value of query_str: '$query_str'"

    dx download "$variants_vcfgzs" -o input.vcf.gz
    dx download "$variants_vcfgztbis" -o input.vcf.gz.tbi
    
	OUT_DIR=$(mktemp -d)
	PREFIX=$(dx describe --name "$variants_vcfgzs" | sed 's/\.vcf.\(gz\)*$//')

	# TODO: parallelize smartly using -r chr:from-to and output accordingly
	bcftools query input.vcf.gz -f "$query_str"  -H | sed -e'1 s/\[[0-9]*\]//g' -e '1 s/  *//g' |  bgzip -c > $OUT_DIR/$PREFIX.query.gz
	query_gzs=$(dx upload "$OUT_DIR/$PREFIX.query.gz" --brief)

    dx-jobutil-add-output query_gzs "$query_gzs" --class=file
}

cat_results() {

	echo "Value of query_in: '$query_in'"
	echo "Value of prefix: '$prefix'"

	CAT_ARGS=""
	for i in "${!query_in[@]}"; do
		CAT_ARGS="$CAT_ARGS -V $(dx describe --json "${query_in[$i]}" | jq -r .id)"
	done
	
	# get the dict
	WKDIR=$(mktemp -d)
	cd $WKDIR
	
	OUTDIR=$(mktemp -d)
	
	cat_vcf.py -D /usr/bin/human_g1k_v37_decoy.dict $CAT_ARGS -o $OUTDIR/$prefix.query.gz
	query_gzs=$(dx upload $OUTDIR/$prefix.query.gz --brief)
	dx-jobutil-add-output query_gzs $query_gzs
}
	
