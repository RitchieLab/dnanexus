#!/bin/bash

set -e -x -o pipefail

main() {
	# switch to a temp directory and download all user input files
	
	TMPDIR="$(mktemp -d)"
	cd "$TMPDIR"
	mkdir input
	VCF_ARG="input/input.vcf"
	dx download "$vcf_file" -o "$VCF_ARG"
	CAT_CMD="cat"
	if [[ "$(dx describe "$vcf_file" --name)" == *".gz" ]]; then
		CAT_CMD="zcat"
	fi
	POSITION_ARG="."
	if [[ -n "$position_file" ]]; then
		POSITION_ARG="input/input.positions"
		dx download "$position_file" -o "$POSITION_ARG"
	fi
	REGION_ARG="."
	if [[ -n "$region_file" ]]; then
		REGION_ARG="input/input.regions"
		dx download "$region_file" -o "$REGION_ARG"
	fi
	GENE_ARG="."
	if [[ -n "$gene_file" ]]; then
		GENE_ARG="input/input.genes"
		dx download "$gene_file" -o "$GENE_ARG"
	fi
	
	
	# parse other user options
	
	ANNOS_ARG=""
	if [[ ${#annotations[*]} -gt 0 ]]; then
		ANNOS_ARG="$(IFS=" " ; echo "${annotations[*]}")"
	fi
	EXTRA_ARG=""
	if [[ "$extra" == "true" ]]; then
		EXTRA_ARG="extra"
	fi
	
	
	# extract and format annotations
	
	mkdir output
	$CAT_CMD "$VCF_ARG" | \
	vcf_anno_format.py \
		"output/$prefix" \
		$POSITION_ARG \
		$REGION_ARG \
		$GENE_ARG \
		$ANNOS_ARG \
		$EXTRA_ARG
	
	
	# upload outputs
	
	for f in output/*.txt ; do
		anno_file=$(dx upload "$f" --brief)
		dx-jobutil-add-output annotation_files --class="array:file" "$anno_file"
        done
	
}
