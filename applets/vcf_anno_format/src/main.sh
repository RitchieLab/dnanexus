#!/bin/bash

set -e -x -o pipefail

main() {
	# switch to a temp directory and download all user input files
	
	TMPDIR="$(mktemp -d)"
	cd "$TMPDIR"
	mkdir input
	VCF_FILE="input/input.vcf.gz"
	dx download "$vcf_file" -o "$VCF_FILE"
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
	
	echo "annotations=$annotations"
	echo "extra=$extra"
	
}
