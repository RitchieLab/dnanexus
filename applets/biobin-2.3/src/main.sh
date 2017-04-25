#!/bin/bash

set -e -x -o pipefail

main() {
	
	ulimit -c unlimited
	
	
	# switch to a temp directory and download all user input files
	
	NUM_CORES="$(nproc --all)"
	TMPDIR="$(mktemp -d)"
	cd "$TMPDIR"
	mkdir input
	BIOBIN_ROLE_ARG=""
	if [[ -n "$role_file" ]]; then
		dx download "$role_file" -o input/input.role
		BIOBIN_ROLE_ARG="--role-file input/input.role"
	fi
	BIOBIN_PHENO_ARG=""
	BIOBIN_TEST_ARG=""
	if [[ -n "$phenotype_file" ]]; then
		dx download "$phenotype_file" -o input/input.phenotype
		BIOBIN_PHENO_ARG="--phenotype-file input/input.phenotype"
		if [[ ${#regression_type[*]} -gt 0 ]]; then
			BIOBIN_TEST_ARG="--test $(IFS="," ; echo "${regression_type[*]}")"
		fi
	fi
	BIOBIN_COVAR_ARG=""
	if [[ -n "$covariate_file" ]]; then
		dx download "$covariate_file" -o input/input.covariate
		BIOBIN_COVAR_ARG="--covariates input/input.covariate"
	fi
	BIOBIN_WEIGHT_ARG=""
	if [[ -n "$weight_file" ]]; then
		dx download "$weight_file" -o input/input.weight
		BIOBIN_WEIGHT_ARG="--weight-file input/input.weight"
	fi
	BIOBIN_REGION_ARG=""
	if [[ -n "$region_file" ]]; then
		dx download "$region_file" -o input/input.region
		BIOBIN_REGION_ARG="--region-file input/input.region"
	fi
	BIOBIN_INCLUDE_REGION_ARG=""
	if [[ -n "$include_region_file" ]]; then
		dx download "$include_region_file" -o input/input.gene
		BIOBIN_INCLUDE_REGION_ARG="--include-region-file input/input.gene"
	fi
	VCF_FILE="input/input.vcf.gz"
	TBI_FILE="$VCF_FILE.tbi"
	dx download "$vcf_file" -o "$VCF_FILE"
	if [[ -n "$vcf_tbi_file" ]]; then
		dx download "$vcf_tbi_file" -o "$TBI_FILE"
	else
		tabix -p vcf "$VCF_FILE" 2>&1 | tee -a output.log
	fi
	
	
	# run biobin
	
	loki_db="$(ls -1 /opt/biobin/LOKI*.db | sort | tail -n 1)"
	mkdir biobin
	/opt/biobin/biobin \
		--threads "$NUM_CORES" \
		--settings-db "$loki_db" \
		--vcf-file "$VCF_FILE" \
		$BIOBIN_ROLE_ARG \
		$BIOBIN_PHENO_ARG \
		$BIOBIN_COVAR_ARG \
		$BIOBIN_TEST_ARG \
		$BIOBIN_WEIGHT_ARG \
		$BIOBIN_REGION_ARG \
		$BIOBIN_INCLUDE_REGION_ARG \
		--report-prefix "biobin/$output_prefix" \
		$biobin_args \
	2>&1 | tee -a output.log
	ls -laR biobin
	
	
	# run summary script
	
	RGX="--force-all-control[[:space:]]+(y|Y)"
	ALL_CONTROL=""
	if [[ $biobin_args =~ $RGX ]]; then
		ALL_CONTROL="--all-control"
	fi
	python /opt/biobin/biobin-summary.py \
		--prefix="biobin/$output_prefix" \
		$ALL_CONTROL \
		> "biobin/${output_prefix}-summary.tsv"
	
	
	# upload output files
	
	for f in biobin/*-bins.csv ; do
		bin_file=$(dx upload "$f" --brief)
		dx-jobutil-add-output bins_files --class="array:file" "$bin_file"
	done
	
	summary_file=$(dx upload biobin/*-summary.tsv --brief)
	dx-jobutil-add-output summary_file --class="file" "$summary_file"
	
	locus_file=$(dx upload biobin/*-locus.csv --brief)
	dx-jobutil-add-output locus_file --class="file" "$locus_file"
	
	mv output.log "${output_prefix}.log"
	log_file=$(dx upload "${output_prefix}.log" --brief)
	dx-jobutil-add-output log_file --class="file" "$log_file"
	
	for f in biobin/*-unlifted.csv ; do
		if [[ -f "$f" ]]; then
			unlifted_file=$(dx upload "$f" --brief)
			dx-jobutil-add-output unlifted_files --class="array:file" "$unlifted_file"
		fi
	done
	
}
