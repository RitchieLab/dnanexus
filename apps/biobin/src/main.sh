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
		BIOBIN_TEST_ARG="--test $(IFS="," ; echo "${regression_type[*]}")"
	fi
	BIOBIN_COVAR_ARG=""
	if [[ -n "$covariate_file" ]]; then
		dx download "$covariate_file" -o input/input.covariate
		BIOBIN_COVAR_ARG="--covariates input/input.covariate"
	fi
	VCF_FILE="input/input.vcf.gz"
	TBI_FILE="$VCF_FILE.tbi"
	dx download "$vcf_file" -o "$VCF_FILE"
	if [[ -n "$vcf_tbi_file" ]]; then
		dx download "$vcf_tbi_file" -o "$TBI_FILE"
	else
		tabix -p vcf "$VCF_FILE" 2>&1 | tee -a output.log
	fi
	
	
	# download shared resource files
	
#	DX_RESOURCES_ID="$(dx find projects --name "App Resources" --brief)"
	DX_RESOURCES_ID="project-BYpFk1Q0pB0xzQY8ZxgJFv1V"
	mkdir shared
	dx download \
		"$(dx find data \
			--name "loki-20150427-nosnps.db" \
			--project "$DX_RESOURCES_ID" \
			--folder /LOKI \
			--brief \
		)" -o shared/loki.db
	
	
	# run biobin
	
	mkdir biobin
	biobin \
		--threads "$NUM_CORES" \
		--settings-db shared/loki.db \
		--vcf-file "$VCF_FILE" \
		$BIOBIN_ROLE_ARG \
		$BIOBIN_PHENO_ARG \
		$BIOBIN_COVAR_ARG \
		$BIOBIN_TEST_ARG \
		--report-prefix biobin/output \
		$biobin_args \
	2>&1 | tee -a output.log
	ls -laR biobin
	
	
	# archive output
	
	tar cjvf biobin.tar.bz biobin 2>&1 | tee -a output.log
	
	
	# move files into position for upload
	mkdir -p "$HOME/out/log_file"
	mv output.log "$HOME/out/log_file"
	mkdir -p "$HOME/out/biobin_output"
	mv biobin.tar.bz "$HOME/out/biobin_output"
	
	
	# return to the home dir and upload all files
	cd "$HOME"
	dx-upload-all-outputs
}
