#!/bin/bash

set -e -x -o pipefail

main() {
	# update to newer boost (for biobin)
	# Don't need this now!!
	#sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	#sudo add-apt-repository ppa:boost-latest/ppa --yes
	#sudo apt-get update --yes
	#sudo apt-get install --yes libboost1.54-all-dev
	
	
	# switch to a temp directory and download all non-VCF user input files
	
	NUM_CORES="$(nproc --all)"
	TMPDIR="$(mktemp -d)"
	cd "$TMPDIR"
	mkdir input
	BIOBIN_PHENO_ARG=""
	if [[ -n "$phenotype_file" ]]; then
		dx download "$phenotype_file" -o input/input.phenotype
		BIOBIN_PHENO_ARG="--phenotype-file input/input.phenotype --test $regression_type,SKAT-$regression_type"
	fi
	BIOBIN_COVAR_ARG=""
	if [[ -n "$covariate_file" ]]; then
		dx download "$covariate_file" -o input/input.covariate
		BIOBIN_COVAR_ARG="--covariates input/input.covariate"
	fi
	if [[ -n "$region_file" ]]; then
		dx download "$region_file" -o input/input.regions
	fi
	if [[ -n "$sample_file" ]]; then
		dx download "$sample_file" -o input/input.samples
	fi
	
	
	# download shared resource files
	
	DX_RESOURCES_ID="$(dx find projects --name "GATK Resources" --brief)"
	DX_RESOURCES_ID="project-BYpFk1Q0pB0xzQY8ZxgJFv1V"
	mkdir shared
	dx download \
		"$(dx find data \
			--name "loki-20150427-nosnps.db" \
			--project "$DX_RESOURCES_ID" \
			--folder /resources \
			--brief \
		)" -o shared/loki.db
	
	
	# if there's a regions file, make sure it's not too big
	
	if [[ -f input/input.regions ]]; then
		if [[ $(cat input/input.regions | wc -l) -gt 100 ]]; then
			mv input/input.regions input/input.list
		fi
	fi
	
	
	# if there's still a (small) regions file, run download_part.py to grab only the parts we need
	
	if [[ -f input/input.regions ]]; then
		VCF_HASH="$(dx describe "$vcf_file" --json | jq -r .id)"
		TBI_HASH="$(dx describe "$vcf_tbi_file" --json | jq -r .id)"
		
		VCF_FILE="input/input.vcf.gz"
		TBI_FILE=""
		
		PART_ARG="--header --keep-open"
		for INTERVAL in $(head -n -1 input/input.regions) ; do
			python /usr/share/download_part.py \
				--interval "$INTERVAL" \
				--vcf "$VCF_HASH" \
				--index "$TBI_HASH" \
				$PART_ARG \
				--output "$VCF_FILE" \
			2>&1 | tee -a output.log
			PART_ARG="--append --keep-open"
		done
		INTERVAL="$(tail -n 1 input/input.regions)"
		if [[ -f "$VCF_FILE" ]]; then
			PART_ARG="--append"
		else
			PART_ARG="--header"
		fi
		python /usr/share/download_part.py \
			--interval "$INTERVAL" \
			--vcf "$VCF_HASH" \
			--index "$TBI_HASH" \
			$PART_ARG \
			--output "$VCF_FILE" \
		2>&1 | tee -a output.log
	else
		if [[ "$(dx describe "$vcf_file" --name)" == *".gz" ]]; then
			VCF_FILE="input/input.vcf.gz"
		else
			VCF_FILE="input/input.vcf"
		fi
		dx download "$vcf_file" -o "$VCF_FILE"
		
		if [[ -n "$vcf_tbi_file" ]]; then
			TBI_FILE="$VCF_FILE.tbi"
			dx download "$vcf_tbi_file" -o "$TBI_FILE"
		else
			TBI_FILE=""
		fi
	fi
	
	
	# if there's a samples or list file, use it to further filter the VCF with GATK
	
	GATK_SAMPLES_ARG=""
	if [[ -f input/input.samples ]]; then
		GATK_SAMPLES_ARG="--sample_file input/input.samples"
	fi
	GATK_LIST_ARG=""
	if [[ -f input/input.list ]]; then
		GATK_LIST_ARG="--intervals input/input.list"
	fi
	if [[ -f input/input.samples || -f input/input.list ]]; then
		# generate the TBI if we don't have one already
		if [[ -z "$TBI_FILE" ]]; then
			tabix -p vcf "$VCF_FILE" 2>&1 | tee -a output.log
			TBI_FILE="$VCF_FILE.tbi"
		fi
		
		# download GATK files
		dx download "$(dx find data --brief \
			--name "GenomeAnalysisTK-3.3-0.jar" \
			--project "$DX_RESOURCES_ID" \
			--folder /jar \
		)" -o shared/GATK.jar
		dx download "$(dx find data --brief \
			--name "human_g1k_v37_decoy.fasta" \
			--project "$DX_RESOURCES_ID" \
			--folder /resources \
		)" -o shared/reference.fasta
		dx download "$(dx find data --brief \
			--name "human_g1k_v37_decoy.fasta.fai" \
			--project "$DX_RESOURCES_ID" \
			--folder /resources \
		)" -o shared/reference.fasta.fai
		dx download "$(dx find data --brief \
			--name "human_g1k_v37_decoy.dict" \
			--project "$DX_RESOURCES_ID" \
			--folder /resources \
		)" -o shared/reference.dict
		
		# run GATK
		mkdir gatk
		java -jar shared/GATK.jar \
			--analysis_type SelectVariants \
			--variant "$VCF_FILE" \
			--reference_sequence shared/reference.fasta \
			--excludeFiltered \
			--ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES \
			$GATK_SAMPLES_ARG \
			$GATK_LIST_ARG \
			--out gatk/output.vcf.gz \
		2>&1 | tee output.log
		ls -laR gatk
		VCF_FILE="gatk/output.vcf.gz"
		TBI_FILE=""
	fi
	
	
	# run biobin
	
	# decmpress on the fly
#	CAT_CMD="cat"
#	if [[ "$VCF_FILE" == *".gz" ]]; then
#		CAT_CMD="zcat"
#	fi
	
	mkdir biobin
	biobin \
		--threads "$NUM_CORES" \
		--settings-db shared/loki.db \
		--vcf-file "$VCF_FILE" \
		$BIOBIN_PHENO_ARG \
		$BIOBIN_COVAR_ARG \
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
