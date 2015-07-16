#!/bin/bash

set -e -x -o pipefail

main() {
	
	# switch to a temp directory and download all non-VCF user input files
	
	TMPDIR="$(mktemp -d)"
	cd "$TMPDIR"
	mkdir input
	if [[ -n "$region_file" ]]; then
		dx download "$region_file" -o input/input.regions
	fi
	if [[ -n "$sample_file" ]]; then
		dx download "$sample_file" -o input/input.samples
	fi
	
	
	# locate the shared resource project
	
	DX_RESOURCES_ID="$(dx find projects --name "GATK Resources" --brief)"
	DX_RESOURCES_ID="project-BYpFk1Q0pB0xzQY8ZxgJFv1V"
	mkdir shared
	
	
	# if there's a region_file, run download_part.py to grab only the parts we need
	
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
			TBI_FILE="input/input.vcf.gz.tbi"
		else
			VCF_FILE="input/input.vcf"
			TBI_FILE="input/input.vcf.tbi"
		fi
		dx download "$vcf_file" -o "$VCF_FILE"
		dx download "$vcf_tbi_file" -o "$TBI_FILE"
	fi
	
	
	# if there's a sample_file, use it to further filter the VCF with GATK
	
	if [[ -f input/input.samples ]]; then
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
			--sample_file input/input.samples \
			--excludeFiltered \
			--out gatk/output.vcf.gz \
		2>&1 | tee output.log
		ls -laR gatk
		VCF_FILE="gatk/output.vcf.gz"
		TBI_FILE="gatk/output.vcf.gz.tbi"
	fi
	
	
	# generate the TBI if we don't have one already
	if [[ -z "$TBI_FILE" ]]; then
		tabix -p vcf "$VCF_FILE" 2>&1 | tee -a output.log
		TBI_FILE="$VCF_FILE.tbi"
	fi
	
	
	# convert the VCF to plink files
	
	mkdir vcftools
	if [[ "$VCF_FILE" == *".gz" ]]; then
		VCFTOOLS_OPT="--gzvcf"
	else
		VCFTOOLS_OPT="--vcf"
	fi
	vcftools \
		$VCFTOOLS_OPT "$VCF_FILE" \
		--plink-tped \
		--out vcftools/output \
	2>&1 | tee -a output.log
	ls -laR vcftools

    # convert BIM file to add chromosome to unlabeled variants
    mv vcftools/output.tped vcftools/original.tped
    gawk \
        '{ if ($2==$4) $2=("" $1 ":" $2); print }' \
        vcftools/original.tped \
    > vcftools/output.tped


	mkdir plink
	if [[ -n "$sample_file" ]]; then
		paste input/input.samples input/input.samples > plink/input.samples
		PLINK_KEEP_ARG="--keep plink/input.samples"
	else
		PLINK_KEEP_ARG=""
	fi
	/usr/lib/plink/plink \
		--tfile vcftools/output \
		$PLINK_KEEP_ARG \
		--make-bed \
		--out plink/output \
	2>&1 | tee -a output.log
	ls -laR plink
	
	
	# archive the output so far
	
	tar cjvf vcftools.tar.bz vcftools 2>&1 | tee -a output.log
	tar cjvf plink.tar.bz plink 2>&1 | tee -a output.log
	
	
	# create output directories
	
	
	# move files into position for upload
	mkdir -p "$HOME/out/log_file"
	mv output.log "$HOME/out/log_file"
	mkdir -p "$HOME/out/vcftools_output"
	mv vcftools.tar.bz "$HOME/out/vcftools_output"
	mkdir -p "$HOME/out/plink_output"
	mv plink.tar.bz "$HOME/out/plink_output"
	
	
	# return to the home dir and upload all files
	cd "$HOME"
	dx-upload-all-outputs
}
