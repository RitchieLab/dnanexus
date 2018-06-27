#!/bin/bash
# vcf_annotate 0.1.0
#
# * App runs on a single machine from beginning to end.
# * Job input variables are loaded as environment variables before this script runs.
# Any array inputs will be loaded as bash arrays.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.
#
# TODO:
# * dxapp.json defines $build_version, but this script is hardcoded to use build 37.
# * It would be simpler, easier, and faster to use vcfanno.
#
# Environment variables (defined in dxapp.json)
# ---------------------------------------------
# vcfidx_fn: The index file(s) to download; may be an array
# vcf_fn: The VCF file(s) to download; may be an array of the same length as vcfidx_fn
# VEP: annotate using Variant Effect Predictor (VEP)
# annotate_header: If true, downloads VEP and all the supporting VEP data files and
#   annotates the VCF by calling VEP. Otherwise, downloads a VCF with pre-computed
#   annotations (stored in the INFO/CSQ field) and applies them to the query VCF using
#   bcftools annotate. [TODO: maybe name this variable something more informative]
# dbnsfp: annotate with dbNSFP (compendium of functional prediction scores for
#   non-synonymous variants
# HGMD: annotate with Human Gene Mutation Database (HGMD)
# clinvar: annotate with ClinVar (database of clinically relevant variants)
#
# Result
# ------
# The output VCF and index files that were successfully annotated are appended to the
# vcf_out and vcfidx_out list variables.

# Set default values for environment variables
: "${VEP:=false}"
: "${annotate_header:=false}"
: "${dbnsfp:=true}"
: "${HGMD:=true}"
: "${clinvar:=true}"

# TODO: are there any reasonable defaults for these?
: "${vcf_fn:=()}"
: "${vcfidx_fn:=()}"
: "${DX_RESOURCES_ID:=}"

export SHELL="/bin/bash"

procCount=$(nproc --all)

set -x -e -o pipefail

# build bcftools
cd /home/dnanexus/bcftools-1.3.1
make
make prefix=/usr/local/ install

# create temporary working directories
WKDIR=$(mktemp -d)
DXVCF_LIST=$(mktemp)
DXIDX_LIST=$(mktemp)

cd ${WKDIR}


# This function downloads a file (first argument) to a directory (second argument)
# such that multiple downloads can be run at the same time (by gnu parallel).
function parallel_download() {
	cd $2
	dx download "$1"
	cd - >/dev/null
}
export -f parallel_download


# This function downloads source databases for annotation based on which environment
# variables are set.
function download_resources(){

	cd ${WKDIR}

	if test "$VEP" = "true"; then

		if "$annotate_header" = "true"; then

		    # Download VEP dependencies, executable, and database

			cpanm Archive::Extract
			cpanm File::Copy::Recursive
			cpanm Module::Build
			cpanm Archive::Zip
			cpanm Bundle::DBI

			dx download "$DX_RESOURCES_ID:/VEP/ensembl-tools-release-87.zip" -o /usr/share/ensembl-tools-release.zip
			cd /usr/share/
			unzip ensembl-tools-release.zip
			cd ensembl-tools-release-87/scripts/variant_effect_predictor/
			yes n | perl INSTALL.pl

			dx download "$DX_RESOURCES_ID:/VEP/vep_b38.tar.gz" -o $HOME/vep_b38.tar.gz
			cd $HOME

			tar -xzvf vep_b38.tar.gz

			sudo mkdir -p /usr/share/GATK/resources
			sudo chmod -R a+rwX /usr/share/GATK

			dx download "$DX_RESOURCES_ID:/GATK/resources/human_g1k_v37_decoy.fasta" -o /usr/share/GATK/resources/build.fasta
			dx download "$DX_RESOURCES_ID:/GATK/resources/human_g1k_v37_decoy.fasta.fai" -o /usr/share/GATK/resources/build.fasta.fai
			dx download "$DX_RESOURCES_ID:/GATK/resources/human_g1k_v37_decoy.dict" -o /usr/share/GATK/resources/build.dict

            cd ${WKDIR}

		else

            # Just download the pre-computed VEP VCF
			dx download "$DX_RESOURCES_ID:/VEP/combined.padded.header.vcf.gz" -o VEP.vcf.gz
			dx download "$DX_RESOURCES_ID:/VEP/combined.padded.header.vcf.gz.tbi" -o VEP.vcf.gz.tbi
			dx download "$DX_RESOURCES_ID:/VEP/VEP_header.txt" -o VEP_header.txt

		fi

	fi

	if test "$dbnsfp" = "true"; then

		dx download "$DX_RESOURCES_ID:/dbNSFP/3.x/dbNSFP3.3a.vcf.gz" -o dbNSFP.vcf.gz
		dx download "$DX_RESOURCES_ID:/dbNSFP/3.x/dbNSFP3.3a.vcf.gz.tbi" -o dbNSFP.vcf.gz.tbi
		dx download "$DX_RESOURCES_ID:/dbNSFP/3.x/dbNSFP3.3a_header.txt" -o dbNSFP_header.txt

	fi

	if test "$HGMD" = "true"; then

		dx download "$DX_RESOURCES_ID:/HGMD/HGMD_PRO_2016.4_hg38.vcf.gz" -o HGMD.vcf.gz
		dx download "$DX_RESOURCES_ID:/HGMD/HGMD_PRO_2016.4_hg38.vcf.gz.tbi" -o HGMD.vcf.gz.tbi
		dx download "$DX_RESOURCES_ID:/dbNSFP/3.x/dbNSFP3.3a_header.txt" -o HGMD_header.txt

	fi

	if test "$clinvar" = "true"; then

		dx download "$DX_RESOURCES_ID:/CLINVAR/variant_summary.Jan2017.b38.vcf.gz" -o variant_summary.vcf.gz
		dx download "$DX_RESOURCES_ID:/CLINVAR/variant_summary.Jan2017.b38.vcf.gz.tbi" -o variant_summary.vcf.gz.tbi
		dx download "$DX_RESOURCES_ID:/CLINVAR/ClinVar_Jan2017_header.txt" -o ClinVar_header.txt

	fi
}
export -f download_resources


# This function downloads a VCF file and applies annotations. Annotations are applied
# serially, with the output of the previous step (OUT_VCF) being used as the input
# to the next step (IN_VCF). All the intermediate files are deleted.
function parallel_download_and_annotate() {

    # Download the VCF file
	parallel_download $1 $2

    # Apply the annotations

    procCount=$(nproc --all)

	IN_VCF=$(dx describe "$1" --name)

	OUT_VCF=${IN_VCF}

	if test "$VEP" = "true"; then

		if "$annotate_header" = "true"; then

			OUT_VCF=${IN_VCF%.vcf.gz}.VEP.u.vcf

			perl /usr/share/ensembl-tools-release-87/scripts/variant_effect_predictor/variant_effect_predictor.pl -i ${IN_VCF} --everything --cache --offline --vcf --fasta /usr/share/GATK/resources/build.fasta --merged -o ${OUT_VCF} --no_progress --no_stats --force_overwrite

			rm ${IN_VCF}*

			vcf-sort ${OUT_VCF} | bgzip > ${OUT_VCF%.u.vcf}.vcf.gz

			rm ${OUT_VCF}*

			OUT_VCF=${OUT_VCF%.u.vcf}.vcf.gz

		else

			OUT_VCF=${IN_VCF%.vcf.gz}.VEP.vcf.gz

			tabix -p vcf -f ${IN_VCF}

			echo "VEP"
			echo ${OUT_VCF}
			echo ${IN_VCF}

			bcftools annotate -a VEP.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO/CSQ -h VEP_header.txt

			rm ${IN_VCF}*

		fi

	fi

	if test "$HGMD" = "true"; then

		IN_VCF=${OUT_VCF}
		tabix -p vcf -f ${IN_VCF}

		OUT_VCF=${IN_VCF%.vcf.gz}.HGMD.vcf.gz

		echo "HGMD"
		echo ${OUT_VCF}
		echo ${IN_VCF}

		bcftools annotate -a HGMD.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO -h HGMD_header.txt

		rm ${IN_VCF}*

	fi

	if test "$dbnsfp" = "true"; then

		IN_VCF=${OUT_VCF}
		tabix -p vcf -f ${IN_VCF}

		OUT_VCF=${IN_VCF%.vcf.gz}.dbNSFP.vcf.gz

		echo "dbnsfp"
		echo ${OUT_VCF}
		echo ${IN_VCF}

		bcftools annotate -a dbNSFP.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO -h dbNSFP_header.txt

		rm ${IN_VCF}*

	fi

	if test "$clinvar" = "true"; then

		IN_VCF=${OUT_VCF}
		tabix -p vcf -f ${IN_VCF}

		OUT_VCF=${OUT_VCF%.vcf.gz}.ClinVar.vcf.gz

		echo "clinvar"
		echo ${OUT_VCF}
		echo ${IN_VCF}

		bcftools annotate -a variant_summary.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO -h ClinVar_header.txt
		rm ${IN_VCF}*
		tabix -p vcf -f ${OUT_VCF}

	fi

    # Upload the results

	VCF_UP=$(dx upload --brief ${OUT_VCF})
	IDX_UP=$(dx upload --brief${OUT_VCF}.tbi)

    dx-jobutil-add-output vcf_out "$VCF_UP" --class=array:file
	dx-jobutil-add-output vcfidx_out "$IDX_UP" --class=array:file

	TO_RM=$(dx describe "$1" --name)

	rm ${TO_RM%.vcf.gz}*

}
export -f parallel_download_and_annotate


main() {

	export SHELL="/bin/bash"

	cd ${HOME}

    echo "Value of vcf_fn: '$vcf_fn'"
    echo "Value of vcfidx_fn: '$vcfidx_fn'"

    cd ${WKDIR}

    # Download the VCF index files (in parallel)
    for i in "${!vcfidx_fn[@]}"; do
      echo "${vcfidx_fn[$i]}" >> "$DXIDX_LIST"
    done

    parallel -j $(nproc --all) -u --gnu parallel_download :::: ${DXIDX_LIST} ::: ${WKDIR}

    # Download the annotation databases
	download_resources

    echo "VCF indexes downloaded"

	procCount=$(nproc --all)
	halfCount=$(($procCount/2))

	# Download the VCF files (in parallel)
    for i in "${!vcf_fn[@]}"; do
      echo "${vcf_fn[$i]}" >> ${DXVCF_LIST}
    done

    parallel -j ${halfCount} -u --gnu parallel_download_and_annotate :::: ${DXVCF_LIST}::: ${WKDIR}

    echo "VCF files downloaded and annotated"

}
