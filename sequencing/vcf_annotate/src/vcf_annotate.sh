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
# * The ability to install and run VEP is left in-tact, but it it strongly recommended
#   to instead make VEP a separate app. Sources for a VEP app are available here:
#   https://github.com/MelbourneGenomics/dx-vep.
# * It would be simpler, easier, and faster to use vcfanno.
set -x -e -o pipefail

export SHELL="/bin/bash"

procCount=$(nproc --all)

# build bcftools
cd /home/dnanexus/bcftools-1.3.1
make
make prefix=/usr/local/ install

# create temporary working directories
WKDIR=$(mktemp -d)
DXVCF_LIST=$(mktemp)
DXIDX_LIST=$(mktemp)

cd ${WKDIR}


#######################################
# Download a file to a directory such that multiple downloads can be run at the same
# time (by gnu parallel).
# Globals:
#   None
# Arguments:
#   1. The file to download
#   2. The directory in which to save the downloaded file
# Returns:
#   None
#######################################
function parallel_download() {
    echo "parallel_download"
    echo "$1"
    echo "$2"
	cd $2
	dx download "$1"
	cd - > /dev/null
}
export -f parallel_download


#######################################
# Download source databases for annotation based on which environment variables are set.
# Globals:
#   vep_vcfgz, vep_vcfgztbi, vep_header: annotation files for Variant Effect Predictor
#   dbnsfp_vcfgz, dbnsfp_vcfgztbi: anntation files for dbNSFP
#   hgmd_vcfgz, hgmd_vcfgztbi: annotation files for Human Gene Mutation Database
#   clinvar_vcfgz, clinvar_vcfgztbi: annotation files for ClinVar
#   HOME: home directory
#   WKDIR: working directory
# Arguments:
#   None
# Returns:
#   None
#######################################
function download_resources(){

	cd ${WKDIR}

	if [ -n "$vep_vcfgz" ]; then

        dx download "$vep_vcfgz" -o VEP.vcf.gz
        dx download "$vep_vcfgztbi" -o VEP.vcf.gz.tbi
        if [ -n "$vep_header_txt" ]; then
          dx download "$vep_header_txt" -o VEP_header.txt
        fi

    elif [ -n "$vep_source_zip" ]; then

        # Download VEP dependencies, executable, and database
        cpanm Archive::Extract
        cpanm File::Copy::Recursive
        cpanm Module::Build
        cpanm Archive::Zip
        cpanm Bundle::DBI

        dx download "$vep_source_zip" -o /usr/share/ensembl-tools-release.zip
        cd /usr/share/
        unzip ensembl-tools-release.zip
        cd ensembl-tools-release-87/scripts/variant_effect_predictor/
        yes n | perl INSTALL.pl

        dx download "$vep_cache_targz" -o $HOME/vep_b38.tar.gz
        cd $HOME
        tar -xzvf vep_b38.tar.gz

        sudo mkdir -p /usr/share/GATK/resources
        sudo chmod -R a+rwX /usr/share/GATK

        dx download "$DX_RESOURCES_ID:/GATK/resources/human_g1k_v37_decoy.fasta" -o /usr/share/GATK/resources/build.fasta
        dx download "$DX_RESOURCES_ID:/GATK/resources/human_g1k_v37_decoy.fasta.fai" -o /usr/share/GATK/resources/build.fasta.fai
        dx download "$DX_RESOURCES_ID:/GATK/resources/human_g1k_v37_decoy.dict" -o /usr/share/GATK/resources/build.dict

        cd ${WKDIR}

	fi

	if [ -n "$dbnsfp_vcfgz" ]; then

		dx download "$dbnsfp_vcfgz" -o dbNSFP.vcf.gz
		dx download "$dbnsfp_vcfgztbi" -o dbNSFP.vcf.gz.tbi
		if [ -n "$dbnsfp_header_txt" ]; then
		  dx download "$dbnsfp_header_txt" -o dbNSFP_header.txt
		fi

	fi

	if [ -n "$hgmd_vcfgz" ]; then

		dx download "$hgmd_vcfgz" -o HGMD.vcf.gz
		dx download "$hgmd_vcfgztbi" -o HGMD.vcf.gz.tbi
		if [ -n "$hgmd_header_txt" ]; then
		  dx download "$hgmd_header_txt" -o HGMD_header.txt
        fi
	fi

	if [ -n "$clinvar_vcfgz" ]; then

		dx download "$clinvar_vcfgz" -o variant_summary.vcf.gz
		dx download "$clinvar_vcfgztbi" -o variant_summary.vcf.gz.tbi
		if [ -n "$clinvar_header_txt" ]; then
		  dx download "$clinvar_header_txt" -o ClinVar_header.txt
		fi

	fi
}
export -f download_resources


#######################################
# Download a VCF file and apply annotations. Annotations are applied serially, with the
# output of the previous step (OUT_VCF) being used as the input to the next step
# (IN_VCF). All the intermediate files are deleted.
# Globals:
#   vep_vcfgz, vep_vcfgztbi, vep_header: annotation files for Variant Effect Predictor
#   dbnsfp_vcfgz, dbnsfp_vcfgztbi: anntation files for dbNSFP
#   hgmd_vcfgz, hgmd_vcfgztbi: annotation files for Human Gene Mutation Database
#   clinvar_vcfgz, clinvar_vcfgztbi: annotation files for ClinVar
# Arguments:
#   1. The VCF file to download
#   2. The directory in which to save the downloaded file
# Returns:
#   None
#######################################
function parallel_download_and_annotate() {

    # Download the VCF file
	parallel_download "$1" "$2"

    # Apply the annotations

    procCount=$(nproc --all)

	IN_VCF=$(dx describe "$1" --name)

	OUT_VCF=${IN_VCF}

	if [ -n "$vep_vcfgz" ]; then

        OUT_VCF=${IN_VCF%.vcf.gz}.VEP.vcf.gz

        tabix -p vcf -f ${IN_VCF}

        echo "VEP ${IN_VCF} -> ${OUT_VCF}"

        if [ -n "$vep_header_txt" ]; then

          bcftools annotate -a VEP.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO/CSQ -h VEP_header.txt

        else

          bcftools annotate -a VEP.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO/CSQ

        fi

        rm ${IN_VCF}*

    elif [ -n "$vep_source_zip" ]; then

        OUT_VCF=${IN_VCF%.vcf.gz}.VEP.u.vcf

        perl /usr/share/ensembl-tools-release-87/scripts/variant_effect_predictor/variant_effect_predictor.pl -i ${IN_VCF} --everything --cache --offline --vcf --fasta /usr/share/GATK/resources/build.fasta --merged -o ${OUT_VCF} --no_progress --no_stats --force_overwrite

        rm ${IN_VCF}*

        vcf-sort ${OUT_VCF} | bgzip > ${OUT_VCF%.u.vcf}.vcf.gz

        rm ${OUT_VCF}*

        OUT_VCF=${OUT_VCF%.u.vcf}.vcf.gz

	fi

	if [ -n "$dbnsfp_vcfgz" ]; then

		IN_VCF=${OUT_VCF}
		tabix -p vcf -f ${IN_VCF}

		OUT_VCF=${IN_VCF%.vcf.gz}.dbNSFP.vcf.gz

		echo "dbnsfp ${IN_VCF} -> ${OUT_VCF}"

        if [ -n "$dbnsfp_header_txt" ]; then

		  bcftools annotate -a dbNSFP.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO -h dbNSFP_header.txt

        else

          bcftools annotate -a dbNSFP.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO

        fi

		rm ${IN_VCF}*

	fi

	if [ -n "$hgmd_vcfgz" ]; then

		IN_VCF=${OUT_VCF}
		tabix -p vcf -f ${IN_VCF}

		OUT_VCF=${IN_VCF%.vcf.gz}.HGMD.vcf.gz

		echo "HGMD ${IN_VCF} -> ${OUT_VCF}"

        if [ -n "$hgmd_header_txt" ]; then

		  bcftools annotate -a HGMD.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO -h HGMD_header.txt

        else

          bcftools annotate -a HGMD.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO

        fi

		rm ${IN_VCF}*

	fi

	if [ -n "$clinvar_vcfgz" ]; then

		IN_VCF=${OUT_VCF}
		tabix -p vcf -f ${IN_VCF}

		OUT_VCF=${OUT_VCF%.vcf.gz}.ClinVar.vcf.gz

		echo "clinvar ${IN_VCF} -> ${OUT_VCF}"

        if [ -n "$clinvar_header_txt" ]; then

		  bcftools annotate -a variant_summary.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO -h ClinVar_header.txt

        else

          bcftools annotate -a variant_summary.vcf.gz -o ${OUT_VCF} -Oz ${IN_VCF} -c +INFO

        fi

		rm ${IN_VCF}*

	fi

    tabix -p vcf -f ${OUT_VCF}

    # Upload the results

	VCF_UP=$(dx upload --brief ${OUT_VCF})
	IDX_UP=$(dx upload --brief ${OUT_VCF}.tbi)

    dx-jobutil-add-output out_variants_vcfgz "$VCF_UP" --class=array:file
	dx-jobutil-add-output out_variants_vcfgztbi "$IDX_UP" --class=array:file

	TO_RM=$(dx describe "$1" --name)

	rm ${TO_RM%.vcf.gz}*

}
export -f parallel_download_and_annotate


#######################################
# Main function.
# Globals:
#   variants_vcfgz: The index file(s) to download; may be an array
#   variants_vcfgztbi: The VCF file(s) to download; may be an array of the same length as
#     variants_vcfgz
#   HOME: home directory
#   WKDIR: working directory
# Arguments:
#   None
# Returns:
#   None; The output VCF and index files that were successfully annotated are appended
#   to the out_variants_vcfgz and out_variants_vcfgztbi list variables.
#######################################
main() {

    echo "Value of variants_vcfgz: '$variants_vcfgz'"
    echo "Value of variants_vcfgztbi: '$variants_vcfgztbi'"
    echo "Value of vep_vcfgz: '$vep_vcfgz'"
    echo "Value of dbnsfp_vcfgz: '$dbnsfp_vcfgz'"
    echo "Value of hgmd_vcfgz: '$hgmd_vcfgz'"
    echo "Value of clinvar_vcfgz: '$clinvar_vcfgz'"

	export SHELL="/bin/bash"

    # TODO: I'm not sure why cd HOME and then WKDIR, unless WKDIR might be relative
	cd ${HOME}
    cd ${WKDIR}

    # Download the VCF index files (in parallel)
    for i in "${!variants_vcfgztbi[@]}"; do
      echo "${variants_vcfgztbi[$i]}" >> "$DXIDX_LIST"
    done

    parallel -j $(nproc --all) -u --gnu parallel_download :::: ${DXIDX_LIST} ::: ${WKDIR}

    # Download the annotation databases
	download_resources

    echo "VCF indexes downloaded"

	procCount=$(nproc --all)
	halfCount=$(($procCount/2))

	# Download the VCF files (in parallel)
    for i in "${!variants_vcfgz[@]}"; do
      echo "${variants_vcfgz[$i]}" >> ${DXVCF_LIST}
    done

    parallel -j ${halfCount} -u --gnu parallel_download_and_annotate :::: ${DXVCF_LIST} ::: ${WKDIR}

    echo "VCF files downloaded and annotated"

}
