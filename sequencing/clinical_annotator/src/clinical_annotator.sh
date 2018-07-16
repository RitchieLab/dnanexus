#!/bin/bash
# clinical_annotator 1.0.0
set -e -x -o pipefail

#################################################
# Downloads input VCF index files
# Arguments: 1) VCF Index File 2) /home/dnanexus
# Returns: None
#################################################
function parallel_download() {
	cd $2
	dx download "$1"
	cd - >/dev/null
}
export -f parallel_download

###########################################
# Downloads and reformats resource files
# Arguments: None
# Returns: None
############################################
function reformat_resources() {

    cd $HOME

    python ClinVar_tsvTOvcf.py variant_summary.txt.gz
    
    if [ "$build_version" = "b37" ];
    then
        vcf-sort -c variant_summary.b37.vcf | bgzip -c > variant_summary.b37.vcf.gz
        tabix -p vcf -f variant_summary.b37.vcf.gz
    else
        vcf-sort -c variant_summary.b38.vcf | bgzip -c > variant_summary.b38.vcf.gz
        tabix -p vcf -f variant_summary.b38.vcf.gz
    fi

    if test "${hgmd_pro_file}"; then

        dx download "${hgmd_pro_file}"
        python reformatHGMD.py $(dx describe "$hgmd_pro_file" --name) | bgzip -c > $HOME/HGMD_PRO.reformated.vcf.gz
        tabix -p vcf -f HGMD_PRO.reformated.vcf.gz
	
    fi

}
export -f reformat_resources

###########################################
# Downloads input VCF files and Annotates
# Arguments: 1) VCF File 2) /home/dnanexus
# Returns: None
############################################
function parallel_download_and_annotate() {
	
    cd $2
    dx download "$1"

    IN_VCF=$(dx describe "$1" --name)

    OUT_VCF=$IN_VCF

    # CLINVAR
    if [ ! -f $IN_VCF.tbi ]; then
        
        tabix -p vcf -f $IN_VCF

    fi

    OUT_VCF=${OUT_VCF%.vcf.gz}.ClinVar.vcf.gz

    echo "Annotating ClinVar"
    echo $OUT_VCF
    echo $IN_VCF

    if [ "$build_version" = "b37" ];
    then
        bcftools annotate -a variant_summary.b37.vcf.gz -o $OUT_VCF -Oz $IN_VCF  -c +INFO
    else
        bcftools annotate -a variant_summary.b38.vcf.gz -o $OUT_VCF -Oz $IN_VCF  -c +INFO
    fi
    rm $IN_VCF*

    # test whether the reformatted HGMD VCF exists and wether its associated index doesn't
    # then we know whether to annotate the HGMD file or not
    if [ ! -f "$HGMD_PRO.reformated.vcf.gz.tbi" ] && [ -f "$HGMD_PRO.reformated.vcf.gz" ]; then

        IN_VCF=$OUT_VCF

        if [ ! -f $IN_VCF.tbi ]; then
            tabix -p vcf -f $IN_VCF
        fi

        OUT_VCF=${OUT_VCF%.vcf.gz}.HGMD.vcf.gz

        echo "Annotating HGMD"
        echo $OUT_VCF
        echo $IN_VCF

        bcftools annotate -a $HOME/HGMD_PRO.reformated.vcf.gz -o $OUT_VCF -Oz $IN_VCF  -c +INFO

        rm $IN_VCF*

    fi

    if [ ! -f $OUT_VCF.tbi ]; then
        tabix -p vcf -f $OUT_VCF
    fi

    VCF_UP=$(dx upload --brief $OUT_VCF)
    IDX_UP=$(dx upload --brief $OUT_VCF.tbi)

    dx-jobutil-add-output out_variants_vcfgzs "$VCF_UP" --class=array:file
    dx-jobutil-add-output out_variants_vcfgztbis "$IDX_UP" --class=array:file

    TO_RM=$(dx describe "$1" --name)

    rm $OUT_VCF*

}
export -f parallel_download_and_annotate


##########################
# Runs main code pipeline
# Arguments: None
# Returns: None
##########################
main() {

    export SHELL="/bin/bash"

    # download variant summary file 
    dx download "${variant_summary}" -o variant_summary.txt.gz

    # take security measures with temporary filenames 
    # for VCF file(s) and indices
    DXVCF_LIST=$(mktemp)
    DXIDX_LIST=$(mktemp)

    # reformat the resources in the variant summary file
    reformat_resources

    # Download all the VCF index files
    for i in "${!variants_vcfgztbis[@]}"; do
        
        echo "${variants_vcfgztbis[$i]}" >> "$DXIDX_LIST"
    
    done
    parallel -j $(nproc --all) -u --gnu parallel_download :::: $DXIDX_LIST ::: $HOME

    #Download and annotate VCF Files
    procCount=$(nproc --all)
    quarterCount=$(($procCount * 3 / 4))
    for i in "${!variants_vcfgzs[@]}"; do
      
        echo "${variants_vcfgzs[$i]}" >> $DXVCF_LIST
    
    done
    parallel -j $quarterCount -u --gnu parallel_download_and_annotate :::: $DXVCF_LIST ::: $HOME

}
