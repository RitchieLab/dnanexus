#!/bin/bash
set -ex -o pipefail

function slice_vcf () {
    set -x
    VCFFN=$(echo "$1" | cut -f3)
    PREFIX=$(echo "$1" | cut -f1)
    VCFIDXFN=$(echo "$1" | cut -f5)
    REGIONFN="$2"
    OUTDIR=$3
    SUFFIX="$4"
    if [ -n "$SUFFIX" ]; then
        SUFFIX="${SUFFIX}."
    fi
    OUTFN="$OUTDIR/$PREFIX.sliced.${SUFFIX}vcf.gz"

    download_intervals.py -H -f $VCFFN -i $VCFIDXFN -L "$REGIONFN" -o "$OUTFN"
    tabix -p vcf "$OUTFN"
}
export -f slice_vcf

main() {

    export SHELL=/bin/bash

    echo "Value of vcf_fn: '${vcf_fn[@]}'"
    echo "Value of vcfidx_fn: '${vcfidx_fn[@]}'"
    echo "Value of region_fn: '$region_fn'"
    echo "Value of suffix: '$suffix'"

    WKDIR=$(mktemp -d)
    cd $WKDIR
    OUTDIR=$(mktemp -d)

    VCF_NAMES=$(mktemp)
    VCFIDX_NAMES=$(mktemp)

    if test "${#vcf_fn[@]}" -ne "${#vcfidx_fn[@]}"; then
        dx-jobutil-report-error "ERROR: Number of VCF files and index files do not match"
    fi


    for i in ${!vcf_fn[@]}
    do
        echo -e "${vcf_fn_prefix[$i]}\t${vcf_fn_name[$i]}\t$(echo ${vcf_fn[$i]} | jq -r .["$dnanexus_link"])\t"
    done | tee $VCF_NAMES | wc -l

    # check the output
    less $VCF_NAMES

    for i in ${!vcfidx_fn[@]}
    do
        echo -e "${vcfidx_fn_prefix[$i]}\t${vcfidx_fn_name[$i]}\t$(echo ${vcfidx_fn[$i]} | jq -r .["$dnanexus_link"])\t"
    done | tee $VCFIDX_NAMES | wc -l
    
    # check the output
    less $VCFIDX_NAMES

    # OK, join the 2 files
    VCF_ALLINFO=$(mktemp)
    NUMCOMB=$(join -t$'\t' -j1 <(sort -t$'\t' -k1,1 $VCF_NAMES) <(sort -t$'\t' -k1,1 $VCFIDX_NAMES) | tee $VCF_ALLINFO | wc -l)

    # check the output
    less $VCF_ALLINFO

    if test "${#vcf_fn[@]}" -ne $NUMCOMB; then
        dx-jobutil-report-error "ERROR: Could not match VCF files to indexes"
    fi

    dx download "$region_fn" -o region.list

    # Run the slice_vcf function in parallel
    parallel --gnu -j $(nproc --all) slice_vcf :::: $VCF_ALLINFO ::: $WKDIR/region.list ::: $OUTDIR ::: $suffix

    # Upload all the results
    while read f; do

        vcf_out=$(dx upload "$f" --brief)
        dx-jobutil-add-output vcf_out "$vcf_out" --class=array:file
        vcfidx_out=$(dx upload "$f.tbi" --brief)
        dx-jobutil-add-output vcfidx_out "$vcfidx_out" --class=array:file

    done < <(ls -1 $OUTDIR/*.vcf.gz)
}
