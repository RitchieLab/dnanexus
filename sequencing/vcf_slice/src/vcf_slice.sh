#!/bin/bash
# vcf_slice 0.0.1
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.

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

    echo "Value of variants_vcfgzs: '${variants_vcfgzs[@]}'"
    echo "Value of variants_vcfgztbis: '${variants_vcfgztbis[@]}'"
    echo "Value of region_list: '$region_list'"
    echo "Value of suffix: '$suffix'"

    WKDIR=$(mktemp -d)
    cd $WKDIR
    OUTDIR=$(mktemp -d)

    VCF_NAMES=$(mktemp)
    VCFIDX_NAMES=$(mktemp)

    if test "${#variants_vcfgzs[@]}" -ne "${#variants_vcfgztbis[@]}"; then
        dx-jobutil-report-error "ERROR: Number of VCF files and index files do not match"
    fi


    for i in ${!variants_vcfgzs[@]}
    do
        echo -e "${variants_vcfgzs_prefix[$i]}\t${variants_vcfgzs_name[$i]}\t$(echo ${variants_vcfgzs[$i]} | jq -r .["$dnanexus_link"])"
    done | tee $VCF_NAMES | wc -l

    # check the output
    less $VCF_NAMES

    for i in ${!variants_vcfgztbis[@]}
    do
        echo -e "${variants_vcfgztbis_prefix[$i]}\t${variants_vcfgztbis_name[$i]}\t$(echo ${variants_vcfgztbis[$i]} | jq -r .["$dnanexus_link"])"
    done | tee $VCFIDX_NAMES | wc -l
    
    # check the output
    less $VCFIDX_NAMES

    # OK, join the 2 files
    VCF_ALLINFO=$(mktemp)
    NUMCOMB=$(join -t$'\t' -j1 <(sort -t$'\t' -k1,1 $VCF_NAMES) <(sort -t$'\t' -k1,1 $VCFIDX_NAMES) | tee $VCF_ALLINFO | wc -l)

    # check the output
    less $VCF_ALLINFO

    if test "${#variants_vcfgzs[@]}" -ne $NUMCOMB; then
        dx-jobutil-report-error "ERROR: Could not match VCF files to indices."
    fi

    dx download "$region_list" -o region.list

    # Check that the region.list is in the correct format
    format_check=$(head region.list -n 1 | sed -e '/^[a-z0-9]*:[0-9]*-[0-9]*$/d')
    if ! [[ -z $format_check ]]; then
        dx-jobutil-report-error "ERROR: Regions list not in proper format. Regions should be listed: chrom:123-456"
    fi

    # Run the slice_vcf function in parallel
    parallel --gnu -j $(nproc --all) slice_vcf :::: $VCF_ALLINFO ::: $WKDIR/region.list ::: $OUTDIR ::: $suffix

    # Upload all the results
    while read f; do

        out_variants_vcfgz=$(dx upload "$f" --brief)
        dx-jobutil-add-output out_variants_vcfgz "$out_variants_vcfgz" --class=array:file
        out_variants_vcfgztbi=$(dx upload "$f.tbi" --brief)
        dx-jobutil-add-output out_variants_vcfgztbi "$out_variants_vcfgztbi" --class=array:file

    done < <(ls -1 $OUTDIR/*.vcf.gz)
}
