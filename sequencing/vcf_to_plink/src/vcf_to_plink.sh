#!/bin/bash

set -e -x -o pipefail

main() {
    echo "Value of vcf_fn: '${vcf_fn[@]}'"
    echo "Value of vcfidx_fn: '${vcfidx_fn[@]}'"
    echo "Value of excl_region: '${excl_region[@]}'"
    echo "Value of merge: '$merge'"
    echo "Value of snp_only: '$snp_only'"
    echo "Value of sel_args: '$sel_args'"
    echo "Value of merge_args: '$merge_args'"
    echo "Value of prefix: '$prefix'"

    # Sanity check - make sure vcf + vcfidx have same # of elements
    if test "${#vcfidx_fn[@]}" -ne "${#vcf_fn[@]}"; then
        dx-jobutil-report-error "ERROR: Number of VCFs and VCF indexes do NOT match!"
    fi

    # first, we need to match up the VCF and tabix index files
    # To do that, we'll get files of filename -> dxfile ID
    VCF_LIST=$(mktemp)
    for i in "${!vcf_fn[@]}"; do    
        dx describe --json "${vcf_fn[$i]}" | jq -r ".name,.id" | tr '\n' '\t' | sed 's/\t$/\n/' >> $VCF_LIST
    done
    
    VCFIDX_LIST=$(mktemp)
    for i in "${!vcfidx_fn[@]}"; do    
        dx describe --json "${vcfidx_fn[$i]}" | jq -r ".name,.id" | tr '\n' '\t' | sed -e 's/\t$/\n/' -e 's/\.tbi\t/\t/' >> $VCFIDX_LIST
    done
    
    # Now, get the prefix (strip off any .tbi) and join them
    JOINT_LIST=$(mktemp)
    join -t$'\t' -j1 <(sort -k1,1 $VCF_LIST) <(sort -k1,1 $VCFIDX_LIST) > $JOINT_LIST
        
    # Ensure that we still have the same number of files; throw an error if not
    if test $(cat $JOINT_LIST | wc -l) -ne "${#vcf_fn[@]}"; then
        dx-jobutil-report-error "ERROR: VCF files and indexes do not match!"
    fi
    
    SUBJOB_ARGS=""
    if test "$strict"; then
        SUBJOB_ARGS="$SUBJOB_ARGS -istrict:boolean=$strict"
    fi
    
    SUBJOB_ARGS="$SUBJOB_ARGS -isnp_only:boolean=$snp_only -ibiallelic:boolean=$biallelic -isel_args:string=\"$sel_args\""
    
    # and loop through the file, submitting sub-jobs
    while read VCF_LINE; do
        VCF_DXFN=$(echo "$VCF_LINE" | cut -f2)
        VCFIDX_DXFN=$(echo "$VCF_LINE" | cut -f3)
    
        SUBJOB=$(eval dx-jobutil-new-job convert_vcf "$SUBJOB_ARGS" -ivcf_fn:file="$VCF_DXFN" -ivcfidx_fn:file="$VCFIDX_DXFN" -iprefix_in:string="$prefix")
        
        dx-jobutil-add-output bed_out --array "$SUBJOB:bed" --class=jobref
        dx-jobutil-add-output bim_out --array "$SUBJOB:bim" --class=jobref
        dx-jobutil-add-output fam_out --array "$SUBJOB:fam" --class=jobref
        
    done < $JOINT_LIST

}

convert_vcf() {

    WKDIR=$(mktemp -d)
    OUTDIR=$(mktemp -d)
    
    cd $WKDIR
    # Now, get our VCF and VCF Idx file
    dx download "$vcf_fn" -o input.vcf.gz
    dx download "$vcfidx_fn" -o input.vcf.gz.tbi

    PREFIX=$(dx describe --name "$vcf_fn" | sed 's/.vcf\.gz$//')
    if test "$prefix_in"; then
        PREFIX="$prefix_in.$PREFIX"
    fi

    PLINK_ARGS=""
    
    if test "$snp_only" = "true"; then
        PLINK_ARGS="$PLINK_ARGS --snps-only"
    fi
    
    if test "$biallelic" = "true"; then
        PLINK_ARGS="$PLINK_ARGS --biallelic-only"
        if test "$strict" = "true"; then
            PLINK_ARGS="$PLINK_ARGS strict"
        fi
    fi
    
    # Now, convert the VCF into a PLINK file
    eval plink2 --vcf input.vcf.gz --double-id --id-delim "' '" --set-missing-var-ids @:#:\$1 --vcf-filter --make-bed "$sel_args" $PLINK_ARGS --out $OUTDIR/$PREFIX -allow-no-sex --allow-extra-chr --threads $(nproc --all)
    
    # upload all 3 bed/bim/fam files
    for ext in bed bim fam; do
        dxfn=$(dx upload --brief $OUTDIR/$PREFIX.$ext)
        dx-jobutil-add-output $ext $dxfn
    done
}

run_merge() {
    
    echo "Value of merge_args: '$merge_args'"
    echo "Value of fast_pca: '$fast_pca'"
    echo "Value of twstats: '$twstats'"
    echo "Value of num_evec: '$num_evec'"
    echo "Value of ldregress: '$ldregress'"
    echo "Value of numoutlier: '$numoutlier'"
    echo "Value of pca_opts: '$pca_opts'"
    
    WKDIR=$(mktemp -d)
    cd $WKDIR
    
    # First, we need to download all of the bed/bim/fam files
    # I'll assume that they are in order
    FAM_OVERALL=$(mktemp)
    N_F=0
    MAX_N=0
    FIRST_PREF=""
    MERGE_FILE=$(mktemp)
    
    for i in "${!bed[@]}"; do
        dx download "${bed[$i]}" -o f_$i.bed
        dx download "${bim[$i]}" -o f_$i.bim
        dx download "${fam[$i]}" -o f_$i.fam
        if test $N_F -eq 0; then
            FIRST_PREF="f_$i"
            MAX_N=$(cat f_$i.fam | wc -l)
            sed 's/[ \t][ \t]*/\t/g' f_$i.fam | cut -f1-2 | sort -t'\0' > $FAM_OVERALL
        else
            echo "f_$i" >> $MERGE_FILE
            TEST_N=$(cat f_$i.fam | wc -l)
            MAX_N=$((TEST_N > MAX_N ? TEST_N : MAX_N))
            TMPFAM=$(mktemp)
            join -t'\0' $FAM_OVERALL <(sed 's/[ \t][ \t]*/\t/g' f_$i.fam | cut -f1-2 | sort -t'\0') > $TMPFAM
            mv $TMPFAM $FAM_OVERALL            
        fi
        
        N_F=$((N_F + 1))
    done
    
    # Make sure the number of samples is identical for each
    if test "$(cat $FAM_OVERALL | wc -l)" -lt $MAX_N; then
        dx-jobutil-report-error "ERROR: Samples from parallel VCF conversions do not overlap!"
    fi
    
    OUTDIR=$(mktemp -d)
    
    # Allow some sample-level dropping to happen here (i.e. geno).
    eval plink2 --bfile "$FIRST_PREF" --merge-list $MERGE_FILE "$plink_args" --out $OUTDIR/$prefix --make-bed -allow-no-sex
    
    # get a list of those dropped
    join -v1 -t'\0' $FAM_OVERALL <(sed 's/[ \t][ \t]*/\t/g' $OUTDIR/$prefix.fam | cut -f1-2 | sort -t'\0') > $OUTDIR/$prefix.excluded
    
    # And upload results
    for ext in bed bim fam excluded; do
        dxfn=$(dx upload --brief $OUTDIR/$prefix.$ext)
        dx-jobutil-add-output $ext $dxfn
    done    
    
}
