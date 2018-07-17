#!/bin/bash
# vcf_pca 0.0.1
# Generated by dx-app-wizard.
#
# Parallelized execution pattern: Your app will generate multiple jobs
# to perform some computation in parallel, followed by a final
# "postprocess" stage that will perform any additional computations as
# necessary.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() or any other entry point is ALWAYS
# executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.

set -x -e -o pipefail

#echo "deb http://us.archive.ubuntu.com/ubuntu xenial main restricted universe multiverse " >> /etc/apt/sources.list

if [[ "$DX_RESOURCES_ID" != "" ]]; then
  DX_ASSETS_ID="$DX_RESOURCES_ID"
else
  DX_ASSETS_ID="$DX_PROJECT_CONTEXT_ID"
fi

main() {

    echo "Value of variants_vcfgz: '${variants_vcfgz[@]}'"
    echo "Value of variants_vcfgztbi: '${variants_vcfgztbi[@]}'"
    echo "Value of bed_fn: '${bed_fn[@]}'"
    echo "Value of bim_fn: '${bim_fn[@]}'"
    echo "Value of fam_fn: '${fam_fn[@]}'"

    echo "Value of excl_region: '${excl_region[@]}'"
    echo "Value of sel_args: '$sel_args'"
    echo "Value of maf: '$maf'"
    echo "Value of ld_args: '$ld_args'"
    echo "Value of merge_args: '$merge_args'"
    echo "Value of fast_pca: '$fast_pca'"
    echo "Value of twstats: '$twstats'"
    echo "Value of num_evec: '$num_evec'"
    echo "Value of ldregress: '$ldregress'"
    echo "Value of numoutlier: '$numoutlier'"
    echo "Value of pca_opts: '$pca_opts'"
    echo "value of DX_ASSETS_ID: $DX_ASSETS_ID"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

	# Sanity checks:

	USE_VCF=1
	if test -z "${variants_vcfgz[@]}" -o -z "${variants_vcfgztbi[@]}"; then
		# If here, either VCF/TBI is empty
		USE_VCF=0
		if test -z "${bed_fn[@]}" -o -z "${bim_fn[@]}" -o -z "${fam_fn[@]}" ; then
			dx-jobutil-report-error "ERROR: You must provide VCF/TBI files or BED/BIM/FAM files"
		fi
	fi

	PCA_ARGS="-inum_evec:int=$num_evec -imerge_args:string=\"$merge_args\" -ifast_pca:boolean=$fast_pca -itwstats:boolean=$twstats -ildregress:int=$ldregress -inumoutlier:int=$numoutlier -ipca_opts:string=\"$pca_opts\" -iprefix:string=$prefix"
	SUBJOB_ARGS="-imaf:float=\"$maf\" -ild_args:string=\"$ld_args\" -isel_args:string=\"$sel_args\""

	if test $USE_VCF -gt 0; then
		# - make sure vcf + vcfidx have same # of elements
		if test "${#variants_vcfgztbi[@]}" -ne "${#variants_vcfgz[@]}"; then
			dx-jobutil-report-error "ERROR: Number of VCFs and VCF indexes do NOT match!"
		fi

		# first, we need to match up the VCF and tabix index files
		# To do that, we'll get files of filename -> dxfile ID
		VCF_LIST=$(mktemp)
		for i in "${!variants_vcfgz[@]}"; do
			dx describe --json "${variants_vcfgz[$i]}" | jq -r ".name,.id" | tr '\n' '\t' | sed 's/\t$/\n/' >> $VCF_LIST
		done

		VCFIDX_LIST=$(mktemp)
		for i in "${!variants_vcfgztbi[@]}"; do
			dx describe --json "${variants_vcfgztbi[$i]}" | jq -r ".name,.id" | tr '\n' '\t' | sed -e 's/\t$/\n/' -e 's/\.tbi\t/\t/' >> $VCFIDX_LIST
		done

		# Now, get the prefix (strip off any .tbi) and join them
		JOINT_LIST=$(mktemp)
		join -t$'\t' -j1 <(sort -k1,1 $VCF_LIST) <(sort -k1,1 $VCFIDX_LIST) > $JOINT_LIST

		# Ensure that we still have the same number of files; throw an error if not
		if test $(cat $JOINT_LIST | wc -l) -ne "${#variants_vcfgz[@]}"; then
			dx-jobutil-report-error "ERROR: VCF files and indexes do not match!"
		fi


		if test "${excl_region[@]}"; then
			for i in "${!excl_region[@]}"; do
				SUBJOB_ARGS="$SUBJOB_ARGS -iexcl_region:array:file=$(dx describe --json "${excl_region[$i]}" | jq -r .id)"
			done
		fi

		# and loop through the file, submitting sub-jobs
		while read VCF_LINE; do
			VCF_DXFN=$(echo "$VCF_LINE" | cut -f2)
			VCFIDX_DXFN=$(echo "$VCF_LINE" | cut -f3)

			SUBJOB=$(eval dx-jobutil-new-job downsample_vcf "$SUBJOB_ARGS" -ivcf_fn:file="$VCF_DXFN" -ivcfidx_fn:file="$VCFIDX_DXFN" -iproject_1kg:boolean="$project_1kg")

			PCA_ARGS="$PCA_ARGS -ibed:array:file=${SUBJOB}:bed -ibim:array:file=${SUBJOB}:bim -ifam:array:file=${SUBJOB}:fam"

		done < $JOINT_LIST

	else
		# We're using bed/bim/fam here!

		# - make sure vcf + vcfidx have same # of elements
		if test "${#bed_fn[@]}" -ne "${#bim_fn[@]}" -o "${#bed_fn[@]}" -ne "${#fam_fn[@]}"; then
			dx-jobutil-report-error "ERROR: Number of BED/BIM/FAM files do NOT match!"
		fi

		# first, we need to match up the VCF and tabix index files
		# To do that, we'll get files of filename -> dxfile ID
		BED_LIST=$(mktemp)
		for i in "${!bed_fn[@]}"; do
			dx describe --json "${bed_fn[$i]}" | jq -r ".name,.id" | tr '\n' '\t' | sed -e 's/\t$/\n/' -e 's/\.bed\t/\t/' >> $BED_LIST
		done

		BIM_LIST=$(mktemp)
		for i in "${!bim_fn[@]}"; do
			dx describe --json "${bim_fn[$i]}" | jq -r ".name,.id" | tr '\n' '\t' | sed -e 's/\t$/\n/' -e 's/\.bim\t/\t/' >> $BIM_LIST
		done

		FAM_LIST=$(mktemp)
		for i in "${!fam_fn[@]}"; do
			dx describe --json "${fam_fn[$i]}" | jq -r ".name,.id" | tr '\n' '\t' | sed -e 's/\t$/\n/' -e 's/\.fam\t/\t/' >> $FAM_LIST
		done


		# Now, get the prefix (strip off any .tbi) and join them
		JOINT_LIST=$(mktemp)
		join -t$'\t' -j1 <(sort -k1,1 $BED_LIST) <(join -t$'\t' -j1 <(sort -k1,1 $BIM_LIST) <(sort -k1,1 $FAM_LIST) ) > $JOINT_LIST

		# Ensure that we still have the same number of files; throw an error if not
		if test $(cat $JOINT_LIST | wc -l) -ne "${#bed_fn[@]}"; then
			dx-jobutil-report-error "ERROR: BED/BIM/FAM files do not match!"
		fi

		while read PLINK_LINE; do
			BED_DXFN=$(echo "$PLINK_LINE" | cut -f2)
			BIM_DXFN=$(echo "$PLINK_LINE" | cut -f3)
			FAM_DXFN=$(echo "$PLINK_LINE" | cut -f4)

			SUBJOB=$(eval dx-jobutil-new-job downsample_plink "$SUBJOB_ARGS" -ibed_fn:file="$BED_DXFN" -ibim_fn:file="$BIM_DXFN" -ifam_fn:file="$FAM_DXFN" -iproject_1kg:boolean="$project_1kg")

			PCA_ARGS="$PCA_ARGS -ibed:array:file=${SUBJOB}:bed -ibim:array:file=${SUBJOB}:bim -ifam:array:file=${SUBJOB}:fam"

		done < $JOINT_LIST

	fi

    pcarun=$(eval dx-jobutil-new-job run_pca -iproject_1kg:boolean="$project_1kg" -iproject_superpop:boolean="$project_superpop" "$PCA_ARGS")

    dx-jobutil-add-output evec_out "$pcarun:evec_out" --class=jobref
	dx-jobutil-add-output samp_excl "$pcarun:samp_excl" --class=jobref
    if test "$fast_pca" = "false"; then
	    dx-jobutil-add-output eval_out "$pcarun:eval_out" --class=jobref
	    if test "$twstats" = "true"; then
		    dx-jobutil-add-output twstats_out "$pcarun:twstats_out" --class=jobref
		fi
	fi
}

function run_ld(){
	OUTDIR="$1"
	PREFIX="$2"
	INBASE="$3"
	ld_args="$4"

	# And LD prune (if needed)
	if test "$ld_args"; then
		eval plink2 --bfile $INBASE "$ld_args" --threads $(nproc --all) -allow-no-sex --out ld_list
		plink2 --bfile $INBASE --extract ld_list.prune.in --make-bed --out $OUTDIR/$PREFIX
	else
		for ext in bed bim fam; do
			mv $INBASE.$ext $OUTDIR/$PREFIX.$ext
		done
	fi
}


downsample_plink(){
	WKDIR=$(mktemp -d)
	OUTDIR=$(mktemp -d)

	cd $WKDIR
	# Now, get our VCF and VCF Idx file
	dx download "$bed_fn" -o input.bed
	dx download "$bim_fn" -o input.bim
	dx download "$fam_fn" -o input.fam

	eval plink2 --bfile input --maf $maf "$sel_args" --out preld -allow-no-sex --threads $(nproc --all) --make-bed || touch preld.bed preld.bim preld.fam

	PREFIX=$(dx describe --name "$bed_fn" | sed 's/.bed$//')

	if test -s preld.bed; then
		PRELD_NAME="preld"

		if test "$project_1kg" = "true"; then
			NCHR=0
			PLINK_FN=""
			MERGE_FILE=$(mktemp)
			SNP_LIST=$(mktemp)
			GEN_DIR=$(mktemp -d)

			cut -f2 preld.bim > $SNP_LIST

            sed -i '/rs10761581\|rs11204210\|rs3167875\|rs35515471\|rs3814160\|rs10902343\|rs11246606\|rs11246607\|rs11246608\|rs34285763\|rs61890334\|rs61890420\|rs61890422\|rs200029677\|rs201794505\|rs7170838\|rs367595809\|rs10910824\|rs11590105\|rs1628172\|rs2274616\|rs28391411\|rs28602496\|rs28639473\|rs28676508\|rs41302235\|rs41315701\|rs61742539\|rs6670984\|rs142831593\|rs199675524\|rs201703456\|rs201906478\|rs202089732\|rs62561229\|rs78487056/d' $SNP_LIST

			for c in $(sed 's/  */\t/g' preld.bim | cut -f1 | sort -u); do
				# download the bed/bim/fam from dnanexus
				for ext in bed bim fam; do
                    dx download "$DX_ASSETS_ID:/1K_genomes/b38/ALL.chr${c}_GRCh38.genotypes.20170504.genotypes.$ext" -o $GEN_DIR/ALL.chr$c.$ext
                    #ALL.chr$c.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz
                    done
				# extract the markers in preld
				plink2 --bfile $GEN_DIR/ALL.chr$c --extract $SNP_LIST --out $GEN_DIR/ALL.chr$c.extracted --make-bed --allow-no-sex

				# add the bed/bim/fam to the MERGE_FILE
				echo -e "$GEN_DIR/ALL.chr$c.extracted.bed\t$GEN_DIR/ALL.chr$c.extracted.bim\t$GEN_DIR/ALL.chr$c.extracted.fam" >> $MERGE_FILE
			done

			# now, extract the markers overlapping the 1kg data
			GEN_SNPS=$(mktemp)
			for f in $(cut -f2 $MERGE_FILE); do
				cut -f2 $f >> $GEN_SNPS
			done

			plink2 --bfile preld --extract $GEN_SNPS --out premerge --make-bed --allow-no-sex

			PRELD_NAME="postmerge"
			# and finally, merge the files together
			plink2 --bfile premerge --merge-list $MERGE_FILE --out $PRELD_NAME --allow-no-sex --make-bed

		fi

		run_ld "$OUTDIR" "$PREFIX" $PRELD_NAME "$ld_args"
	else
		for ext in bed bim fam; do
			mv preld.$ext $OUTDIR/$PREFIX.$ext
		done
	fi

	# upload all 3 bed/bim/fam files
	for ext in bed bim fam; do
		dxfn=$(dx upload --brief $OUTDIR/$PREFIX.$ext)
	    dx-jobutil-add-output $ext $dxfn
	done
}


downsample_vcf() {
    # Fill in your process code here

    WKDIR=$(mktemp -d)
	OUTDIR=$(mktemp -d)

	cd $WKDIR
	# Now, get our VCF and VCF Idx file
	dx download "$vcf_fn" -o input.vcf.gz
	dx download "$vcfidx_fn" -o input.vcf.gz.tbi

	PREFIX=$(dx describe --name "$vcf_fn" | sed 's/.vcf\.gz$//')

	# Now, convert the VCF into a PLINK file
	eval plink2 --vcf input.vcf.gz --double-id --id-delim "' '" --vcf-filter --biallelic-only strict --snps-only --set-missing-var-ids @:#:\$1 --make-bed --maf $maf "$sel_args" --out preld -allow-no-sex --threads $(nproc --all) || touch preld.bed preld.bim preld.fam

	if test -s preld.bed; then
		# if we want to project onto 1Kg, let's download only the necessary chromosomes, but make sure only autosomes!

		PRELD_NAME="preld"

		if test "$project_1kg" = "true"; then
			NCHR=0
			PLINK_FN=""
			MERGE_FILE=$(mktemp)
			SNP_LIST=$(mktemp)
			GEN_DIR=$(mktemp -d)

			cut -f2 preld.bim > $SNP_LIST
			for c in $(join <(seq 1 22 | sort) <(tabix -l input.vcf.gz | sort)); do
				# download the bed/bim/fam from dnanexus
				for ext in bed bim fam; do
					dx download "$DX_ASSETS_ID:/1K_genomes/b38/ALL.chr${c}_GRCh38.genotypes.20170504.genotypes.$ext" -o $GEN_DIR/ALL.chr$c.$ext

				done
				# extract the markers in preld
				plink2 --bfile $GEN_DIR/ALL.chr$c --extract $SNP_LIST --out $GEN_DIR/ALL.chr$c.extracted --make-bed --allow-no-sex

				# add the bed/bim/fam to the MERGE_FILE
				echo -e "$GEN_DIR/ALL.chr$c.extracted.bed\t$GEN_DIR/ALL.chr$c.extracted.bim\t$GEN_DIR/ALL.chr$c.extracted.fam" >> $MERGE_FILE
			done

			# now, extract the markers overlapping the 1kg data
			GEN_SNPS=$(mktemp)
			for f in $(cut -f2 $MERGE_FILE); do
				cut -f2 $f >> $GEN_SNPS
			done

			plink2 --bfile preld --extract $GEN_SNPS --out premerge --make-bed --allow-no-sex

			PRELD_NAME="postmerge"
			# and finally, merge the files together
			plink2 --bfile premerge --merge-list $MERGE_FILE --out $PRELD_NAME --allow-no-sex --make-bed

		fi

		run_ld "$OUTDIR" "$PREFIX" $PRELD_NAME "$ld_args"
	else
		for ext in bed bim fam; do
			mv preld.$ext $OUTDIR/$PREFIX.$ext
		done
	fi

	# upload all 3 bed/bim/fam files
	for ext in bed bim fam; do
		dxfn=$(dx upload --brief $OUTDIR/$PREFIX.$ext)
	    dx-jobutil-add-output $ext $dxfn
	done
}

run_pca() {
    # Fill in your postprocess code here

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

		if test -s f_$i.bed; then
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
		fi
	done

	# Make sure the number of samples is identical for each
	if test "$(cat $FAM_OVERALL | wc -l)" -lt $MAX_N; then
		dx-jobutil-report-error "ERROR: Samples from parallel VCF conversions do not overlap!"
	fi

	INPUTDIR=$(mktemp -d)

	# Allow some sample-level dropping to happen here (i.e. geno).
	eval plink2 --bfile "$FIRST_PREF" --merge-list $MERGE_FILE "$merge_args" --out $INPUTDIR/input --make-bed -allow-no-sex

	# get a list of those dropped
	SAMPLE_DROPPED=$(mktemp)
	join -v1 -t'\0' $FAM_OVERALL <(sed 's/[ \t][ \t]*/\t/g' $INPUTDIR/input.fam | cut -f1-2 | sort -t'\0') > $SAMPLE_DROPPED

	PAR_F=$(mktemp)

	# OK, now if projecting, we need to add populations
	if test "$project_1kg" = "true"; then
		FAM_LINENO=$(mktemp)
		nl -ba -nln -w1 $INPUTDIR/input.fam | sed -e 's/  */\t/g' > $FAM_LINENO
		POP_FILE=$(mktemp)
		dx download "$DX_ASSETS_ID:/1K_genomes/integrated_call_samples_v3.20130502.ALL.panel" -o $POP_FILE -f

		POP_COL=2
		if test "$project_superpop" = "true"; then
			POP_COL=3
		fi

		join -t$'\t' -1 3 -2 1 -a 1 <(sort -t$'\t' -k3,3 $FAM_LINENO | cut -f1-6) <(cut -f1,$POP_COL $POP_FILE | sort -t$'\t' -k1,1) \
			| sort -k2,2n -t$'\t' | cut -f1,3- \
			| awk '{if (NF<6) print $0 "\tUNK"; else print $0;}' > $INPUTDIR/input.fam.new

		mv $INPUTDIR/input.fam $INPUTDIR/input.fam.old
		mv $INPUTDIR/input.fam.new $INPUTDIR/input.fam

		POPLIST=$(mktemp)
		cut -f $POP_COL $POP_FILE | tail -n+2 | sort -u > $POPLIST
		echo "poplist: $POPLIST" >> $PAR_F

	else
		# remove the "-9" in the last column of the PED file
		sed -i 's/-9$/UNK/' $INPUTDIR/input.fam

	fi

	OUTDIR=$(mktemp -d)

	echo "genotypename: $INPUTDIR/input.bed" >> $PAR_F
	echo "snpname: $INPUTDIR/input.bim" >> $PAR_F
	echo "indivname: $INPUTDIR/input.fam" >> $PAR_F
	echo "evecoutname: $OUTDIR/$prefix.evec" >> $PAR_F
  #print PARFILE ("evaloutname: $e\n");
  #print PARFILE ("altnormstyle: NO\n");
	echo "numoutevec: $num_evec" >> $PAR_F
  #print PARFILE ("numoutlieriter: $m\n");
  #print PARFILE ("numoutlierevec: $t\n");
	echo "ldregress: $ldregress" >> $PAR_F
	echo "familynames: NO" >> $PAR_F

    if test "$fast_pca" = "false"; then
		echo "evaloutname: $OUTDIR/$prefix.eval" >> $PAR_F
		echo "numthreads: $(nproc --all)" >> $PAR_F
		echo "numoutlieriter: $numoutlier" >> $PAR_F
		echo "outlieroutname: $OUTDIR/$prefix.outlier" >> $PAR_F
		echo "numoutlierevec: $num_evec" >> $PAR_F
	else
		echo "fastmode: YES" >> $PAR_F
	fi

	echo "$pca_opts" | sed 's/"//g' | tr ',' '\n' >> $PAR_F

	ulimit -c unlimited

	smartpca -p $PAR_F > $OUTDIR/$prefix.eval

	# And upload results
	evec_dxfn=$(dx upload --brief $OUTDIR/$prefix.evec)
    dx-jobutil-add-output evec_out "$evec_dxfn" --class=file

    # concatenate the sample exluded lists
    cat $SAMPLE_DROPPED <(awk '{print $3 "\t" $3}' $OUTDIR/$prefix.outlier) > $OUTDIR/$prefix.excluded
    se_dxfn=$(dx upload --brief $OUTDIR/$prefix.excluded)
	dx-jobutil-add-output samp_excl "$se_dxfn" --class=file


    if test "$fast_pca" = "false"; then
    	eval_dxfn=$(dx upload --brief $OUTDIR/$prefix.eval)
	    dx-jobutil-add-output eval_out "$eval_dxfn" --class=file
	    if test "$twstats" = "true"; then
	    	# do the twstats calculation here
	    	twstats -t /usr/share/twtable -i $OUTDIR/$prefix.eval -o $OUTDIR/$prefix.twstats
	    	tw_dxfn=$(dx upload --brief $OUTDIR/$prefix.twstats)
		    dx-jobutil-add-output twstats_out "$tw_dxfn" --class=file
		fi
	fi


}
