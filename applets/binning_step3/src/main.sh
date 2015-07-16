#!/bin/bash

set -e -x -o pipefail


analyze_phenotype() {
	P=$1
	PCOL=$(echo "$P + 1" | bc)
	LCOL=$(echo "$P + 4" | bc)
	PHENO="$( \
		head -n 1 biobin/output-locus.csv \
		| cut -d ',' -f $LCOL \
		| sed -r 's/ Bin Name\(s\)$//g' \
	)"
	
	# symlink shared input files for SKAT
	mkdir -p skat/$P
	cd skat/$P
	ln -s ../input.covariate input.covariate
	ln -s ../../plink/output.bed input.bed
	ln -s ../../plink/output.bim input.bim
	
	# splice the phenotype column into the FAM file
	sed 's/,/\t/g' ../../biobin/output-$PHENO-bins.csv \
	| gawk \
		'{ if (FNR==NR) { p[$1]=$2 } else { if ($2 in p) $6=p[$2]; else $6=-9; print } }' \
		- \
		../../plink/output.fam \
	> input.fam
	
	# prepare Set ID file
	cat ../../biobin/output-locus.csv \
	| cut -d ',' -f 3,$LCOL \
	| sed '1d' \
	| sed -r 's/^([^,-]+)(-?[^,-]*)[^,]*,(.*)$/\1\2,\3/g' \
	| sed 's/chr//g' \
	| sed 's/-/:/g' \
	| gawk '!x[$0]++' \
	> temp-SetID.txt
	
	/usr/bin/SKAT_magic.py2 temp-SetID.txt \
	| grep '. ' \
	| sort -k1,1 \
	> input-SetID.txt
	
	# run SKAT
	if [[ "$regression_type" == "logistic" ]]; then
		SKAT_ARG="D"
	else
		SKAT_ARG="C"
	fi
	Rscript /usr/share/runSKAT.R \
		$SKAT_ARG \
	2>&1 | tee output.log
	cd ../..
	ls -laR skat/$P
	
	# run PLATO
	COVARIATES="$(head -n 1 input/input.covariate | sed -r 's/[ \t]+/,/g' | cut -d ',' -f 2-)"
	mkdir -p plato/$P
	cd plato/$P
	plato \
		load-trait \
			--ignore-error \
			--dummy-samples \
			--no-fid \
			--missing "nan" \
			--file <( sed '2,9 d' ../../biobin/output-$PHENO-bins.csv | tr ',' ' ') \
		load-trait \
			--ignore-error \
			--extra-samples \
			--no-fid \
			--missing "$missing_code" \
			--file ../../input/input.covariate \
		"$regression_type" \
			--outcome Status \
			--covariates "$COVARIATES" \
			--exclude-markers \
			--use-traits \
			--output output-linear.txt \
	2>&1 | tee output.log
	cd ../..
	ls -laR plato/$P
	
	# move files into position for upload
	mv skat/$P/temp.SSD_LOG.txt "$HOME/out/skat_files/temp-$P-SSD_LOG.txt"
	mv skat/$P/output-unweighted-uncorrected.txt "$HOME/out/skat_files/output-$P-unweighted-uncorrected.txt"
	mv skat/$P/output-weighted-uncorrected.txt "$HOME/out/skat_files/output-$P-weighted-uncorrected.txt"
	mv skat/$P/output-defaultweighted-uncorrected.txt "$HOME/out/skat_files/output-$P-defaultweighted-uncorrected.txt"
	mv plato/$P/output-linear.txt "$HOME/out/plato_files/output-linear-$P.txt"
}


main() {
	# install GNU parallel
	
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	sudo sed -i 's/^# *\(deb .*backports.*\)$/\1/' /etc/apt/sources.list 
	sudo apt-get update --yes
	sudo apt-get install --yes parallel
	
	
	# install SKAT-unlimited
	R CMD INSTALL /usr/share/SKAT_unlimited.tgz
	
	
	# switch to a temp directory and download all input files
	
	TMPDIR="$(mktemp -d)"
	cd "$TMPDIR"
	mkdir input
	dx download "$covariate_file" -o input/input.covariate
	ls -laR input
	
	
	# extract checkpoint archives
	dx download "$checkpoint_biobin" -o biobin.tar.bz
	tar xvf biobin.tar.bz 2>&1 | tee -a output.log
	ls -laR biobin
	dx download "$checkpoint_vcftools" -o vcftools.tar.bz
	tar xvf vcftools.tar.bz 2>&1 | tee -a output.log
	ls -laR vcftools
	dx download "$checkpoint_plink" -o plink.tar.bz
	tar xvf plink.tar.bz 2>&1 | tee -a output.log
	ls -laR plink
	
	
	# convert BIM file to add chromosome to unlabeled variants
	mv plink/output.bim plink/original.bim
	gawk \
		'{ if ($2==$4) $2=("" $1 ":" $2); print }' \
		plink/original.bim \
	> plink/output.bim
	
	
	# convert covariate file for SKAT
	mkdir skat
	gawk \
		'{ if (NR==1) { print "FID " $0 } else { print $1 " " $0 } }' \
		input/input.covariate \
	> skat/input.covariate
	
	
	# create output directories
	mkdir -p "$HOME/out/debug_file"
	mkdir -p "$HOME/out/skat_files"
	mkdir -p "$HOME/out/plato_files"
	
	
	# analyze each phenotype in parallel
	NUM_PHENO=$(head -n 1 biobin/output-locus.csv | gawk -F ',' '{ print (NF - 4) }')
	export -f analyze_phenotype
	export SHELL="/bin/bash"
	seq 1 $NUM_PHENO | parallel --gnu --max-procs $(nproc --all) analyze_phenotype
	
	# collect each pheotype's SKAT and plato logs
	for P in $(seq 1 $NUM_PHENO); do
		cat skat/$P/output.log >> output.log
		cat plato/$P/output.log >> output.log
	done
	
	
	# move files into position for upload
	tar cjvf "$HOME/out/debug_file/debug.tar.bz" \
		plink/original.bim \
		plink/output.bim \
		skat/input.covariate \
		skat/*/input.fam \
		skat/*/temp-SetID.txt \
		skat/*/input-SetID.txt
	mkdir -p "$HOME/out/log_file"
	mv output.log "$HOME/out/log_file"
	
	
	# make sure we don't fail just for missing output files
	if [[ $(ls -1 "$HOME/out/skat_files" | wc -l) -lt 1 ]]; then
		touch "$HOME/out/skat_files/NO_SKAT_FILES"
	fi
	if [[ $(ls -1 "$HOME/out/plato_files" | wc -l) -lt 1 ]]; then
		touch "$HOME/out/plato_files/NO_PLATO_FILES"
	fi
	
	
	# return to the home dir and upload all files
	cd "$HOME"
	dx-upload-all-outputs
}
