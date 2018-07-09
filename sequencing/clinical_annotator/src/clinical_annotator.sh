#!/bin/bash
# clinical_annotator 0.0.1
# Generated by dx-app-wizard.
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
export SHELL="/bin/bash"

set -x
# install GNU parallel!
sudo sed -i 's/^# *\(deb .*backports.*\)$/\1/' /etc/apt/sources.list
sudo apt-get update
sudo apt-get install --yes parallel

cd $HOME
wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
tar -xjf bcftools-1.8.tar.bz2
cd /home/dnanexus/bcftools-1.8
./configure --prefix=/usr/local/
make -s
make install

function parallel_download() {
	cd $2
	dx download "$1"
	cd - >/dev/null
}
export -f parallel_download

#Downloads and reformats resource files
function download_resources(){

	cd $HOME

	wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

	python $HOME/ClinVar_tsvTOvcf.py variant_summary.txt.gz

	vcf-sort -c variant_summary.b37.vcf | bgzip -c > variant_summary.b37.vcf.gz
	vcf-sort -c variant_summary.b38.vcf | bgzip -c > variant_summary.b38.vcf.gz
	tabix -p vcf -f variant_summary.b37.vcf.gz
	tabix -p vcf -f variant_summary.b38.vcf.gz

	if test "$hgmd_pro_file"; then
		dx download "$hgmd_pro_file"
		python $HOME/reformatHGMD.py $(dx describe "$hgmd_pro_file" --name) | bgzip -c > $HOME/HGMD_PRO.reformated.vcf.gz
		tabix -p vcf -f HGMD_PRO.reformated.vcf.gz
	fi

}
export -f download_resources

#Downloads input vcf files and Annotates
function parallel_download_and_annotate() {
	set -x
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

	# HGMD
		if [ ! -f $HGMD_PRO.reformated.vcf.gz.tbi ]; then

			IN_VCF=$OUT_VCF

			if [ ! -f $IN_VCF.tbi ]; then
				tabix -p vcf -f $IN_VCF
			fi

			OUT_VCF=${OUT_VCF%.vcf.gz}.HGMD.vcf.gz

			echo "HGMD"
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

	dx-jobutil-add-output vcf_out "$VCF_UP" --class=array:file
	dx-jobutil-add-output vcfidx_out "$IDX_UP" --class=array:file

	TO_RM=$(dx describe "$1" --name)

	rm $OUT_VCF*


}
export -f parallel_download_and_annotate

main() {

	export SHELL="/bin/bash"

	DXVCF_LIST=$(mktemp)
	DXIDX_LIST=$(mktemp)

	download_resources

  # Download all the VCF index files
	  for i in "${!vcfidx_fn[@]}"; do
	    echo "${vcfidx_fn[$i]}" >> "$DXIDX_LIST"
	  done

  	parallel -j $(nproc --all) -u --gnu parallel_download :::: $DXIDX_LIST ::: $HOME

	#Download and annotate VCF Files
		procCount=$(nproc --all)
		quarterCount=$(($procCount * 3 / 4))

    for i in "${!vcf_fn[@]}"; do
      echo "${vcf_fn[$i]}" >> $DXVCF_LIST
    done

    parallel -j $quarterCount -u --gnu parallel_download_and_annotate :::: $DXVCF_LIST ::: $HOME


}
