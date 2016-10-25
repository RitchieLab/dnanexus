#!/bin/bash
# filter_bams 0.0.1
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

echo "Value of bam_files: '${bam_files[@]}'"
echo "Value of bai_files: '${bai_files[@]}'"
echo "Value of bed_file: '$bed_file'"

dx download "$bed_file" -o filter.bed

mkdir -p $HOME/out/filtered_bam_files

WKDIR=$(mktemp -d)
DXBAM_LIST=$(mktemp)
DXBAI_LIST=$(mktemp)

cd $WKDIR


function parallel_download() {
	#set -x
	cd $2
	dx download "$1"
	cd - >/dev/null
}
export -f parallel_download

function parallel_download_and_subset() {
	set -x
	cd $2
	samtools view -b -L filter.bed $1 > $HOME/out/filtered_bam_files/${1%.*}.g76.bam
	samtools index $HOME/out/filtered_bam_files/${1%.*}.g76.bam
	rm ${1%.*}.*

}


main() {

    cd $WKDIR

    # Download the BAM index files (in parallel)
    for i in "${!bai_files[@]}"; do
      echo "${bai_files[$i]}" >> $DXBAI_LIST
    done

    parallel -j $(nproc --all) -u --gnu parallel_download :::: $DXBAI_LIST ::: $WKDIR

    for i in "${!bam_files[@]}"; do
  		echo "${bam_files[$i]}" >> $DXBAM_LIST
  	done

    parallel -j $(nproc --all) -u --gnu parallel_download_and_subset :::: $DXBAM_LIST ::: $WKDIR

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    #dx download "$bed_file" -o bed_file
    #for i in ${!bam_files[@]}
    #do
    #    dx download "${bam_files[$i]}" -o bam_files-$i
    #done

    #for i in ${!bai_files[@]}
    #do
    #    dx download "${bai_files[$i]}" -o bai_files-$i
    #done

    # Fill in your application code here.
    #
    # To report any recognized errors in the correct format in
    # $HOME/job_error.json and exit this script, you can use the
    # dx-jobutil-report-error utility as follows:
    #
    #   dx-jobutil-report-error "My error message"
    #
    # Note however that this entire bash script is executed with -e
    # when running in the cloud, so any line which returns a nonzero
    # exit code will prematurely exit the script; if no error was
    # reported in the job_error.json file, then the failure reason
    # will be AppInternalError with a generic error message.

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    #dx-jobutil-add-output filtered_bam_files "$filtered_bam_files" --class=array:file

    dx-upload-all-outputs --parallel
}
