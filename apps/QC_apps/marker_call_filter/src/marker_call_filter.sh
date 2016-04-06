#!/bin/bash
# sample_call_filter 0.0.1
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

main() {

    echo "Value of bed_file: '$bed_file'"
    echo "Value of bim_file: '$bim_file'"
    echo "Value of fam_file: '$fam_file'"
    echo "Value of output_prefix: '$output_prefix'"
    echo "Value of threshold: '$threshold'"
    echo "Value of old_plink: '$old_plink'"

    PLINK_CMD="plink2"
    if test "$old_plink" == "true"; then
        PLINK_CMD="plink"
    fi

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    dx download "$bed_file" -o bed_file

    dx download "$bim_file" -o bim_file

    dx download "$fam_file" -o fam_file

    # Fill in your application code here.
	if [ -n $output_prefix ]
	then
		OUTPUT_PREFIX="$output_prefix"
	else
		OUTPUT_PREFIX=geno_filter
	fi
	$PLINK_CMD --bed bed_file --bim bim_file --fam fam_file --noweb --allow-no-sex --mind 1 --geno 1 --out $OUTPUT_PREFIX.orig --write-snplist
	$PLINK_CMD --bed bed_file --bim bim_file --fam fam_file --noweb --allow-no-sex --mind 1 --geno $threshold --out $OUTPUT_PREFIX.filtered --write-snplist

	join -v 1 <(sort $OUTPUT_PREFIX.orig.snplist) <(sort $OUTPUT_PREFIX.filtered.snplist) > $OUTPUT_PREFIX.exclude

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

    # The following line(s) use the dx command-line tool to upload your file
    # outputs after you have created them on the local file system.  It assumes
    # that you have used the output field name for the filename for each output,
    # but you can change that behavior to suit your needs.  Run "dx upload -h"
    # to see more options to set metadata.

    drop_list=$(dx upload $OUTPUT_PREFIX.exclude --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output drop_list "$drop_list" --class=file
}