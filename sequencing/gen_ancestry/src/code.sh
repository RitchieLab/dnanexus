set -x -e -o pipefail

main() {

    echo "Value of pca_evec: '$pca_evec'"
    echo "Value of ancestry_txt: '$ancestry_txt'"
    echo "Value of thresh: '$thresh'"

	if test -z "$prefix"; then
		prefix=$(dx describe --name "$pca_evec" | sed 's/\.evec$//')
	fi

	WKDIR=$(mktemp -d)
	OUTDIR=$(mktemp -d)
	cd $WKDIR
	

    dx download "$pca_evec" -o evec
    dx download "$ancestry_txt" -o ancestry
    
    # remove "Other" from the ancestry; it's special!
    grep -v 'Other$' ancestry > ancestry_R
    # Format eigenvectors for loading into R
    tail -n+2 evec | sed -e 's/^  *//' -e 's/  *$//' -e 's/  */\t/g' -e 's/\t[^\t]*$//' > evec_R
    
    if test $(($num_evec + 1)) -gt "$(head -1 evec_R | tr '\t' '\n' | wc -l)"; then
    	dx-jobutil-report-error "ERROR: Desired number of eigenvectors exceeds the number of eigenvectors provided"
    fi
    
    # Assumes that "ancestry" and "evec" are formatted correctly.
    gen_ancestry.R ancestry_R evec_R $thresh $num_evec > $OUTDIR/$prefix.ancestry
    
    output_ancestry=$(dx upload $OUTDIR/$prefix.ancestry --brief)

    dx-jobutil-add-output output_ancestry "$output_ancestry" --class=file
}
