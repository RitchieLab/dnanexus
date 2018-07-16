set -x

main() {

    echo "Value of evec_fn: '$evec_fn'"
    echo "Value of ancestry_fn: '$ancestry_fn'"
    echo "Value of thresh: '$thresh'"

	if test -z "$prefix"; then
		prefix=$(dx describe --name "$evec_fn" | sed 's/\.evec$//')
	fi

	WKDIR=$(mktemp -d)
	OUTDIR=$(mktemp -d)
	cd $WKDIR
	

    dx download "$evec_fn" -o evec
    dx download "$ancestry_fn" -o ancestry
    
    # remove "Other" from the ancestry; it's special!
    grep -v 'Other$' ancestry > ancestry_R
    # Format eigenvectors for loading into R
    tail -n+2 evec | sed -e 's/^  *//' -e 's/  *$//' -e 's/  */\t/g' -e 's/\t[^\t]*$//' > evec_R
    
    if test $(($num_evec + 1)) -gt "$(head -1 evec_R | tr '\t' '\n' | wc -l)"; then
    	dx-jobutil-report-error "ERROR: Desired number of eigenvectors exceeds the number of eigenvectors provided"
    fi
    
    # Assumes that "ancestry" and "evec" are formatted correctly.
    gen_ancestry.R ancestry_R evec_R $thresh $num_evec > $OUTDIR/$prefix.ancestry
    
    ancestry_out=$(dx upload $OUTDIR/$prefix.ancestry --brief)

    dx-jobutil-add-output ancestry_out "$ancestry_out" --class=file
}
