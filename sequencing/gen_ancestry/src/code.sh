#!/bin/bash

set -ex -o pipefail

main() {

    echo "Value of evec_fn: '$evec_fn'"
    echo "Value of ancestry_fn: '$ancestry_fn'"
    echo "Value of thresh: '$thresh'"

	if test -z "$prefix"; then
		prefix=$(dx describe --name "$evec_fn" | sed 's/\.evec$//')
	fi

    dx download "$evec_fn"
    dx download "$ancestry_fn"
    
    # remove "Other" from the ancestry; it's special!
    grep -v 'Other$' "$ancestry_fn_name" > "$ancestry_fn_name"_R.txt
    # Format eigenvectors for loading into R
    tail -n+2 "$evec_fn_name" | sed -e 's/^  *//' -e 's/  *$//' -e 's/  */\t/g' -e 's/\t[^\t]*$//' > "$evec_fn_name"_R.txt
    
    if test $(($num_evec + 1)) -gt "$(head -1 "$evec_fn_name"_R.txt | tr '\t' '\n' | wc -l)"; then
    	dx-jobutil-report-error "ERROR: Desired number of eigenvectors exceeds the number of eigenvectors provided"
    fi
    
    # Now, we will assume that "ancestry" and "evec" are formatted correctly.
    # NOTE: the following should go into an R script...
    
    gen_ancestry.R "$ancestry_fn_name"_R.txt "$evec_fn_name"_R.txt $thresh $num_evec > $OUTDIR/$prefix.ancestry
    
    ancestry_out=$(dx upload $OUTDIR/$prefix.ancestry --brief)

    dx-jobutil-add-output ancestry_out "$ancestry_out" --class=file
}
