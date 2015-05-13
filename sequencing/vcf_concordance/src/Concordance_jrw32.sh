#!/bin/bash

#Debugging only, remove later!
#set -x
set -e

#------------------- EDIT THESE VARIABLES (if desired) -----------------------

# # of processors to use on various GATK steps
N_PROC=1
# Location of the GATK reference needed for basically everything
REF="/data/public/GATK_resources/human_g1k_v37_decoy.fasta"

#------------------- Function Definitions for concordance metrics ------------

# Make script stop on errors and undefined variables
set -u

#What are all these Concordance metrics
#NRS: non-ref sensitivity(based on COMP set),  NRS_Reverse: non-ref sensitivity (based on EVAL set),  NRD: non-ref discrepancy,  OGC: Overall Genotype Conc, OGCM: Overall Genotype Conc with missing as discordant,  HETHET: Het calls over both sets w/o missing, HETMISSING: Het calls over both sets w/ missing  ,  NRC: Non-ref concordance (precision),  RNRC: Reverse non-ref Concordance(precision), Recall: ,RRecall: Reverse Recall 

function calc_NRS(){
	echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($16+$17+$22+$23)/($10+$11+$16+$17+$22+$23+$4+$5+$28+$29)}' 2>/dev/null
}	

function calc_NRS_Reverse(){
	echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($16+$17+$22+$23)/($15+$21+$16+$17+$22+$23+$14+$20+$18+$24)}' 2>/dev/null
}

function calc_NRD(){
	 echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print (1-($23+$16)/($10+$11+$15+$16+$17+$21+$22+$23))}' 2>/dev/null
}

function calc_OGC(){
	echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($9+$16+$23)/($9+$10+$11+$15+$16+$17+$21+$22+$23)}' 2>/dev/null
}	

function calc_OGCM(){
	echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($9+$16+$23+$2)/($9+$10+$11+$8+$15+$16+$17+$14+$21+$22+$23+$20+$3+$4+$5+$2)}' 2>/dev/null
}

function calc_OGCMU(){
	echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($9+$16+$23+$2+$30+$26+$6)/($9+$10+$11+$8+$12+$15+$16+$17+$14+$18+$21+$22+$23+$20+$24+$3+$4+$5+$2+$6+$27+$28+$29+$26+$30)}' 2>/dev/null
}

function calc_HETHET(){
        echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($16)/($16+$17+$10+$22+$15)}' 2>/dev/null
}

function calc_HETMissing(){
        echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($16)/($15+$16+$17+$14+$10+$22+$4)}' 2>/dev/null
}

function calc_NRC(){
        echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($16+$23)/($15+$16+$17+$21+$22+$23)} ' 2>/dev/null
}

function calc_RNRC(){
        echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($16+$23)/($10+$11+$16+$17+$22+$23)} ' 2>/dev/null
}

function calc_Recall(){
        echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($16+$23)/($10+$11+$16+$17+$22+$23+$4+$5+$28+$29)} ' 2>/dev/null
}

function calc_RevRecall(){
        echo $1 | awk 'BEGIN{OFMT = "%.5f"} { print ($16+$23)/($15+$21+$16+$17+$22+$23+$14+$20+$18+$24)} ' 2>/dev/null
}
function print_stats(){
	paste <(calc_NRS "$1") <(calc_NRS_Reverse "$1") <(calc_NRD "$1") <(calc_OGC "$1") <(calc_OGCM "$1") <(calc_OGCMU "$1") <(calc_HETHET "$1") <(calc_HETMissing "$1") <(calc_NRC "$1") <(calc_RNRC "$1") <(calc_Recall "$1") <(calc_RevRecall "$1")
}

####### v---- COMMANDS BELOW ----v
# Abandon all hope, ye who enter here

# inputs will be $1 and $2
# $1 will certainly be a VCF file, $2 may be a VCF file or and ID mapping of dupes for the given VCF

# First, let's see if we have a VCF for $2 or an ID mapping
CAT_CMD="cat"
if test "$(echo $2 | grep '.gz$')"; then
	CAT_CMD="zcat"
fi

CAT1_CMD="cat"
if test "$(echo $1 | grep '.gz$')"; then
	CAT1_CMD="zcat"
fi

VCF1=$PWD/$1
if test "$(echo $1 | grep '^/')"; then
	VCF1=$1
fi

VCF2=$PWD/$2
if test "$(echo $2 | grep '^/')"; then
	VCF2=$2
fi

# creatE an empty working directory and go there
CALL_DIR="$PWD"
WKDIR=$(mktemp -d)
cd $WKDIR

ID_LIST_ALL=$(mktemp)

PARALLEL_ARGS=""
if test "$N_PROC" -gt 1; then
	PARALLEL_ARGS="$PARALLEL_ARGS -nt $N_PROC"
fi

if test "$($CAT_CMD $VCF2 | head -1 | grep '^##fileformat=VCF')"; then
	ID_LIST1=$(mktemp)
	ID_LIST2=$(mktemp)
	
	$CAT1_CMD $VCF1 | head -1000 | grep '^#CHR' | cut -f10- | tr '\t' '\n' | sort > $ID_LIST1
	$CAT_CMD $VCF2 | head -1000 | grep '^#CHR' | cut -f10- | tr '\t' '\n' | sort > $ID_LIST2
	
	join $ID_LIST1 $ID_LIST2 > $ID_LIST_ALL
	
	if test $(cat $ID_LIST_ALL | wc -l) -eq 0; then
		echo "ERROR: No overlapping samples, nothing to compare!!"
		exit 1
	fi
	
	
	# at this point, we need to subset our 2 VCFs to only overlapping IDs
	
	#TODO: only perform this step if we really, truly NEED to!
	# Now, let's SelectVariants up in here!
	GenomeAnalysisTK-3.3-0 -T SelectVariants -R $REF $PARALLEL_ARGS \
    	-V $VCF1 -sf $ID_LIST_ALL -o VCF1.renamed.vcf.gz -l ERROR

	GenomeAnalysisTK-3.3-0 -T SelectVariants -R $REF $PARALLEL_ARGS \
    	-V $VCF2 -sf $ID_LIST_ALL -o VCF2.renamed.vcf.gz -l ERROR

	rm $ID_LIST1
	rm $ID_LIST2
	

	
else
	# First, let's get a listing of all IDs that are duplicated
	$CAT_CMD $VCF2 | cut -f2 | sort | uniq -d > $ID_LIST_ALL
	if test $(cat $ID_LIST_ALL | wc -l) -eq 0; then
		echo "ERROR: No duplicated samples, nothing to compare!!"
		exit 1
	fi
	
	ID_LIST1=$(mktemp)
	ID_LIST2=$(mktemp)
	
	# This command selects one duplicate to go to stderr and one to stdout (the shuf combined w/ stable sort randomizes it!)
	$CAT_CMD $VCF2 | shuf | sort -k2,2 -s | uniq -f1 -d -D | awk 'BEGIN {pv=""; nl=0;} {if($2!=pv){print $0; pv=$2; nl=1;}else if(nl==1){print $0 > "/dev/stderr"; nl=0;}}' >$ID_LIST1 2>$ID_LIST2
	
	ORIG_IDLIST1=$(mktemp)
	cut -f1 $ID_LIST1 > $ORIG_IDLIST1
	
	ORIG_IDLIST2=$(mktemp)
	cut -f1 $ID_LIST2 > $ORIG_IDLIST2
	
	GenomeAnalysisTK-3.3-0 -T SelectVariants -R $REF $PARALLEL_ARGS \
    	-V $VCF1 -sf $ORIG_IDLIST1 -o VCF1_orig.vcf.gz -l ERROR

	GenomeAnalysisTK-3.3-0 -T SelectVariants -R $REF $PARALLEL_ARGS \
    	-V $VCF1 -sf $ORIG_IDLIST2 -o VCF2_orig.vcf.gz -l ERROR

	$CALL_DIR/renameVCFIDs.sh VCF1_orig.vcf.gz $ID_LIST1 VCF1.renamed.vcf.gz
	$CALL_DIR/renameVCFIDs.sh VCF2_orig.vcf.gz $ID_LIST2 VCF2.renamed.vcf.gz
	
	rm $ID_LIST1
	rm $ID_LIST2
	rm $ORIG_IDLIST1
	rm $ORIG_IDLIST2

fi

#RUN GATK GENOTYPE CONCORDANCE ALL
#GATK Genotype Concordance (raw)
#this uses all records, so the filters will be ignored
GenomeAnalysisTK-3.3-0-jrw32  -T GenotypeConcordance -R $REF --comp VCF1.renamed.vcf.gz --eval VCF2.renamed.vcf.gz  --ignoreFilters --out GATK.raw.txt -l ERROR

#GATK Genotype Concordance (all_filter)
#gfe This will apply Genotype filters to EVAL set (anything true in the expression will be a no call
#gfc This will apply Genotype filters to the Comp set (anything true in the expression will be a no call)
#Add FT expression in gfc and gfe
GenomeAnalysisTK-3.3-0-jrw32  -T GenotypeConcordance -R $REF --comp VCF1.renamed.vcf.gz --eval VCF2.renamed.vcf.gz -gfc 'FT!="PASS"' -gfe  'FT!="PASS"' --out GATK.filtered.txt -l ERROR

# get SNP-only VCFs
GenomeAnalysisTK-3.3-0 -T SelectVariants -selectType SNP -V VCF1.renamed.vcf.gz -R $REF -o VCF1.renamed.SNP.vcf.gz -l ERROR
GenomeAnalysisTK-3.3-0 -T SelectVariants -selectType SNP -V VCF2.renamed.vcf.gz -R $REF -o VCF2.renamed.SNP.vcf.gz -l ERROR


#RUN GATK GENOTYPE CONCORDANCE SNP ONLY
#GATK Genotype Concordance (raw)
#this uses all records, so the filters will be ignored
GenomeAnalysisTK-3.3-0-jrw32  -T GenotypeConcordance -R $REF --comp VCF1.renamed.SNP.vcf.gz --eval VCF2.renamed.SNP.vcf.gz  --ignoreFilters --out GATK.raw.SNP.txt -l ERROR

#GATK Genotype Concordance (all_filter)
#gfe This will apply Genotype filters to EVAL set (anything true in the expression will be a no call
#gfc This will apply Genotype filters to the Comp set (anything true in the expression will be a no call)
#Add FT expression in gfc and gfe
GenomeAnalysisTK-3.3-0-jrw32  -T GenotypeConcordance -R $REF --comp VCF1.renamed.SNP.vcf.gz --eval VCF2.renamed.SNP.vcf.gz -gfc 'FT!="PASS"' -gfe 'FT!="PASS"'  --out GATK.filtered.SNP.txt -l ERROR

# get INDEL-only VCF
GenomeAnalysisTK-3.3-0 -T SelectVariants -selectType INDEL -selectType MIXED -selectType MNP -V VCF1.renamed.vcf.gz -R $REF -o VCF1.renamed.INDEL.vcf.gz -l ERROR
GenomeAnalysisTK-3.3-0 -T SelectVariants -selectType INDEL -selectType MIXED -selectType MNP -V VCF2.renamed.vcf.gz -R $REF -o VCF2.renamed.INDEL.vcf.gz -l ERROR

#RUN GATK GENOTYPE CONCORDANCE INDEL ONLY
#GATK Genotype Concordance (raw)
#this uses all records, so the filters will be ignored
GenomeAnalysisTK-3.3-0-jrw32  -T GenotypeConcordance -R $REF --comp VCF1.renamed.INDEL.vcf.gz --eval VCF2.renamed.INDEL.vcf.gz  --ignoreFilters --out GATK.raw.INDEL.txt -l ERROR

#GATK Genotype Concordance (all_filter)
#gfe This will apply Genotype filters to EVAL set (anything true in the expression will be a no call
#gfc This will apply Genotype filters to the Comp set (anything true in the expression will be a no call)
#Add FT expression in gfc and gfe
GenomeAnalysisTK-3.3-0-jrw32  -T GenotypeConcordance -R $REF --comp VCF1.renamed.INDEL.vcf.gz --eval VCF2.renamed.INDEL.vcf.gz -gfc 'FT!="PASS"' -gfe 'FT!="PASS"' --out GATK.filtered.INDEL.txt -l ERROR

#RUN ALL METRICS ON RAW AND FILTERS FOR  ALL AND SNP ONLY 
# Extract only the info and lines from #:GATKTable:GenotypeConcordance_Counts 

#Create Conconordance metrics output files with columns:  sample NRD NRS NRS_Reverse OGC OGCM
N_SAMPLES=$(cat $ID_LIST_ALL | wc -l)

paste <(printf "%-16s" $(echo "Metric:")) <(echo "NRS") <(echo "NRS-R") <(echo "NRD") <(echo "OGC") <( echo "OGC-M") <(echo "OGC-MU") <( echo "HHC") <(echo "HHC-M") <(echo "NRC") <(echo "NRC-R") <(echo "REC") <(echo "REC-R")

for f in filtered raw filtered.SNP raw.SNP filtered.INDEL raw.INDEL; do

	FN="GATK.$f.txt"
	START_LINE=$(grep -n "^ALL" $FN | head -n 2 | tail -1 | cut -d':' -f1)
	SAMPL_LINE=$(grep -n "^Sample" $FN | head -n 2 | tail -1 | cut -d':' -f1)
	
	paste <(printf "%-16s" $(echo "$f:") ) <(print_stats "$(head -n$START_LINE $FN | tail -1)" )
	
	
	touch $FN.by_sample	
	RAW_SAMP=$(mktemp)
	head -n$((SAMPL_LINE+N_SAMPLES+1)) $FN | tail -n $((N_SAMPLES+1)) | grep -v '^ALL' > $RAW_SAMP
	
	while read l; do
		paste <(printf "%-16s" $(echo "$l" | sed -e 's/^[ \t]*//' -e 's/[ \t].*/:/')) <(print_stats "$l") >> $FN.by_sample
	done < $RAW_SAMP
	rm $RAW_SAMP
done
	
echo
echo "Filtered"
echo "========"
cat GATK.filtered.txt.by_sample
echo
echo "Raw"
echo "========"
cat GATK.raw.txt.by_sample
echo
echo "Filtered (SNP only)"
echo "==================="
cat GATK.filtered.SNP.txt.by_sample
echo
echo "Raw (SNP only)"
echo "=============="
cat GATK.raw.SNP.txt.by_sample
echo
echo "Filtered (INDEL only)"
echo "==================="
cat GATK.filtered.INDEL.txt.by_sample
echo
echo "Raw (INDEL only)"
echo "==================="
cat GATK.raw.INDEL.txt.by_sample

# and please clean up after yourself like a good little boy or girl!
rm $ID_LIST_ALL
cd - 2>&1 >/dev/null
rm -rf $WKDIR




