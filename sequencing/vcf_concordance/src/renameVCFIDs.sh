#!/bin/bash

# set to 1 if you want to have only 1 ID ("merging")
# set to 0 to allow for multiple IDs (checking for consistency)
UNIQUE=0

OLD_VCF_FILE="$1"
if [ -z "$3" ]; then
	NEW_VCF_FILE="$(echo $OLD_VCF_FILE | sed 's/\.vcf\.gz$/.renamed.vcf.gz/')"
else
	NEW_VCF_FILE="$3" 
fi

# File of id -> new id (tab separated)
ID_TRANS="$2"


# DO NOT EDIT ANYTHING BELOW THIS LINE
# ====================================================================

vcf_header=$(mktemp)
vcf_newheader=$(mktemp)

id_t=$(mktemp)
if [ -n "$UNIQUE" ] && [ "$UNIQUE" -ne 0 ]; then
	cat $ID_TRANS | sort -k 2,2 -u > $id_t
else
	cat $ID_TRANS > $id_t
fi

#cut -f1 $id_t

# We should filter the VCF to include only those IDs in the idlist
zcat "$OLD_VCF_FILE" | head -5000 | grep -E '^#' > $vcf_header
# In this case, we need to filter
if test $(join -1 2 -2 1 <(grep '#CHROM' $vcf_header | cut -f10- | tr '\t' '\n' | grep -n '.' | tr ':' '\t' | sort -k2) <(cat "$id_t" | sort) | wc -l) -ne $(grep '#CHROM' $vcf_header | cut -f10- | tr '\t' '\n' | wc -l); then
	echo "Filtering VCF..."

	vcf_filtered="$(mktemp).vcf.gz"
	JAVA_OPTIONS="-d64 -Xms512m -Xmx8G" GenomeAnalysisTK -T SelectVariants -R /data/public/GATK_resources/human_g1k_v37_decoy.fasta -V $OLD_VCF_FILE --sample_file <(cut -f1 $id_t) -o $vcf_filtered -nt 4
	tabix -p vcf $vcf_filtered
	OLD_VCF_FILE=$vcf_filtered
	zcat "$OLD_VCF_FILE" | head -5000 | grep -E '^#' > $vcf_header
fi

cat <(grep -Ev '#CHROM' $vcf_header) <(paste <(echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT") <(join -1 2 -2 1 <(grep '#CHROM' $vcf_header | cut -f10- | tr '\t' '\n' | grep -n '.' | tr ':' '\t' | sort -k2) <(cat "$id_t" | sort) | tr ' ' '\t' | sort -g -k2,2 | cut -f 3 | tr '\n' '\t' | sed 's/\t$//') ) > $vcf_newheader

if test "$(tail -1 $vcf_header | tr '\t' '\n' | wc -l)" -ne "$(tail -1 $vcf_newheader | tr '\t' '\n' | wc -l)"; then
	echo "Some columns not mapped.. exiting"
	exit
fi

tabix -r $vcf_newheader "$OLD_VCF_FILE" > "$NEW_VCF_FILE"
tabix -p vcf "$NEW_VCF_FILE"

rm $vcf_header
rm $vcf_newheader
rm $id_t

if test ! -z "$vcf_filtered"; then
	rm $vcf_filtered
	rm $vcf_filtered.tbi
fi
