#!/usr/bin/env python
"""Script to reformat the ClinVar TSV file into two VCF files, one for b37 and
one for b38.
"""
__author__ = "Thomas Nate Person"
__license__ = "GPLv3"
__email__ = "thomas.n.person@gmail.com"

import re
import sys
import time
import gzip

MonthYear =time.strftime("%b")+time.strftime("%Y")
DayMonthYear =time.strftime("%d")+time.strftime("%b")+time.strftime("%Y")

dbFile = gzip.open(sys.argv[1],'r')
db_b37_vcf = open(sys.argv[1].replace("txt.gz","")+'b37.vcf','w')
db_b38_vcf = open(sys.argv[1].replace("txt.gz","")+'b38.vcf','w')


for line in dbFile:
	if line.startswith("#AlleleID"):
		headerFields="|".join(line.lstrip("#").strip().split("\t"))
		headerFields="Allele|"+headerFields
	else:
		break

db_b37_vcf.write("##fileformat=VCFv4.2\n")
db_b37_vcf.write("##INFO=<ID=ClinVar.TSV."+DayMonthYear+",Number=.,Type=String,Description=\"ClinVar.TSV." + DayMonthYear+": '"+headerFields+"'\">\n")
db_b37_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

db_b38_vcf.write("##fileformat=VCFv4.2\n")
db_b38_vcf.write("##INFO=<ID=ClinVar.TSV."+DayMonthYear+",Number=.,Type=String,Description=\"ClinVar.TSV." + DayMonthYear+": '"+headerFields+"'\">\n")
db_b38_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

CHROM =-1
POS = -1
REF = -1
ALT = -1
ASSEMBLY = -1

dbFile.seek(0)

new_vcf=dict()
new_vcf_b37=dict()
new_vcf_b38=dict()
new_vcf["GRCh37"]=new_vcf_b37
new_vcf["GRCh38"]=new_vcf_b38
new_vcf["source"]="ClinVar.TSV."+DayMonthYear+"="

for line in dbFile:
	if line.startswith("#"):
		fields = [x.strip() for x in line.split('\t')]
		for i, f in enumerate(fields):
			if f == 'Chromosome':
				CHROM = i
			elif f == 'Start':
				POS = i
			elif f == 'ReferenceAllele':
				REF = i
			elif f == 'AlternateAllele':
				ALT = i
			elif f == 'Assembly':
				ASSEMBLY = i
	else:
		fields = [x.strip() for x in line.split('\t')]
		if fields[REF]==fields[ALT]:
			continue
		elif fields[ASSEMBLY]=='NCBI36':
			continue

		if bool(re.search('^[AGCT]+$', fields[REF])) and bool(re.search('^[AGCT]+$', fields[ALT])):
			site_id=fields[CHROM]+"\t"+fields[POS]+"\t.\t"+fields[REF]+"\t"
			variant_id = fields[CHROM]+"\t"+fields[POS]+"\t"+fields[REF]+"\t"+fields[ALT]

			INFO = fields[ALT]+"|"+"|".join(fields)
			INFO = INFO.replace(';','/')
			INFO = INFO.replace(',','/')
			INFO = INFO.replace(' ','_')

			if site_id in new_vcf[fields[ASSEMBLY]]:
				site = new_vcf[fields[ASSEMBLY]][site_id]
				site[site_id].add(fields[ALT])

				if variant_id in site:
					site[variant_id]=site[variant_id]+","+INFO
				else:
					site[variant_id]=INFO

				new_vcf[fields[ASSEMBLY]][site_id]=site

			else:
				site = dict()
				site[variant_id]=INFO
				site[site_id]=set()
				site[site_id].add(fields[ALT])
				new_vcf[fields[ASSEMBLY]][site_id]=site

for site_id, site in new_vcf['GRCh37'].viewitems():
	outsite= site_id + ",".join(site[site_id])+"\t.\t.\t"+new_vcf["source"]
	site.pop(site_id,None)
	outsite=outsite+",".join(site.values())+"\n"
	db_b37_vcf.write(outsite)

db_b37_vcf.close()

for site_id, site in new_vcf['GRCh38'].viewitems():
	outsite= site_id + ",".join(site[site_id])+"\t.\t.\t"+new_vcf["source"]
	site.pop(site_id,None)
	outsite=outsite+",".join(site.values())+"\n"
	db_b38_vcf.write(outsite)

db_b38_vcf.close()
exit()
