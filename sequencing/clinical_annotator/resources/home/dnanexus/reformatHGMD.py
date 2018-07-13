#!/usr/bin/env python
"""Script to reformat the HGMD PRO vcf file for one annotation per line.
"""

__author__ = "Thomas Nate Person"
__license__ = "GPLv3"
__email__ = "thomas.n.person@gmail.com"

import sys
import gzip
import collections

vcf_file = None
if sys.argv[1].endswith(".gz"):

	vcf_file = gzip.open(sys.argv[1],'r')

else:
	
	vcf_file = open(sys.argv[1],'r')

source=""
ids = {}
last_site_id = ""
hgmd_annotations=collections.OrderedDict()
new_vcf=collections.OrderedDict()
for line in vcf_file:

	if line.startswith("#"):

		if line.startswith("##source"):

			source = line.rstrip("#").strip().split("=")[-1]
		
		if line.startswith("#C"):

			print "##INFO=<ID={},Number=.,Type=String,Description=\"HGMD Annotations. Format: Allele|CLASS|MUT|GENE|STRAND|DNA|PROT|PHEN|RANKSCORE|ID\">".format(source)

		print line.strip()

	else:
		
		fields=line.strip().split("\t")
		variant_id=fields[0]+":"+fields[1]+":"+fields[3]+":"+fields[4]
		site_id=fields[0]+"\t"+fields[1]+"\t.\t"+fields[3]+"\t"
		hgmd_a=[fields[4],".",".",".",".",".",".",".",".",".",fields[2]]
		fields[2]="."
		info_fields=fields[7].split(";")
		for i_f in info_fields:

			if i_f.startswith("CLASS"):
				hgmd_a[1]=i_f.strip().split("=")[-1]
			elif i_f.startswith("MUT"):
				hgmd_a[2]=i_f.strip().split("=")[-1]
			elif i_f.startswith("GENE"):
				hgmd_a[3]=i_f.strip().split("=")[-1]
			elif i_f.startswith("STRAND"):
				hgmd_a[4]=i_f.strip().split("=")[-1]
			elif i_f.startswith("DNA"):
				hgmd_a[5]=i_f.strip().split("=")[-1]
			elif i_f.startswith("PROT"):
				hgmd_a[6]=i_f.strip().split("=")[-1]
			elif i_f.startswith("DB"):
				hgmd_a[7]=i_f.strip().split("=")[-1]
			elif i_f.startswith("PHEN"):
				hgmd_a[8]=i_f.strip().split("=")[-1]
			elif i_f.startswith("RANKSCORE"):
				hgmd_a[9]=i_f.strip().split("=")[-1]

		if site_id == last_site_id:

			site = new_vcf[site_id]
			if variant_id in site:
				
				site[variant_id]=site[variant_id]+","+"|".join(hgmd_a)
			
			else:

				site[variant_id]="|".join(hgmd_a)

			site[site_id].add(fields[4])
			new_vcf[site_id]=site


		else:

			site = {}
			site[variant_id]="|".join(hgmd_a)
			site[site_id]=set()
			site[site_id].add(fields[4])
			new_vcf[site_id]=site
			last_site_id = site_id


for site_id, site in new_vcf.viewitems():

	outsite= site_id + ",".join(site[site_id])+"\t.\t.\t"+source+"="
	site.pop(site_id,None)
	outsite=outsite+",".join(site.values())
	print outsite

sys.stdout.flush()
exit()
