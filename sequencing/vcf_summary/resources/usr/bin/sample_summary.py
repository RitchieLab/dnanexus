#!/usr/bin/env python

import sys
import re

def parse_genestr(cgstr):
	cpg = cgstr.strip().split()
	# should only happen at EOF
	if not cpg:
		return None
	return (cpg[0], int(cpg[1]), [] if len(cpg) == 2 else cpg[2].split(','))

if __name__ == "__main__":
	
	gene_f = file(sys.argv[1], 'r')
	vcf_f = file(sys.argv[2], 'r')
		
	field_list = []
	
	# read the VCF header
	for l in vcf_f:
		if l.startswith("##"):
			continue
		elif l.startswith("#"):
			# read the sample line
			field_list = l.strip().split('\t')	
			break
	
	# samp_data is a dictionary of sample to a list containing the following:
	# 0: # of called variants (of any genotype)
	# 1: # of non-homref variants
	# 2: # of alt alleles
	# 3: a set of genes containing a non-homref call		
	samp_data_raw = {field_list[i] : [0, 0, 0, set([])] for i in range(9, len(field_list))}
	samp_data_filt = {field_list[i] : [0, 0, 0, set([])] for i in range(9, len(field_list))}
	
	curr_genepos=parse_genestr(gene_f.readline())	
		
	# read the rest of the VCF file
	for l in vcf_f:
		fields = l.strip().split('\t')
		
		format_dict = {v: i for (i, v) in enumerate(fields[8].split(':'))}
		gt_pos = format_dict["GT"]
		ft_pos = format_dict.get("FT", None)
		
		filter_pass = (fields[6] == "PASS" or fields[6] == ".")
		
		curr_genes = []
		move_genes = False
		if curr_genepos is not None and curr_genepos[0] == fields[0] and curr_genepos[1] == int(fields[1]):
			curr_genes = curr_genepos[2]
			move_genes = True
		
		geno_sep = re.compile('/|\|')
		
		for (field_pos, geno) in enumerate(fields[9:]):
			
			cs = field_list[field_pos+9]
			genof = geno.split(':')
			gt = genof[gt_pos]
			
			ft_pass = (ft_pos is None or genof[ft_pos] == "PASS")
			
			gt_calls = geno_sep.split(gt)
			gt_called = all( ( (g != ".") for g in gt_calls) )
			gt_alt = sum((int(g) != 0 for g in gt_calls)) if gt_called else 0
			
			samp_data_raw[cs][0] += gt_called
			samp_data_raw[cs][1] += (gt_alt > 0)
			samp_data_raw[cs][2] += gt_alt
			
			if filter_pass and ft_pass:
				samp_data_filt[cs][0] += gt_called
				samp_data_filt[cs][1] += (gt_alt > 0)
				samp_data_filt[cs][2] += gt_alt
			
			if gt_alt > 0:
				for gn in curr_genes:
					samp_data_raw[cs][3].add(gn)
					if filter_pass and ft_pass:
						samp_data_filt[cs][3].add(gn)
			
		
		if move_genes:
			curr_genepos = parse_genestr(gene_f.readline())
			
	# Now, print the data that we have
	print "Sample\tFilt_Call\tFilt_Vars\tFilt_Alts\tFilt_Num_gene\tFilt_gene\tRaw_Call\tRaw_Vars\tRaw_Alts\tRaw_Num_gene\tRaw_gene"
	
	for k in samp_data_raw.keys():
		rd = samp_data_raw[k]
		fd = samp_data_filt[k]
		print '\t'.join( (k, str(fd[0]), str(fd[1]), str(fd[2]), str(len(fd[3])), ','.join(fd[3]), str(rd[0]), str(rd[1]), str(rd[2]), str(len(rd[3])), ','.join(rd[3])) )
	
