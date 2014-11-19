#!/usr/bin/env python

import sys
import string

def complement(ai):
	return ai.translate(string.maketrans("AaTtCcGg", "TTAAGGCC"))

if __name__ == "__main__":
	
	f = file(sys.argv[1],'r')
	
	for l in f:
		data = l.strip().split()
		if data[0] == 'strand':
			# get a list of the current alleles
			curr_alleles = set(data[3:5])
			ref_alleles = set(data[6:8])
			
			# we are monomorphic AND we have at least one of the alleles in common
			if '0' in curr_alleles and len(curr_alleles.intersection(ref_alleles)) > 0:
				common_allele = curr_alleles.intersection(ref_alleles).pop()
				ref_alleles.remove(common_allele)
				new_allele = ref_alleles.pop()
				
				if len(new_allele) == 1:
					if data[3] == '0':
						print data[2], data[3], data[4], new_allele, common_allele
					else:
						print data[2], data[3], data[4], common_allele, new_allele
			elif set([complement(a) for a in curr_alleles]) == ref_alleles:
				print data[2], data[3], data[4], complement(data[3]), complement(data[4])
				

