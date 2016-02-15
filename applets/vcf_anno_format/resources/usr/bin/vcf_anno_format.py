#!/usr/bin/env python

import collections
import sys


chmPosAlleles = chmZoneRegions = genes = None
if len(sys.argv) <= 4:
	print "usage: %s <prefix> <positions> <regions> <genes> [option] [option] ..." % (sys.argv[0],)
	print "  position/region/gene filter files may be omitted with '.'"
	print "  options are:"
	print "    clinvar , dbnsfp , sift , snpeff"
	print "      include this annotation type in the output"
	print "    extra"
	print "      include a final column to contain any additional annotations"
	print "      of the same type for the same allele"
	sys.exit(1)

empty = list()
prefix = sys.argv[1].strip()
options = set(sys.argv[a].strip().lower() for a in xrange(5,len(sys.argv)))

if sys.argv[2] != '.':
	print "reading position filter from '%s' ..." % (sys.argv[2],)
	chmPosAlleles = collections.defaultdict(lambda: collections.defaultdict(set))
	l = p = a = 0
	with open(sys.argv[2],'rU') as pfile:
		for line in pfile:
			l += 1
			cols = line.split()
			try:
				chm = cols[0].strip().upper()
				if chm.startswith('CHR'):
					chm = chm[3:]
				pos = long(cols[1])
				alleles = empty
				if len(cols) > 2:
					alleles = set(a.strip().upper() for a in cols[2].split(','))
				p -= len(chmPosAlleles[chm])
				a -= len(chmPosAlleles[chm][pos])
				chmPosAlleles[chm][pos].update(alleles)
				p += len(chmPosAlleles[chm])
				a += len(chmPosAlleles[chm][pos])
			except ValueError:
				if l != 1:
					raise
			#try/except
		#for line
	#with pfile
	print "... OK: %d positions, %d specified alleles" % (p,a)
#if arg1

if sys.argv[3] != '.':
	print "reading region filter from '%s' ..." % (sys.argv[3],)
	chmZoneRegions = collections.defaultdict(lambda: collections.defaultdict(list))
	l = r = 0
	with open(sys.argv[3],'rU') as rfile:
		for line in rfile:
			l += 1
			cols = line.split()
			try:
				chm = cols[0].strip().upper()
				if chm.startswith('CHR'):
					chm = chm[3:]
				start = long(cols[1])
				stop = long(cols[2])
				r += 1
				for zone in xrange(start/100000, stop/100000+1):
					chmZoneRegions[chm][zone].append( (start,stop) )
			except ValueError:
				if l != 1:
					raise
			#try/except
		#for line
	#with rfile
	print "... OK: %d regions" % (r,)
#if arg2

if sys.argv[4] != '.':
	print "reading gene filter from '%s' ..." % (sys.argv[4],)
	genes = set()
	l = 0
	with open(sys.argv[4],'rU') as gfile:
		for line in gfile:
			l += 1
			cols = line.split()
			gene = cols[0].strip().upper()
			genes.add(gene)
		#for line
	#with gfile
	print "... OK: %d genes" % (len(genes),)
#if arg3


print "opening output files ..."
outputANN = outputSIFT = outputNSFP = outputCLINVAR = None
headers = ('CHROM','POS','REF','ALT','AF')

if 'snpeff' in options:
	headers2 = ('GeneName','GeneID','NumTransInGene','PctTransAffected')
	headers3 = ('Annotation','Impact','GeneName','GeneID','FeatureType','FeatureID','TranscriptBioType','Rank','HGVS.c','HGVS.p','cDNA_pos','CDS_pos','AA_pos','Distance','ErrWarnInfo')
	outputANN = open(prefix+'.SNPeff.txt','w')
	outputANN.write("#%s\t%s\t%s\t%s\t%s\t%s%s\n" % (
		'\t'.join(headers),
		'\t'.join( ('LOF'+s) for s in headers2 ),
		'\t'.join( ('NMD'+s) for s in headers2 ),
		'\t'.join( (s+'1') for s in headers3 ),
		'\t'.join( (s+'2') for s in headers3 ),
		'\t'.join( (s+'3') for s in headers3 ),
		('\tEXTRA' if ('extra' in options) else ''),
	))
#if snpeff

if 'sift' in options:
	headers2 = ('Annotation','RawSIFTScore','GeneName','EnsemblGeneID','SNPType','EnsemblTranscriptID','SNPPlacement')
	outputSIFT = open(prefix+'.SIFT.txt','w')
	outputSIFT.write("#%s\t%s%s\n" % (
		'\t'.join(headers),
		'\t'.join(headers2),
		('\tEXTRA' if ('extra' in options) else ''),
	))
#if sift

if 'dbnsfp' in options:
	headers2 = ('SIFT','Polyphen2_HDIV','Polyphen2_HVAR','LRT','MutationTaster','MutationAssessor','FATHMM','MetaSVM','MetaLR','PROVEAN','genename','Ensembl_geneid','Ensembl_transcriptid','clinvar_clnsig','clinvar_trait')
	outputNSFP = open(prefix+'.dbNSFP.txt','w')
	outputNSFP.write("#%s\t%s%s\n" % (
		'\t'.join(headers),
		'\t'.join(headers2),
		('\tEXTRA' if ('extra' in options) else ''),
	))
#if dbnsfp

if 'clinvar' in options:
	headers2 = ('Significance','ReviewStatus','GeneID')
	outputCLINVAR = open(prefix+'.CLINVAR.txt','w')
	outputCLINVAR.write("#%s\t%s%s\n" % (
		'\t'.join(headers),
		'\t'.join(headers2),
		('\tEXTRA' if ('extra' in options) else ''),
	))
#if clinvar

print "... OK"


print "reading input VCF ..."
for line in sys.stdin:
	if line.startswith('#'):
		continue
	
	cols = line.rstrip('\r\n').split('\t')
	chm = cols[0].strip().upper()
	pos = long(cols[1])
	zone = pos / 100000
	
	if chmZoneRegions != None:
		skip = True
		if (chm in chmZoneRegions) and (zone in chmZoneRegions[chm]):
			for region in chmZoneRegions[chm][zone]:
				if (region[0] <= pos) and (region[1] >= pos):
					skip = False
					break
				#if region match
			#for each region
		#if chm/zone match
		if skip:
			continue
	#if region filter
	
	ref = cols[3].strip().upper()
	alts = list(allele.strip().upper() for allele in cols[4].split(','))
	freqs = empty
	annoAltData = collections.defaultdict(lambda: collections.defaultdict(list))
	for infotag in cols[7].split(';'):
		infotag = infotag.split('=',1)
		if infotag[0] == 'AF':
			freqs = infotag[1].split(',')
		elif infotag[0] in ('LOF','NMD'):
			data = infotag[1].strip('()').split('|')
			if genes != None:
				if data[0].strip().upper() not in genes:
					continue
			#if genes filter
			annoAltData[infotag[0]][''].append(data)
		elif infotag[0] in ('ANN','SIFT','dbNSFP','CLINVAR'):
			for infosub in infotag[1].split(','):
				data = infosub.split('|')
				if chmPosAlleles != None:
					if (chm not in chmPosAlleles) or (pos not in chmPosAlleles[chm]) or ((len(chmPosAlleles[chm][pos]) > 0) and (data[0].strip().upper() not in chmPosAlleles[chm][pos])):
						continue
				#if pos/allele filter
				if genes != None:
					# it just so happens that the gene is in the 4th subfield of all annotation types except dbNSFP, where it's 12th
					if data[11 if (infotag[0] == 'dbNSFP') else 3].strip().upper() not in genes:
						continue
				#if genes filter
				annoAltData[infotag[0]][data[0]].append(data)
			#for each infosub
		#if infotag matters
	#for each infotag
	altFreq = { alts[a]:(freqs[a] if (a < len(freqs)) else None) for a in xrange(len(alts)) }
	
	
	for alt in alts:
		if outputANN and (len(annoAltData['ANN'][alt]) > 0):
			# 14: Protein_position / Protein_len: Position and number of AA (one based, including START, but not STOP)
			annoAltData['ANN'][alt].sort(key=(lambda data: int((data[13]+'/').split('/',2)[1] or 0)), reverse=True)
			lof = annoAltData['LOF'][''][0] if (len(annoAltData['LOF']['']) >= 1) else empty
			nmd = annoAltData['NMD'][''][0] if (len(annoAltData['NMD']['']) >= 1) else empty
			data1 = annoAltData['ANN'][alt][0] if (len(annoAltData['ANN'][alt]) >= 1) else empty
			data2 = annoAltData['ANN'][alt][1] if (len(annoAltData['ANN'][alt]) >= 2) else empty
			data3 = annoAltData['ANN'][alt][2] if (len(annoAltData['ANN'][alt]) >= 3) else empty
			outputANN.write("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n" % (
				chm, pos, ref, alt, altFreq[alt],
				'\t'.join( ((lof[i] or 'NA') if (i < len(lof)) else 'NA') for i in xrange(0,4) ),
				'\t'.join( ((nmd[i] or 'NA') if (i < len(nmd)) else 'NA') for i in xrange(0,4) ),
				'\t'.join( ((data1[i] or 'NA') if (i < len(data1)) else 'NA') for i in xrange(1,16) ),
				'\t'.join( ((data2[i] or 'NA') if (i < len(data2)) else 'NA') for i in xrange(1,16) ),
				'\t'.join( ((data3[i] or 'NA') if (i < len(data3)) else 'NA') for i in xrange(1,16) ),
				('\t'+','.join( '|'.join(data) for data in annoAltData['ANN'][alt][3:] )) if ('extra' in options) else '',
			))
		#if snpeff
		if outputSIFT and (len(annoAltData['SIFT'][alt]) > 0):
			annoAltData['SIFT'][alt].sort(key=(lambda data: sum((len(s) > 0) for s in data)), reverse=True)
			data1 = annoAltData['SIFT'][alt][0]
			outputSIFT.write("%s\t%d\t%s\t%s\t%s\t%s%s\n" % (
				chm, pos, ref, alt, altFreq[alt],
				'\t'.join( ((data1[i] or 'NA') if (i < len(data1)) else 'NA') for i in xrange(1,8) ),
				('\t'+','.join( '|'.join(data) for data in annoAltData['SIFT'][alt][1:] )) if ('extra' in options) else '',
			))
		#if sift
		if outputNSFP and (len(annoAltData['dbNSFP'][alt]) > 0):
			annoAltData['dbNSFP'][alt].sort(key=(lambda data: sum((len(s) > 0) for s in data)), reverse=True)
			data1 = annoAltData['dbNSFP'][alt][0]
			outputNSFP.write("%s\t%d\t%s\t%s\t%s\t%s%s\n" % (
				chm, pos, ref, alt, altFreq[alt],
				'\t'.join( ((data1[i] or 'NA') if (i < len(data1)) else 'NA') for i in xrange(1,16) ),
				('\t'+','.join( '|'.join(data) for data in annoAltData['dbNSFP'][alt][1:] )) if ('extra' in options) else '',
			))
		#if dbnsfp
		if outputCLINVAR and (len(annoAltData['CLINVAR'][alt]) > 0):
			annoAltData['CLINVAR'][alt].sort(key=(lambda data: sum((len(s) > 0) for s in data)), reverse=True)
			data1 = annoAltData['CLINVAR'][alt][0]
			outputCLINVAR.write("%s\t%d\t%s\t%s\t%s\t%s%s\n" % (
				chm, pos, ref, alt, altFreq[alt],
				'\t'.join( ((data1[i] or 'NA') if (i < len(data1)) else 'NA') for i in xrange(1,4) ),
				('\t'+','.join( '|'.join(data) for data in annoAltData['CLINVAR'][alt][1:] )) if ('extra' in options) else '',
			))
		#if clinvar
	#for each alt
#for line
print "... OK"


if outputANN:
	outputANN.close()
if outputSIFT:
	outputSIFT.close()
if outputNSFP:
	outputNSFP.close()
if outputCLINVAR:
	outputCLINVAR.close()
