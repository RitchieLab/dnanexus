#!/usr/bin/env python

import argparse
import collections
import gzip
import itertools
import string
import struct
import sys
import tempfile
import zlib


class zopen(object):
	
	def __init__(self, fileName, splitChar="\n", chunkSize=1024*1024):
		self._filePtr = open(fileName,'rb')
		self._splitChar = splitChar
		self._chunkSize = chunkSize
		self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
		self._text = ""
		self._lines = list()
	#__init__()
	
	
	def __del__(self):
		if self._filePtr:
			self._filePtr.close()
			self._filePtr = None
	#__del__()
	
	
	def __enter__(self):
		return self
	#__enter__()
	
	
	def __exit__(self, excType, excVal, excTrace):
		pass
	#__exit__()
	
	
	def __iter__(self):
		return self
	#__iter__()
	
	
	def __next__(self):
		# if lines are still cached from the last read, pop one
		if len(self._lines) > 0:
			return self._lines.pop()
		# if there's data left in the source file, read and decompress another chunk
		if self._dc:
			data = self._dc.unused_data
			if data:
				self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
			elif not self._filePtr:
				raise Exception("cannot read a closed file")
			else:
				data = self._filePtr.read(self._chunkSize)
			if data:
				self._text += self._dc.decompress(data)
				data = None
			else:
				self._text += self._dc.flush()
				self._dc = None
		# if there's no text left, we're done
		if not self._text:
			raise StopIteration
		# split the text into lines
		self._lines = self._text.split(self._splitChar)
		self._text = ""
		# if there's more than one line, store the last to combine with the next chunk
		# (but if there's only one line, and more to read, then keep reading until we get a linebreak)
		if len(self._lines) > 1:
			self._text = self._lines.pop()
		elif self._dc:
			self._text = self._lines.pop()
			self._chunkSize *= 2
			return self.__next__()
		# reverse the remaining lines into a stack and pop one to return
		self._lines.reverse()
		return self._lines.pop()
	#__next__()
	
	
	def next(self):
		return self.__next__()
	#next()
	
	
	def seek(self, offset, whence = 0):
		if not self._filePtr:
			raise Exception("cannot seek a closed file")
		if offset != 0:
			raise Exception("zfile.seek() does not support offsets != 0")
		self._filePtr.seek(0, whence)
		self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
		self._text = ""
		self._lines = list()
	#seek()
	
	
	def close(self):
		if self._filePtr:
			self._filePtr.close()
			self._filePtr = None
	#close()
	
	
#zopen


if __name__ == "__main__":
	versMaj,versMin,versRev,versDate = 1,0,0,'2015-01-14'
	versStr = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	versDesc = "impute2-to-plink version %s" % versStr
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=versDesc,
		epilog="""
example: %(prog)s -s my.sample -i my.impute2_info.gz -g my.impute2.gz -m 0.9 -p output -b -t -d $TMPDIR
"""
	)
	parser.add_argument('-s', '--sample', type=str, metavar='file', required=True,
		help="input sample file, plain-text (required)"
	)
	parser.add_argument('-i', '--info', type=str, nargs='+', metavar='file', required=True,
		help="input marker info file(s), compressed (required)"
	)
	parser.add_argument('-g', '--genotype', type=str, nargs='+', metavar='file', required=True,
		help="input genotype file(s), compressed (required)"
	)
	parser.add_argument('-c', '--chromosome', type=str, metavar='label', default='0',
		help="chromosome to fill in for .map/.bim files (default: 0)"
	)
	parser.add_argument('-m', '--minprob', type=float, metavar='decimal', default=0.9,
		help="minimum probability to call a genotype, in decimal (default: 0.9)"
	)
	parser.add_argument('-p', '--prefix', type=str, metavar='prefix', required=True,
		help="prefix for output files (required)"
	)
	parser.add_argument('-b', '--binary', action='store_true',
		help="write binary output files (.bed/.bim/.fam)"
	)
	parser.add_argument('-t', '--text', action='store_true',
		help="write compressed plain-text output files (.ped.gz/.map.gz)"
	)
	parser.add_argument('-d', '--tempdir', type=str, metavar='directory', default='.',
		help="directory to write temporary files while transposing data for plain-text output (default: current directory)"
	)
	parser.add_argument('--version', action='version', version=versDesc)
	
	# parse arguments
	args = parser.parse_args()
	
	# store some oft-used arguments to save dictionary lookups
	args_chromosome = args.chromosome
	args_minprob = args.minprob
	args_binary = args.binary
	args_text = args.text
	
	# verify output format
	if not (args_binary or args_text):
		print "ERROR: no output formats were requested; use --binary and/or --text"
		sys.exit(2)
	
	# read sample file
	print "reading sample file '%s' ..." % args.sample
	samples = list()
	with open(args.sample,'rU') as sampleFile:
		header = sampleFile.next()
		if not header.startswith("ID_1 ID_2 missing father mother sex plink_pheno"):
			print "ERROR: unrecognized file header: %s" % header
			sys.exit(1)
		dummy = sampleFile.next()
		if not dummy.startswith("0 0 0 D D D B"):
			print "ERROR: unrecognized dummy line #2: %s" % dummy
			sys.exit(1)
		for line in sampleFile:
			words = line.rstrip().split()
			if len(words) < 7:
				print "ERROR: invalid sample on line #%d; expected 7+ columns, found %d" % (len(samples)+2, len(words))
				sys.exit(1)
			samples.append(words)
	print "... OK: %d samples" % len(samples)
	
	# read info file(s) for marker positions, frequencies and types
	markers = list()
	markerIndex = collections.defaultdict(set)
	for infoPath in args.info:
		print "reading genotype info file '%s' ..." % infoPath
		with zopen(infoPath) as infoFile:
			for line in infoFile:
				if line.startswith("#"):
					continue
				if line.startswith("snp_id rs_id position exp_freq_a1 info certainty type"):
					continue
				snpid,rsid,pos,freq,info,certainty,imptype = line.split(None,7)[:7]
				markerIndex[rsid].add( len(markers) )
				markers.append( [rsid,pos,freq,imptype] )
			#foreach line in infoFile
		#with infoFile
		print "... OK: %d markers (%d duplicate)" % (len(markers),len(markers)-len(markerIndex))
	#foreach args.info
	
	# read genotype file(s) for marker alleles
	m = 0
	for genoPath in args.genotype:
		print "reading genotype file '%s' ..." % genoPath
		with zopen(genoPath) as genoFile:
			for line in genoFile:
				if m >= len(markers):
					print "ERROR: genotype file contains too many markers"
					sys.exit(1)
				snpid,rsid,pos,a1,a2 = line.split(None,5)[:5]
				if rsid != markers[m][0]:
					print "ERROR: genotype marker #%d is '%s', expected '%s'" % (m+1,rsid,markers[m][0])
					sys.exit(1)
				markers[m].append(a1)
				markers[m].append(a2)
				m += 1
			#foreach line in genoFile
		#with genoFile
		print "... OK"
	#foreach args.genotype
	if m < len(markers):
		print "ERROR: genotype file(s) end after marker #%d, expected %d more markers" % (m,len(markers)-m)
		sys.exit(1)
	
	# markers=[ [rsid,pos,freq,type,a1,a2], ... ]
	
	# choose between duplicate markers
	markerDrop = set()
	if True:
		# the original logic was already complicated, and only covered the case of two versions;
		# in order to handle 3 versions, we moved to a simpler priority system for which ones to rename
		print "annotating duplicate markers ..."
		for marker,indecies in markerIndex.iteritems():
			if marker == ".": # special tag, rename to chr#:pos label
				for i in indecies:
					assert(marker == markers[i][0])
					markers[i][0] = "chr%s:%s" % (args_chromosome,markers[i][1])
				#foreach index
			elif len(indecies) > 1:
				typeMs = {'0':set(), '2':set(), '3':set()}
				for i in indecies:
					assert(marker == markers[i][0])
					if markers[i][3] not in typeMs:
						print "ERROR: unknown type '%s' for marker '%s' at index %d" % (markers[i][3],marker,i)
						sys.exit(1)
					typeMs[markers[i][3]].add(i)
				#foreach index
				if len(typeMs['0']) != len(set(str(markers[i][1]).strip() for i in typeMs['0'])):
					print "ERROR: duplicate type-0 positions for marker '%s' at indecies %s" % (marker,','.join(str(i) for i in typeMs['0']))
					sys.exit(1)
				if len(typeMs['2']) > 1:
					print "ERROR: multiple type-2 records for marker '%s' at indecies %s" % (marker,','.join(str(i) for i in typeMs['2']))
					sys.exit(1)
				if len(typeMs['3']) > 1:
					print "ERROR: multiple type-3 records for marker '%s' at indecies %s" % (marker,','.join(str(i) for i in typeMs['3']))
					sys.exit(1)
				if typeMs['2']:
					for i in typeMs['0']:
						markers[i][0] += "_t0"
					for i in typeMs['3']:
						markers[i][0] += "_t3"
				elif typeMs['0']:
					for i in typeMs['3']:
						markers[i][0] += "_t3"
			#if marker is special or duplicate
		#foreach marker
		print "... OK"
	else: # original dupe handling logic
		print "writing .drop.txt file '%s.drop.txt' ..." % args.prefix
		strandFlipTrans = string.maketrans('AaCcGgTt','TTGGCCAA')
		# s.translate(strandFlipTrans) is a shorter equivalent of s.replace('A','T').replace('C','G').replace('G','C').replace('T','A')
		# which is used to check if one pair of alleles is the strand-flipped complement of the other
		with open(args.prefix+'.drop.txt','wb') as dropFile:
			dropFile.write("rs_id position exp_freq_a1 type allele1 allele2 status reason\n")
			for marker,indecies in markerIndex.iteritems():
				# our logic only covers markers appearing twice; thrice or more is an error
				if len(indecies) > 2:
					print "ERROR: %d occurrences of marker '%s'; only 2 are supported" % (len(indecies),marker)
					sys.exit(1)
				if len(indecies) == 2:
					index0 = indecies.pop()
					index23 = indecies.pop()
					marker0 = markers[index0]
					marker23 = markers[index23]
					# our logic assumes the two versions of a marker will have different types (0 and non-0, usually 2 or 3)
					if marker0[3] == marker23[3]:
						print "ERROR: 2 occurrences of marker '%s' are both type %s" % (marker,marker0[3])
						sys.exit(1)
					if marker23[3] == '0':
						index0,index23 = index23,index0
						marker0,marker23 = marker23,marker0
					# if the two versions have different positions, drop the type 2/3 and keep the type 0
					if marker0[1] != marker23[1]:
						markerDrop.add(index23)
						dropFile.write("\t".join(marker0))
						dropFile.write("\t+\t\n")
						dropFile.write("\t".join(marker23))
						dropFile.write("\t-\tposition_change\n")
					# if either version has a 0 allele, drop the type 0 and keep the type 2/3
					elif marker0[4] == '0' or marker0[5] == '0' or marker23[4] == '0' or marker23[5] == '0':
						markerDrop.add(index0)
						dropFile.write("\t".join(marker0))
						dropFile.write("\t-\tfixed_allele\n")
						dropFile.write("\t".join(marker23))
						dropFile.write("\t+\t\n")
					# if both versions have the same frequency but the alleles of one are the strand-flipped complements of the other, drop the type 2/3 and keep the type 0
					elif marker0[2] == marker23[2] and sorted([marker0[4].translate(strandFlipTrans),marker0[5].translate(strandFlipTrans)]) == sorted([marker23[4].upper(),marker23[5].upper()]):
						markerDrop.add(index23)
						dropFile.write("\t".join(marker0))
						dropFile.write("\t+\t\n")
						dropFile.write("\t".join(marker23))
						dropFile.write("\t-\tstrand_flip\n")
					# if the two versions have different alleles, drop the type 0 and keep the type 2/3
					elif sorted(marker0[4:6]) != sorted(marker23[4:6]):
						markerDrop.add(index0)
						dropFile.write("\t".join(marker0))
						dropFile.write("\t-\tallele_mismatch\n")
						dropFile.write("\t".join(marker23))
						dropFile.write("\t+\t\n")
					# otherwise, keep both versions
					else:
						dropFile.write("\t".join(marker0))
						dropFile.write("\t+\tundecided\n")
						dropFile.write("\t".join(marker23))
						dropFile.write("\t+\tundecided\n")
				#if marker is duplicate
			#foreach marker
		#with dropFile
		print "... OK"
	#which dupe resolution algorithm
	markerIndex = None
	
	# (1,0,0) -> (1,1)
	# (0,1,0) -> (1,2)
	# (0,0,1) -> (2,2)
	# (0,0,0) -> (0,0)
	
	# read genotype file(s)
	sRange = range(0,len(samples))
	sWord = range(5,5+3*len(samples),3)
	m = 0
	if args_binary:
		print "writing .bed file '%s.bed' ..." % args.prefix
		bedFileKeep = open(args.prefix+'.bed','wb')
		bedFileKeep.write("\x6c\x1b") # binary plink file magic number
		bedFileKeep.write("\x01") # SNP-major order (row per snp, columns per sample)
		if markerDrop:
			print "... and .drop.bed file '%s.drop.bed' ..." % args.prefix
			bedFileDrop = open(args.prefix+'.drop.bed','wb')
			bedFileDrop.write("\x6c\x1b") # binary plink file magic number
			bedFileDrop.write("\x01") # SNP-major order (row per snp, columns per sample)
		padding = [0] * (-len(samples) % 4)
		bytes = range(0,(len(samples)+len(padding))/4)
		packformat = "%dB"%len(bytes)
	if args_text:
		mTemp = 0
		tempFilesKeep = list()
		pedOutKeep = list( [] for sample in samples )
		if markerDrop:
			tempFilesDrop = list()
			pedOutDrop = list( [] for sample in samples )
	for genoPath in args.genotype:
		print "reading genotype file '%s' ..." % genoPath
		with zopen(genoPath) as genoFile:
			for line in genoFile:
				words = line.split()
				a1 = a2 = 0
				if args.binary:
					if m in markerDrop:
						bedFile = bedFileDrop
					else:
						bedFile = bedFileKeep
				if args.text:
					if m in markerDrop:
						tempFiles = tempFilesDrop
						pedOut = pedOutDrop
					else:
						tempFiles = tempFilesKeep
						pedOut = pedOutKeep
				# some code is duplicated here for performance;
				# doing the if checks for every sample for every marker is a drag
				if args_binary and args_text:
					bedOut = list()
					geno11 = "%s %s" % (words[3],words[3])
					geno12 = "%s %s" % (words[3],words[4])
					geno22 = "%s %s" % (words[4],words[4])
					geno00 = "0 0"
					for s,w in enumerate(sWord):
						if float(words[w]) >= args_minprob:
							a1 += 2
							bedOut.append(0b00)
							pedOut[s].append(geno11)
						elif float(words[w+1]) >= args_minprob:
							a1 += 1
							a2 += 1
							bedOut.append(0b10)
							pedOut[s].append(geno12)
						elif float(words[w+2]) >= args_minprob:
							a2 += 2
							bedOut.append(0b11)
							pedOut[s].append(geno22)
						else:
							bedOut.append(0b01)
							pedOut[s].append(geno00)
					#foreach sample
				elif args_binary:
					bedOut = list()
					for s,w in enumerate(sWord):
						if float(words[w]) >= args_minprob:
							a1 += 2
							bedOut.append(0b00)
						elif float(words[w+1]) >= args_minprob:
							a1 += 1
							a2 += 1
							bedOut.append(0b10)
						elif float(words[w+2]) >= args_minprob:
							a2 += 2
							bedOut.append(0b11)
						else:
							bedOut.append(0b01)
					#foreach sample
				elif args_text:
					geno11 = "%s %s" % (words[3],words[3])
					geno12 = "%s %s" % (words[3],words[4])
					geno22 = "%s %s" % (words[4],words[4])
					geno00 = "0 0"
					for s,w in enumerate(sWord):
						if float(words[w]) >= args_minprob:
							a1 += 2
							pedOut[s].append(geno11)
						elif float(words[w+1]) >= args_minprob:
							a1 += 1
							a2 += 1
							pedOut[s].append(geno12)
						elif float(words[w+2]) >= args_minprob:
							a2 += 2
							pedOut[s].append(geno22)
						else:
							pedOut[s].append(geno00)
					#foreach sample
				#if binary/text
				if args_binary:
					bedOut.extend(padding)
				if a1 > a2:
					markers[m][2] = "%g" % (1.0-float(markers[m][2]))
					markers[m][4],markers[m][5] = markers[m][5],markers[m][4]
					if args_binary:
						for s in sRange:
							if bedOut[s] == 0b00:
								bedOut[s] = 0b11
							elif bedOut[s] == 0b11:
								bedOut[s] = 0b00
				if args_binary:
					bedFile.write(struct.pack(packformat, *((bedOut[b*4] | (bedOut[b*4+1] << 2) | (bedOut[b*4+2] << 4) | (bedOut[b*4+3] << 6)) for b in bytes)))
				if args_text:
					mTemp += 1
					if mTemp >= 100000:
						print "  writing temp file(s) #%d with %d markers ..." % (len(tempFilesKeep)+1,mTemp)
						tempFile = tempfile.TemporaryFile(mode='w+b', prefix='.impute2-to-plink.tmp%d.'%(len(tempFilesKeep)+1), dir=args.tempdir)
						tempFilesKeep.append(tempFile)
						for out in pedOutKeep:
							tempFile.write(" ".join(out))
							tempFile.write("\n")
						pedOutKeep = list( [] for sample in samples )
						if markerDrop:
							tempFile = tempfile.TemporaryFile(mode='w+b', prefix='.impute2-to-plink.tmp%d.drop.'%(len(tempFilesDrop)+1), dir=args.tempdir)
							tempFilesDrop.append(tempFile)
							for out in pedOutDrop:
								tempFile.write(" ".join(out))
								tempFile.write("\n")
							pedOutDrop = list( [] for sample in samples )
						mTemp = 0
						print "  ... OK"
					#if mTemp > x
				m += 1
			#foreach line in genoFile
		#with genoFile
		print "... OK"
	#foreach args.genotype
	if args_binary:
		bedFileKeep.close()
		if markerDrop:
			bedFileDrop.close()
	
	if args_text:
		print "writing .ped.gz file '%s.ped.gz' ..." % args.prefix
		for tempFile in tempFilesKeep:
			tempFile.seek(0)
		with gzip.GzipFile(filename=args.prefix+'.ped.gz', mode='wb', compresslevel=6) as pedFile:
			for s,sample in enumerate(samples):
				pedFile.write("%s %s %s %s %s %s" % (sample[0],sample[1],sample[3],sample[4],sample[5],sample[6]))
				for tempFile in tempFilesKeep:
					line = tempFile.next().rstrip("\n")
					if line:
						pedFile.write(" ")
						pedFile.write(line)
				if pedOutKeep[s]:
					pedFile.write(" ")
					pedFile.write(" ".join(pedOutKeep[s]))
				pedFile.write("\n")
			#foreach sample
		#with pedFile
		pedOutKeep = None
		for tempFile in tempFilesKeep:
			tempFile.close()
		print "... OK"
		
		if markerDrop:
			print "writing .drop.ped.gz file '%s.drop.ped.gz' ..." % args.prefix
			for tempFile in tempFilesDrop:
				tempFile.seek(0)
			with gzip.GzipFile(filename=args.prefix+'.drop.ped.gz', mode='wb', compresslevel=6) as pedFile:
				for s,sample in enumerate(samples):
					pedFile.write("%s %s %s %s %s %s" % (sample[0],sample[1],sample[3],sample[4],sample[5],sample[6]))
					for tempFile in tempFilesDrop:
						line = tempFile.next().rstrip("\n")
						if line:
							pedFile.write(" ")
							pedFile.write(line)
					if pedOutDrop[s]:
						pedFile.write(" ")
						pedFile.write(" ".join(pedOutDrop[s]))
					pedFile.write("\n")
				#foreach sample
			#with pedFile
			pedOutDrop = None
			for tempFile in tempFilesDrop:
				tempFile.close()
			print "... OK"
		#if markerDrop
		
		print "writing .map.gz file '%s.map.gz' ..." % args.prefix
		with gzip.GzipFile(filename=args.prefix+'.map.gz', mode='wb', compresslevel=6) as mapFile:
			for m,marker in enumerate(markers):
				if m not in markerDrop:
					mapFile.write("%s\t%s\t0\t%s\n" % (args_chromosome,marker[0],marker[1]))
		#with mapFile
		print "... OK"
		
		if markerDrop:
			print "writing .drop.map.gz file '%s.drop.map.gz' ..." % args.prefix
			with gzip.GzipFile(filename=args.prefix+'.drop.map.gz', mode='wb', compresslevel=6) as mapFile:
				for m,marker in enumerate(markers):
					if m in markerDrop:
						mapFile.write("%s\t%s\t0\t%s\n" % (args_chromosome,marker[0],marker[1]))
			#with mapFile
			print "... OK"
		#if markerDrop
	#if text
	
	if args_binary:
		print "writing .bim file '%s.bim' ..." % args.prefix
		with open(args.prefix+'.bim','wb') as bimFile:
			for m,marker in enumerate(markers):
				if m not in markerDrop:
					bimFile.write("%s\t%s\t0\t%s\t%s\t%s\n" % (args_chromosome,marker[0],marker[1],marker[4],marker[5]))
			#foreach marker
		#with bimFile
		print "... OK"
		
		if markerDrop:
			print "writing .drop.bim file '%s.drop.bim' ..." % args.prefix
			with open(args.prefix+'.drop.bim','wb') as bimFile:
				for m,marker in enumerate(markers):
					if m in markerDrop:
						bimFile.write("%s\t%s\t0\t%s\t%s\t%s\n" % (args_chromosome,marker[0],marker[1],marker[4],marker[5]))
				#foreach marker
			#with bimFile
			print "... OK"
		#if markerDrop
		
		print "writing .fam file '%s.fam' ..." % args.prefix
		with open(args.prefix+'.fam','wb') as famFile:
			for sample in samples:
				famFile.write("%s %s %s %s %s %s\n" % (sample[0],sample[1],sample[3],sample[4],sample[5],sample[6]))
		#with famFile
		print "... OK"
	#if binary
#__main__

