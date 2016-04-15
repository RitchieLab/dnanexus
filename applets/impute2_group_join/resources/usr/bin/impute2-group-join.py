#!/usr/bin/env python

import argparse
import collections
import gzip
import itertools
import os
import sys
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
	versMaj,versMin,versRev,versDate = 0,10,11,'2014-11-19'
	versStr = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	versDesc = "impute2-group-join version %s" % versStr
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=versDesc,
		epilog="""
example: %(prog)s -i a_chr22 b_chr22 -f my.samples -m chr22.markers -o ab_chr22 -d dupe_chr22

The script requires ~1.5 hours and ~2GB RAM per million markers to be merged,
but if resource limits are strictly enforced you should add ~500MB-1GB extra.
"""
	)
	parser.add_argument('-i', '--input', action='append', nargs='+', type=str, metavar='prefix', required=True,
		help="prefix(es) of impute2 .phased.sample and .impute2_info files to be joined (required)"
	)
	parser.add_argument('-f', '--filter', action='store', type=str, metavar='file',
		help="a file listing the sample IDs to retain (default: none)"
	)
	parser.add_argument('-m', '--markers', action='store', type=str, metavar='file',
		help="a file listing the expected order of all markers (default: none)"
	)
	parser.add_argument('-o', '--output', action='store', type=str, metavar='prefix', required=True,
		help="prefix for joined output and log files"
	)
	parser.add_argument('-d', '--dupes', action='store', type=str, metavar='prefix',
		help="prefix for duplicate sample output files"
	)
	parser.add_argument('--version', action='version', version=versDesc)
	
	# parse arguments
	args = parser.parse_args()
	
	# open input file(s)
	print "finding input files ..."
	prefixList = list(itertools.chain(*args.input))
	sampleFile = list()
	genoFile = list()
	infoFile = list()
	for prefix in prefixList:
		if os.path.exists(prefix+'.phased.sample.gz'):
			sampleFile.append(zopen(prefix+'.phased.sample.gz'))
		elif os.path.exists(prefix+'.phased.sample'):
			sampleFile.append(open(prefix+'.phased.sample','rU'))
		else:
			exit("ERROR: %s.phased.sample(.gz) not found" % prefix)
		
		if os.path.exists(prefix+'.impute2.gz'):
			genoFile.append(zopen(prefix+'.impute2.gz'))
		elif os.path.exists(prefix+'.impute2'):
			genoFile.append(open(prefix+'.impute2','rU'))
		else:
			exit("ERROR: %s.impute2(.gz) not found" % prefix)
		
		if os.path.exists(prefix+'.impute2_info.gz'):
			infoFile.append(zopen(prefix+'.impute2_info.gz'))
		elif os.path.exists(prefix+'.impute2_info'):
			infoFile.append(open(prefix+'.impute2_info','rU'))
		else:
			exit("ERROR: %s.impute2_info(.gz) not found" % prefix)
		
		print "  #%d: %s" % (len(sampleFile),prefix)
	#foreach prefixList
	iRange0 = range(0,len(sampleFile))
	iRange1 = range(1,len(sampleFile))
	print "... OK: %d sets of input files" % len(sampleFile)
	
	# read the marker file, if any
	markerIndex = collections.OrderedDict()
	if args.markers and ((args.markers == '-') or os.path.exists(args.markers)):
		print "reading markers file ..."
		with (sys.stdin if (args.markers == '-') else open(args.markers,'rU')) as markerFile:
			for line in markerFile:
				if not line.startswith('#'):
					words = line.rstrip("\r\n").split()
					marker = (words[2].lower(), min(words[3],words[4]).lower(), max(words[3],words[4]).lower())
					if marker in markerIndex:
						exit("ERROR: duplicate marker: %s" % (" ".join(words),))
					markerIndex[marker] = (len(markerIndex),words[1])
		#with markerFile
		print "... OK: %d markers" % (len(markerIndex),)
	else:
		print "building markers index from input files ..."
		markerList = list()
		markerOrder = dict()
		markerExpect = set()
		markerGeno = dict()
		markerInfo = dict()
		markerLabels = collections.defaultdict(set)
		markerDupe = collections.defaultdict(set)
		markerSwap = collections.defaultdict( lambda:collections.defaultdict(set) ) # {m1:{m2:{i}}}
		
		for i in iRange0:
			m = mCur = mPrev = 0
			header = infoFile[i].next().rstrip("\r\n")
			while header.startswith('#'):
				header = header[1:]
			if header != "snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0":
				exit("ERROR: invalid header on info file #%d: %s" % (i+1,header))
			markerExpect = set(markerOrder)
			while True:
				# make sure geno/info files agree with eachother
				try:
					geno = genoFile[i].next().rstrip("\r\n").split(None,5)[:-1]
				except StopIteration:
					try:
						line = "#"
						while line.startswith("#") or (line == header):
							line = infoFile[i].next().rstrip("\r\n")
						exit("ERROR: input genotype file #%d ended after %d markers, but info file continues" % (i+1,m))
					except StopIteration:
						break
				try:
					line = "#"
					while line.startswith("#") or (line == header):
						line = infoFile[i].next().rstrip("\r\n")
					info = line.split()
				except StopIteration:
					exit("ERROR: input info file #%d ended after %d markers, but genotype file continues" % (i+1,m))
				m += 1
				if (geno[1].lower() != info[1].lower()) or (geno[2].lower() != info[2].lower()):
					exit("ERROR: marker #%d mismatch in input files #%d: '%s %s' vs '%s %s'" % (m,i+1,geno[1],geno[2],info[1],info[2]))
				marker = (geno[2].lower(), min(geno[3],geno[4]).lower(), max(geno[3],geno[4]).lower())
				
				if i == 0:
					# for the first input, just check for duplicates and store metadata
					if marker in markerOrder:
						markerDupe[marker].add(i)
					else:
						markerOrder[marker] = len(markerList)
						markerList.append(marker)
						markerGeno[marker] = geno
						markerInfo[marker] = info
						for lbl in geno[1].lower().split(';'):
							markerLabels[marker].add(lbl)
				elif marker in markerExpect:
					# for subsequent inputs, verify relative order
					for lbl in geno[1].lower().split(';'):
						markerLabels[marker].add(lbl)
					mCur = markerOrder[marker]
					if mCur >= mPrev:
						markerExpect.remove(marker)
						mPrev = mCur
					while mCur < mPrev:
						mCur += 1
						markerSwap[marker][markerList[mCur]].add(i+1) # +1 here so we can .join() later
				elif marker in markerOrder:
					markerDupe[marker].add(i)
				#if i
			#while next()
			genoFile[i].seek(0)
			infoFile[i].seek(0)
			if i == 0:
				print "  #%d: %d markers" % (i+1,len(markerOrder))
			else:
				for marker in markerExpect:
					del markerOrder[marker]
				print "  #%d: %d markers (%d matching)" % (i+1,m,len(markerOrder))
		#foreach input
		
		# apply highest RS# labels to all markers
		for marker,labels in markerLabels.iteritems():
			rses = set(int(l[2:]) for l in labels if l.startswith('rs'))
			markerGeno[marker][1] = ('rs%d' % max(rses)) if rses else min(labels)
		
		# check for marker dupe warnings
		if markerDupe:
			print "WARNING: %d markers are duplicated in one or more .impute2(.gz) files" % (len(markerDupe),)
			if args.dupes:
				print "writing duplicate markers to '%s' ..." % (args.dupes+'.markers',)
				with open(args.dupes+'.markers','wb') as dupesFile:
					for marker in sorted(markerDupe, key=markerOrder.get):
						dupesFile.write("%s %s\n" % (" ".join(markerGeno[marker]), " ".join(prefixList[i] for i in sorted(markerDupe[marker]))))
				print "... OK"
		
		# check for marker swap errors
		iSwaps = None
		for marker1 in markerSwap:
			for marker2 in markerSwap[marker1]:
				if (marker1 in markerOrder) and (marker2 in markerOrder):
					iSwaps = [str(iSwap) for iSwap in sorted(markerSwap[marker1][marker2])]
					print "ERROR: marker positions %s and %s order swapped in .impute2(.gz) file(s) #%s" % (marker1[0],marker2[0],",#".join(iSwaps))
		if iSwaps:
			exit(1)
		
		# compile matched markers
		mPrev = 0
		for marker in markerList:
			if (marker in markerOrder) and (markerOrder[marker] >= mPrev):
				if marker in markerIndex:
					exit("ERROR: duplicate marker: %s" % (" ".join(marker),))
				markerIndex[marker] = (len(markerIndex),markerGeno[marker][1])
				mPrev = markerOrder[marker]
		print "... OK: %d matched markers" % (len(markerIndex),)
		
		# write final marker index to file
		if args.markers:
			print "writing markers index file ..."
			with open(args.markers,'wb') as markerFile:
				markerFile.write("\n".join( " ".join(markerGeno[marker]) for marker in markerIndex.iterkeys() ))
				markerFile.write("\n")
			print "... OK: %d markers written" % (len(markerIndex),)
		markerList = markerOrder = markerExpect = markerGeno = markerInfo = markerLabels = markerDupe = markerSwap = None
	#if args.markers
	
	# read the sample filter file, if any
	sampleFilter = None
	if args.filter:
		sampleFilter = set()
		print "reading sample filter file ..."
		with (sys.stdin if (args.filter == '-') else open(args.filter,'rU')) as filterFile:
			for line in filterFile:
				if not line.startswith('#'):
					sample = tuple(line.rstrip("\r\n").lower().split(None,2))[0:2]
					if sample in sampleFilter:
						exit("ERROR: duplicate sample %s" % sample)
					sampleFilter.add(sample)
		#with filterFile
		print "... OK: %d samples" % (len(sampleFilter),)
	
	# initialize buffers
	sampleOut = open(args.output+'.phased.sample', 'wb')
	sampleDupe = None
	genoOut = gzip.open(args.output+'.impute2.gz', 'wb', compresslevel=6)
	genoDupe = None
	genoCols = [ None for i in iRange0 ]
	genoUniq = [ list() for i in iRange0 ]
	genoLine = [ None for i in iRange0 ]
	genoMarker = [ None for i in iRange0 ]
	genoSkip = [ 0 for i in iRange0 ]
	infoOut = gzip.open(args.output+'.impute2_info.gz', 'wb', compresslevel=6)
	infoOut.write("snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0\n")
	infoLine = [ None for i in iRange0 ]
	logOut = open(args.output+'.log', 'wb')
	logOut.write("#snp_id\trs_id\tposition\tallele1\tallele2\tstatus\tnote\n")
	
	# check samples
	print "joining samples ..."
	sampleHeader1 = sampleHeader2 = None
	sampleFirst = dict()
	sampleDupes = list()
	inputDrop = set()
	for i in iRange0:
		line = sampleFile[i].next().rstrip("\r\n")
		if not line.startswith("ID_1 ID_2 missing"):
			exit("ERROR: invalid header in .sample input file #%d: %s" % (i+1,line))
		if not sampleHeader1:
			sampleHeader1 = line
			sampleOut.write("%s\n" % sampleHeader1)
		elif line != sampleHeader1:
			exit("ERROR: mismatched header in .sample input file #%d: %s" % (i+1,line))
		
		line = sampleFile[i].next().rstrip("\r\n")
		if not line.startswith("0 0 0"):
			exit("ERROR: invalid subheader in .sample input file #%d: %s" % (i+1,line))
		if not sampleHeader2:
			sampleHeader2 = line
			sampleOut.write("%s\n" % sampleHeader2)
		elif line != sampleHeader2:
			exit("ERROR: mismatched subheader in .sample input file #%d: %s" % (i+1,line))
		
		samples = list()
		for line in sampleFile[i]:
			samples.append(tuple(line.rstrip("\r\n").split()))
		sampleFile[i].close()
		
		# identify duplicate samples
		numFilter = numDupe = 0
		for s,sample in enumerate(samples):
			sampleID = (sample[0].lower(),sample[1].lower())
			if sampleFilter and (sampleID not in sampleFilter):
				numFilter += 1
			elif sampleID in sampleFirst:
				numDupe += 1
				sampleDupes.append(sampleFirst[sampleID]+(i,s))
				if args.dupes:
					if not sampleDupe:
						sampleDupe = open(args.dupes+'.phased.sample', 'wb')
						sampleDupe.write("%s\n" % sampleHeader1)
						sampleDupe.write("%s\n" % sampleHeader2)
					if not genoDupe:
						genoDupe = open(args.dupes+'.impute2.gz', 'wb')
					sampleDupe.write("(%d/%d)%s\n" % (sampleFirst[sampleID][0],i,(" ".join(sample))))
			else:
				sampleFirst[sampleID] = (i,s)
				genoUniq[i].extend(xrange(5+s*3,5+s*3+3))
				sampleOut.write("%s\n" % (" ".join(sample),))
		#foreach samples
		
		# store expected column count
		print "  #%d: %d samples (%d unique, %d duplicate, %d filtered)" % (i+1,len(samples),len(genoUniq[i])/3,numDupe,numFilter)
		genoCols[i] = 5 + len(samples) * 3
		
		# if any samples will be dropped from this input due to filter or dupe, check if *all* samples were dropped;
		if numFilter or numDupe:
			if not genoUniq[i]:
				genoUniq[i] = False
		else:
			genoUniq[i] = True
	#foreach input
	print "... OK: %d unique samples, %d duplicates" % (len(sampleFirst),len(sampleDupes))
	if sampleDupes and not args.dupes:
		print "WARNING: no duplicate output prefix was specified; duplicate samples will be silently dropped!"
	
	# join lines
	print "joining data ..."
	nextPctP = 10
	nextPctM = int(len(markerIndex) * (nextPctP / 100.0))
	numMatch = 0
	numSkip = 0
	markerSkip = set()
	try:
		# validate info headers
		for i in iRange0:
			header = infoFile[i].next().rstrip("\r\n")
			while header.startswith('#'):
				header = header[1:]
			if header != "snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0":
				exit("ERROR: invalid header on info file #%d: %s" % (i+1,header))
		
		# join each marker in index order
		m = 0
		for marker,indexlabel in markerIndex.iteritems():
			index,label = indexlabel
			if m > nextPctM:
				print "  ... %d%% ..." % nextPctP
				nextPctP += 10
				nextPctM = int(len(markerIndex) * (nextPctP / 100.0))
			m += 1
			
			# try to read forward to this marker in all inputs
			for i in iRange0:
				genoMarker[i] = None
			match = True
			for i in iRange0:
				while markerIndex.get(genoMarker[i],(-1,None))[0] < index:
					if genoMarker[i]:
						if genoMarker[i] not in markerSkip:
							markerSkip.add(genoMarker[i])
							logOut.write("%s\t%s\t%s\t-\tnot matched\n" % (genoLine[i][0], genoLine[i][1], "\t".join(genoMarker[i])))
						genoSkip[i] += 1
					genoLine[i] = line = genoFile[i].next().rstrip("\r\n").split()
					genoMarker[i] = (line[2].lower(), min(line[3],line[4]).lower(), max(line[3],line[4]).lower())
					line = "#"
					while line.startswith("#") or (line == header):
						line = infoFile[i].next().rstrip("\r\n")
					infoLine[i] = line.split()
					if (genoLine[i][1].lower() != infoLine[i][1].lower()) or (genoLine[i][2].lower() != infoLine[i][2].lower()):
						exit("ERROR: marker #%d mismatch in input files #%d: '%s %s' vs '%s %s'" % (1,i+1,genoLine[i][1],genoLine[i][2],infoLine[i][1],infoLine[i][2]))
				match = match and (genoMarker[i] == marker)
			#foreach input
			
			# if the expected marker wasn't found in all inputs, move on to the next
			if not match:
				numSkip += 1
				continue
			numMatch += 1
			
			# extract marker details, but use the preferred label
			snp = genoLine[0][0]
			pos = genoLine[0][2]
			a1 = genoLine[0][3]
			a2 = genoLine[0][4]
			aliases = set(genoLine[i][1].lower() for i in iRange0 if genoLine[i][1].lower() != label.lower())
			if aliases:
				logOut.write("%s\t%s\t%s\t%s\t%s\t+\t%s\n" % (snp,label,pos,a1,a2,";".join(sorted(aliases))))
			genoLine[0][1] = label
			
			# validate column counts
			for i in iRange0:
				if len(genoLine[i]) != genoCols[i]:
					exit("ERROR: expected %d columns in input .impute2.gz file #%d, but found %d for marker '%s'" % (genoCols[i],i+1,len(genoLine[i]),label))
			#foreach input
			
			# for the first input, store the allele order and then write the data through directly
			values = list()
			if genoUniq[0] == True:
				genoOut.write(" ".join(genoLine[0]))
			elif genoUniq[0] != False:
				genoOut.write("%s %s %s %s %s " % (snp,label,pos,a1,a2))
				genoOut.write(" ".join(genoLine[0][c] for c in genoUniq[0]))
			else:
				genoOut.write("%s %s %s %s %s" % (snp,label,pos,a1,a2))
			if infoLine[0][3] != "-1":
				values.append(float(infoLine[0][3]))
			
			# for other inputs, compare allele order to input 1
			for i in iRange1:
				if genoLine[i][3] == a2 and genoLine[i][4] == a1:
					# swap all the probabilities
					genoLine[i][3] = a1
					genoLine[i][4] = a2
					for c in xrange(5,len(genoLine[i]),3):
						genoLine[i][c],genoLine[i][c+2] = genoLine[i][c+2],genoLine[i][c]
					if infoLine[i][3] != "-1":
						values.append(1.0 - float(infoLine[i][3]))
					print "  WARNING: swapped allele order for .impute2(.gz) #%d marker '%s'" % (i+1,label)
				elif genoLine[i][3] != a1 or genoLine[i][4] != a2:
					exit("ERROR: .impute2(.gz) #%d marker '%s' allele mismatch (%s/%s expected, %s/%s found)" % (i+1,label,a1,a2,genoLine[i][3],genoLine[i][4]))
				elif infoLine[i][3] != "-1":
					values.append(float(infoLine[i][3]))
				if genoUniq[i] == True:
					genoOut.write(" ")
					genoOut.write(" ".join(genoLine[i][5:]))
				elif genoUniq[i] != False:
					genoOut.write(" ")
					genoOut.write(" ".join(genoLine[i][c] for c in genoUniq[i]))
			#foreach input
			genoOut.write("\n")
			
			# merge info data (snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0)
			infoOut.write("%s %s %s %s" % (snp,label,pos,("%1.3f" % (sum(values)/len(values))) if values else "-1"))
			for c in (4,5):
				values = list(float(infoLine[i][c]) for i in iRange0 if infoLine[i][c] != "-1")
				infoOut.write(" %s" % (("%1.3f" % (sum(values)/len(values))) if values else "-1",))
			values = list(int(infoLine[i][6]) for i in iRange0)
			infoOut.write(" %d" % (min(values),))
			for c in (7,8,9):
				values = list(float(infoLine[i][c]) for i in iRange0 if infoLine[i][c] != "-1")
				infoOut.write(" %s" % (("%1.3f" % (sum(values)/len(values))) if values else "-1",))
			infoOut.write("\n")
			
			# write dupe lines from various inputs, if any
			if sampleDupes and args.dupes:
				genoDupe.write("%s %s %s %s %s " % (snp,label,pos,a1,a2))
				genoDupe.write(" ".join(("%s %s %s" % tuple(genoLine[dupe[0]][(5+3*dupe[1]):(8+3*dupe[1])])) for dupe in sampleDupes))
				genoDupe.write("\n%s %s %s %s %s " % (snp,label,pos,a1,a2))
				genoDupe.write(" ".join(("%s %s %s" % tuple(genoLine[dupe[2]][(5+3*dupe[3]):(8+3*dupe[3])])) for dupe in sampleDupes))
				genoDupe.write("\n")
			#if dupes
		#foreach marker
	except StopIteration:
		pass
	genoOut.close()
	infoOut.close()
	logOut.close()
	if genoDupe:
		genoDupe.close()
	for i in iRange0:
		if genoSkip[i] > 0:
			print "  WARNING: input .impute2(.gz) file #%d had %d extra markers skipped during processing" % (i+1,genoSkip[i])
		n = 0
		try:
			while True:
				genoFile[i].next()
				n += 1
		except StopIteration:
			pass
		genoFile[i].close()
		infoFile[i].close()
		if n > 0:
			print "  WARNING: input .impute2(.gz) file #%d has %d leftover lines" % (i+1,n)
	#foreach input
	print "... OK: joined %d markers (%d matched, %d incomplete)" % (len(markerIndex),numMatch,numSkip)
#__main__
