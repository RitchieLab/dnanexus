#!/usr/bin/env python

import sys
import dxpy
import argparse
import zlib
import urllib
import cStringIO
import struct
import gzip
import bisect
from collections import deque

_bgzf_magic = b"\x1f\x8b\x08\x04"
_bgzf_header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
_bytes_BC = b"BC"


def convert_vfp(vfp):
	'''
	convert the virtual file pointer into a tuple of 
	(block offset, within block offset)
	'''
	return (vfp >> 16, vfp & (2**16 -1))
	
def convert_offsets(fo, bo):
	'''
	Convert a pair of offsets to a single virtual file pointer
	'''
	return (fo << 16) + bo

def reg2bins(rbeg, rend=None):
	'''
	Copied directly from tabix spec
	'''
	if rend is None:
		rend = rbeg
	else:
		rend = rend-1
	binlist = []
	
	for shift, offset in [(26,1),(23,9),(20,73),(17,585),(14,4681)]:
		k= offset+(rbeg>>shift)
		while k<= offset + (rend>>shift):
			binlist.append(k)
			k+=1
	
	return binlist

class bgzopen(object):

	def __init__(self, filePtr, splitChar=None, chunkSize=64*1024):
		'''
		NOTE: The chunkSize should be the maximum size of a bgzip block
		'''
	
		self._filePtr = filePtr
		self._foq = deque([])
		self._boq = deque([])
		self._fo = filePtr.tell()
		self._bo = 0
		self._dcFilePtr = cStringIO.StringIO()
		self._splitChar = splitChar
		
		self._delimChar = self._splitChar
		if self._delimChar is None:
			self._delimChar = '\n\r'
			
		
		self._chunkSize = chunkSize
		self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
		self._text = ""
		# do not store the lines, use filePtr.readline() instead
		# self._lines = list()
	#__init__()


	def __del__(self):
		# DO NOTHING to the file pointer, please!!
		pass
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
		# read a line, stripping the newline
		return self.readline().rstrip(self._delimChar)
	#__next__()
	
	def _cacheChunk(self):
	
		newText = ''
		newLen = 0
		procDataLen = 1
		
		# While I've downloaded additional data and haven't added anything to 
		# my cached text, keep trying
		while newLen == 0 and procDataLen > 0:
			
			data = None
			
			# If this is true, we likely have a partial block
			if self._dc and self._dc.unused_data:
				data = self._dc.unused_data
				#print >> sys.stderr, "looking at unused data", repr(data[:10])
			else:
				#print >> sys.stderr, "current ptr:", self._filePtr.tell()
				data = self._filePtr.read(self._chunkSize)
			
			datalen=len(data)
			if data:
				self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
								
			if len(self._foq) == 0:							
				# Append the file offset to the file offset queue
				self._foq.append(self._filePtr.tell() - len(data)) 
				
			# If we're here and we have no data, we've likely hit the end of 
			# the compressed file
			if data:
				try:
					newText = self._dc.decompress(data)
				except:
					print >> sys.stderr, "ERROR decompressing:"
					#print >> sys.stderr, repr(data[:50])
					raise
				
				dlMore = True
				while len(self._dc.unused_data) == 0 and dlMore:
					newDat = self._filePtr.read(self._chunkSize)
					dlMore = len(newDat) > 0
					datalen += len(newDat)
					newText += self._dc.decompress(newDat)
								
				procDataLen = datalen - len(self._dc.unused_data)
				newLen = len(newText)
				self._boq.append(newLen)
				self._foq.append(self._foq[-1] + procDataLen)
				self._text += newText
				data = None
				#print >> sys.stderr, "decompressing!, added", newLen, "bytes from", procDataLen, "data"
			elif self._dc:
				newText = self._dc.flush()
				newLen = len(newText)
				procDataLen = datalen
				self._boq.append(newLen)
				self._foq.append(self._foq[-1] + procDataLen)
				self._text += newText
				self._dc = None
				#print "No data to be had, flushed and got", newLen, "bytes from", procDataLen, "data"
			else:
				# break out with no data read!
				# also kill the _foq added previously
				self._foq.popleft()
				newText = ""
				newLen = 0
				procDataLen = 0
	
		# end while loop
		
		#print >> sys.stderr, "Cached chunk"
		#print >> sys.stderr, "FOQ:", self._foq
		#print >> sys.stderr, "BOQ:", self._boq
				
		# return the # of bytes decompressed
		return newLen
		
	def _returnString(self, strlen):
		retlen = min(len(self._text), strlen)
		retstr = self._text[:retlen]
		self._text = self._text[retlen:]
		
		##print >> sys.stderr, "returning string of length", retlen
		
		# add the block offset
		#self._bo += retlen
		popped = False
		
		while len(self._boq) > 0 and self._bo + retlen >= self._boq[0]:
			retlen -= self._boq.popleft() - self._bo
			self._bo = 0	
			self._fo = self._foq.popleft()
			if len(self._foq) > 0:
				self._fo = self._foq[0]
		
		self._bo += retlen
		
		return retstr
	
	def tellvfp(self):
		'''
		Gives the current virtual file pointer of the current stream
		'''
		return convert_offsets(self._fo, self._bo)
	
	def seek(self, vfp):
		'''
		Seeks to the given virtual file pointer by seeking to the beginning of
		the block and reading the given number of characters
		'''
		# If we're seeking, we need to throw away our decompression object
		self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
		
		#print >> sys.stderr, "Seeking to", convert_vfp(vfp)
		
		# clear the foq and boq
		self._foq.clear()
		self._boq.clear()
		
		(fo, bo) = convert_vfp(vfp)
		if self._filePtr.tell() != fo:
			self._filePtr.seek(fo)
		else:
			pass
			#print >> sys.stderr, "No need to seek!"

			#print >> sys.stderr, "let's peek, though:"
			#print >> sys.stderr, repr(self._filePtr.read(10))
			#self._filePtr.seek(fo)
		
		self._fo = fo
		self._bo = 0
		self._text=""
		
		txt = self.read(bo)
		#print >> sys.stderr, "Seeking, throwing away:", (len(txt) == bo)
		#print >> sys.stderr, txt
		
		# return "True" if successful
		return (len(txt) == bo)
		
	def readUntil(self, vfp):
		'''
		Reads until the given virtual file pointer offset
		'''
		fo, bo = convert_vfp(vfp)
		text = ''
		while fo > self._fo:
			
			text += self.readChunk()
		
		# OK, now we're in the correct block, so read only 'bo' bytes more
		text += self.read(bo)
		
		return text
		
	def readChunk(self):
		'''
		Reads until the end of the chunk and no more
		Helpful when reading bgzipped data
		'''
		if len(self._text) == 0:
			self._cacheChunk()
			
		
		toret = self._returnString(self._boq[0] - self._bo)
		#self._cacheChunk()
		return toret
		
	def read(self, nbytes):
		'''
		Read at most nbytes, decompressed
		'''
		ret_text = ""
		addlDataLen = 1
		# Note, _cacheChunk() will change the length of self._text
		while len(self._text) + len(ret_text) < nbytes and addlDataLen > 0:
			if len(self._text) != 0:
				# read until end of chunk, clearing out self._text
				ret_text += self.readChunk()
			# now, fill self._text with a new block
			addlDataLen = self._cacheChunk()
			#print "cached an additional", addlDataLen, "chars"

		return ret_text + self._returnString(nbytes-len(ret_text))
		
	def readline(self, breakChar=None, nbytes=-1):
		'''
		Read until a cahracter in 'breakChars' is found.  Return the string
		read, including the character
		
		Note: if breakChars is None, read until any newline, so check to make
		sure if we find '\r' to also check for '\n' immediately following
		'''
		if breakChar is None:
			breakChar=self._delimChar

		firstpos = min([len(self._text)] + [i for i in (self._text.find(c) for c in breakChar) if i >= 0])
		
		while firstpos==len(self._text) and self._cacheChunk() > 0:
			firstpos = min([len(self._text)] + [i for i in (self._text.find(c) for c in breakChar) if i >= 0])
		
		if breakChar == self._delimChar and self._splitChar is None and firstpos < len(self._text) and self._text[firstpos] == '\r':
			# if our last character was carriage return, read another chunk
			if firstpos == len(self._text) - 1:
				self._cacheChunk() > 0
				
			if firstpos +1 < len(self._text) and self._text[firstpos+1] == '\n':
				firstpos+=1
		
		
		pbo = self._bo
		#print >> sys.stderr, "Current bo:", self._bo
		txt = self._returnString(firstpos+1)
		#print >> sys.stderr, "New bo:", self._bo
		
		#print >>sys.stderr, "Diff and len:", self._bo - pbo, len(txt)
		
		return txt

	def next(self):
		return self.__next__()
	#next()
#zopen

class tbi_bin(object):

	def __init__(self, filePtr):
		self.bin = struct.unpack('<L', filePtr.read(4))[0]
		self.n_chunk = struct.unpack('<l', filePtr.read(4))[0]
		self.chunks = [(struct.unpack('<Q', filePtr.read(8))[0], struct.unpack('<Q', filePtr.read(8))[0] ) for i in range(self.n_chunk)]

class tbi_ref(object):

	def __init__(self, filePtr):
		self.n_bin = struct.unpack('<l', filePtr.read(4))[0]
		self.bins = [tbi_bin(filePtr) for i in range(self.n_bin)]
		self.bin_map = {b.bin: b for b in self.bins}
		self.n_intv = struct.unpack('<l', filePtr.read(4))[0]
		self.ioff = [struct.unpack('<Q', filePtr.read(8))[0] for i in range(self.n_intv)]	
		self.first_pos = self.ioff[next((i for i, x in enumerate(self.ioff) if x), None)]
		
	def getBins(self, rbeg, rend=None):
		'''
		Gets all the bins actually in this reference data that overlap the region
		'''
		bin_opts = reg2bins(rbeg, rend)
		return [z for z in [self.bin_map.get(b, None) for b in bin_opts] if z is not None]
		
class tbi_data(object):

	def __init__(self, filePtr):
		'''
		Initialize the tabix index
		'''
		self.magic = filePtr.read(4)
		if self.magic != "TBI\1":
			raise IOError("Not a valid TBI File")

		self.n_ref = struct.unpack('<l', filePtr.read(4))[0]
		self.format = struct.unpack('<l', filePtr.read(4))[0]
		self.col_seq = struct.unpack('<l', filePtr.read(4))[0]
		self.col_beg = struct.unpack('<l', filePtr.read(4))[0]
		self.col_end = struct.unpack('<l', filePtr.read(4))[0]
		self.meta = struct.unpack('<l', filePtr.read(4))[0]
		self.skip = struct.unpack('<l', filePtr.read(4))[0]
		self.l_nm = struct.unpack('<l', filePtr.read(4))[0]
		
		name_str = filePtr.read(self.l_nm)
		self.names = name_str.rstrip('\0').split('\0')
		self.name_map = {v : i for i,v in enumerate(self.names)}
		
		self.ref = [tbi_ref(filePtr) for i in range(self.n_ref)]
		
		self.n_no_coor = None
		qword = filePtr.read(8)
		if len(qword) == 8:
			self.n_no_coor = struct.unpack('<Q', qword)[0]
			
		self.linearIndex=[0] + [i for o in self.ref for i in o.ioff if i>0]
		self.linearIndex.sort()

	def getNextChrom(self, vfp):
		'''
		Gets the offset of the first record in the next chromosome after the 
		given offset.  Returns None if the given offset is in the final 
		chromosome
		'''
		gt_offsets=[r.first_pos for r in self.ref if r.first_pos>vfp]
		if len(gt_offsets) == 0:
			return None
		else:
			return min(gt_offsets)
		
def getFile(file_str):
	retfile = None
	try:
		retfile = file(file_str, 'r')
	except IOError:
		try:
			retfile = urllib.urlopen(file_str)
		except IOError:
			try:
				retfile = dxpy.bindings.dxfile.DXFile(file_str, mode='r')
			except:
				# OK, I give up, just return nothing!
				pass
		
	return retfile


if __name__ =="__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("-L", "--interval", help="Single interval of the file to download")
	parser.add_argument("-H", "--header", help="Print VCF header", action="store_true", default=False)
	parser.add_argument("-i", "--index", help="The tabix index of the VCF file.  Can be HTTP/FTP link, DNANexus file string or local filename")
	
	args = parser.parse_args()
	
	intv_chr_reg=args.interval.split(':')
	intv_chr = intv_chr_reg[0]
	intv_start = None
	intv_end = None
	if len(intv_chr_reg) > 1:
		intv_bounds=intv_chr_reg[1].split('-')
		intv_start=int(intv_bounds[0])
		intv_end=intv_start
		if len(intv_bounds) > 1:
			intv_end=int(intv_bounds[1])
		
		intv_start, intv_end = sorted([intv_start, intv_end])
			
	# OK, we now have intv_chr, intv_start and intv_end
	
	# OK, now I have idx_f, and a list of intervals
	# First, print the header of the file
	vcfidx_f = getFile(args.index)
	vcfidx_zf=bgzopen(vcfidx_f)
	vcfidx_data = tbi_data(vcfidx_zf)
	
	# there's a 28-byte EOF block always
	num_bytes=28
	
	# So, we're going to do this with the vcf_zfile and vcfidx_data objects
	# if we want to print the header, 
	if args.header:

		# get the virtual fp of the first position, which gives us the size of the header	
		start_vfp = vcfidx_data.getNextChrom(0) 
		start_fo, start_bo = convert_vfp(start_vfp)
		
		# add the compressed portion + uncompressed portion.  unless something
		# is REALLY wacky about the last block of the header, the uncompressed
		# portion should compress to a smaller size.  We want an estimate, but an 
		# upper bound is most important
		num_bytes += start_fo + start_bo
		

	if intv_chr not in vcfidx_data.name_map:
		print >> sys.stderr, "ERROR: Chromosome not in target VCF, exiting"
		print num_bytes
		sys.exit(0)
		
	chrom_ref = vcfidx_data.ref[vcfidx_data.name_map[intv_chr]]
	
	if intv_start is None:
		# we want a whole chromosome
		data_start_vfp = chrom_ref.first_pos
	else:
		# get just the linear index?
		st_index = intv_start/(16*1024)
		if st_index >= len(chrom_ref.ioff):
			print >> sys.stderr, "WARNING: start position not in linear index; interval not in VCF?"
			data_start_vfp = chrom_ref.ioff[-1]
		else:
			print >> sys.stderr, convert_vfp(chrom_ref.ioff[st_index])
			data_start_vfp = max(chrom_ref.first_pos, chrom_ref.ioff[st_index])

	# if data_start_vfp is None, then we can't find the interval AT ALL!
	if data_start_vfp is None:
		print >> sys.stderr, "ERROR: Interval not in target VCF, exiting"
		print num_bytes
		sys.exit(0)
		
	
	# Also, let's get the start position of the bin containing the end of our
	# interval
	data_end_vfp=None
	if intv_start is None:
		data_end_vfp = vcfidx_data.getNextChrom(data_start_vfp)
	elif data_start_vfp is not None:
		end_index = intv_end/(16*1024)
		if end_index < len(chrom_ref.ioff) + 1:
			data_end_vfp = max(chrom_ref.first_pos, chrom_ref.ioff[end_index+1])
			print >> sys.stderr, convert_vfp(data_end_vfp)
		else:
			data_end_vfp = vcfidx_data.getNextChrom(data_start_vfp)
		
	# this told us "until EOF", so let's find a reasonable approximation for the 'n+1'-th file offset
	if data_end_vfp is None:
		data_end_vfp = vcfidx_data.linearIndex[-1]
		if len(vcfidx_data.linearIndex) > 2:
			data_next_vfp = vcfidx_data.linearIndex[-2]
			dn_fo, dn_bo = convert_vfp(data_next_vfp)
			de_fo, de_bo = convert_vfp(data_end_vfp)
				
			# and create a "phony" vfp
			data_end_vfp = convert_offsets(de_fo + (de_fo - dn_fo), 0)
			
	ds_fo, ds_bo = convert_vfp(data_start_vfp)
	de_fo, de_bo = convert_vfp(data_end_vfp)
	
	# return the block offset of the start of the beginning 16KB block to the 
	# block offset of the end of the ending 16KB block
	num_bytes += (de_fo - ds_fo)
	
	print num_bytes
