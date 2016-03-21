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
import os
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
		
	
def getFile(file_str):
	retfile = None
	try:
		retfile = file(file_str, 'rb')
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

def get_bgzf_block(block, block_len=(2**16-1)):
    # print("Saving %i bytes" % len(block))
    
    # if the block is too large, recursively write smaller blocks
    if len(block) > block_len:
    	return get_bgzf_block(block[:block_len], block_len) + get_bgzf_block(block[block_len:], block_len)
    
    # Giving a negative window bits means no gzip/zlib headers, -15 used in samtools
    c = zlib.compressobj(6,
                         zlib.DEFLATED,
                         -15,
                         zlib.DEF_MEM_LEVEL,
                         0)
    compressed = c.compress(block) + c.flush()
    del c
    # If we didn't compress enough, split the block in two and try again
    if len(compressed) > block_len:
    	return get_bgzf_block(filePtr, block[:len(block)/2], block_len) + get_bgzf_block(filePtr, block[len(block)/2:], block_len)

    crc = zlib.crc32(block)
    # Should cope with a mix of Python platforms...
    if crc < 0:
        crc = struct.pack("<i", crc)
    else:
        crc = struct.pack("<I", crc)
    bsize = struct.pack("<H", len(compressed) + 25)  # includes -1
    crc = struct.pack("<I", zlib.crc32(block) & 0xffffffff)
    uncompressed_length = struct.pack("<I", len(block))
    # Fixed 16 bytes,
    # gzip magic bytes (4) mod time (4),
    # gzip flag (1), os (1), extra length which is six (2),
    # sub field which is BC (2), sub field length of two (2),
    # Variable data,
    # 2 bytes: block length as BC sub field (2)
    # X bytes: the data
    # 8 bytes: crc (4), uncompressed data length (4)
    data = _bgzf_header + bsize + compressed + crc + uncompressed_length
    
    return data


def block_gzip(filePtr, data_in, block_len=(2**16-1)):
	'''
	Block gzips the input data into at most block_len characters of decompressed text
	'''
	cpos=0
	# make sure to write at least 1 block (for empty block @ end!)
	nblocks=0
	while nblocks == 0 or cpos<len(data_in):
		nblocks+=1
		epos=min(len(data_in), cpos+block_len) 
		filePtr.write(get_bgzf_block(data_in[cpos:epos], block_len))
		cpos=epos
	
	# If we're here, set the special EOF bit
	if nblocks==0:
		filePtr.write(_bgzf_eof)

	
if __name__ =="__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("-D", "--dict", help="Sequence dictionary (to determine sort order).  Can be HTTP/FTP link, DNANexus file or a local file")
	parser.add_argument("-V", "--vcf", help="The VCF file to split.  Can be either a HTTP/FTP link or a DNANexus file string or a local filename", action='append')
	parser.add_argument("-o", "--output", help="Output file", default=None)
	
	args = parser.parse_args()
	
	# make sure the intervals are sorted and non-overlapping
	#intervals.sort()
	#intervals = getNonOverlap(intervals)
	
	# OK, now I have vcf_f, idx_f, and a list of intervals
	# First, print the header of the file
	dict_f = getFile(args.dict)
	dict_idx=0;
	dict_order={}
	for l in dict_f:
		if l.startswith("@SQ"):
			for f in l.strip().split():
				if f.startswith("SN:"):
					sn = f.split(":")
					dict_order[sn[1]] = dict_idx
					dict_idx += 1

	# Now, dict_idx is a mapping of chromosome names -> position in the dictionary
	#dict_f.close()

	if args.output is not None:
		out_f = file(args.output, 'wb')
	else:
		out_f = sys.stdout
	
	# OK, I can now write to out_f as needed
	
	# get a list of the files
	vcf_list = [bgzopen(getFile(f)) for f in args.vcf]
	
	# OK, now let's read the header from the first VCF file
	hdr = []
	nl = 0
	vcfidx=0
	
	while not hdr and vcfidx < len(vcf_list):
		for l in vcf_list[vcfidx]:
			nl+=1
			hdr.append(l)
			if l.startswith("#CHROM"):
				break
		vcfidx+=1
	
	# And consume the headers for all other VCF files
	for f in vcf_list[vcfidx:]:
		for l in f:
			if l.startswith("#CHROM"):
				break
	
	# Now, let's generate a list of (chrom, basepair, position in vcf_list, <rest of chunk>) 
	# tuples for each file in vcf_list
	vcf_order=[]
	for i, f in enumerate(vcf_list):
		chrom = f.readline("\t").strip()
		# if chrom is empty, then the VCF is empty, so ignore!
		if chrom:
			pos = int(f.readline("\t").strip())
			chunk = f.readChunk()
			vcf_order.append((chrom, pos, i, chunk))
		
	# Now, sort the vcf_order
	# stable sort by chrom/pos
	vcf_order.sort(key = lambda val: (dict_order[val[0]], val[1], val[2]) )
	
	# Write the header (block gzipped, of course)
	block_gzip(out_f,'\n'.join(hdr) + "\n")
	
	# try to match a DXFile read buffer size, else fall back to 64MB
	try:
		dl_block=vcf_list[0]._read_bufsize
	except:
		dl_block=64*1024*1024
	
	# For each file, block gzip the start, then write all but the ending EOF block
	for t in vcf_order:
		
		block_gzip(out_f, t[0] + "\t" + str(t[1]) + "\t" + t[3])
		f = vcf_list[t[2]]
		cfp = f._filePtr
		# where am i?
		cvfp = f.tellvfp()
		fo, bo = convert_vfp(cvfp)
		cpos = fo
		# how big am i?
		# NOTE: this invalidates the f pointer (until next f.seek)
		cfp.seek(0, os.SEEK_END)
		csz = cfp.tell()
		# And now, f is re-validated (and cfp also moves to the correct location)
		f.seek(cvfp)
		
		
		# I need to download nbytes (size - current position - EOF block)
		nbytes = csz - cpos - len(_bgzf_eof)	
		newDat = cfp.read(min(dl_block, nbytes))

		while len(newDat) > 0:
			out_f.write(newDat)
			nbytes -= len(newDat)
			newDat = cfp.read(min(dl_block, nbytes))
		
		# OK, we're done with that file
		
	# write the ending EOF block
	out_f.write(_bgzf_eof)
	
	# and done!

