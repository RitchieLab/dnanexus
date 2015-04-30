#!/usr/bin/env python

import sys
import tabix
import dxpy
import argparse
import zlib
import urllib

class zopen(object):

    def __init__(self, filePtr, splitChar="\n", chunkSize=16*1024):
        self._filePtr = filePtr
        self._splitChar = splitChar
        self._chunkSize = chunkSize
        self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
        self._text = ""
        self._lines = list()
    #__init__()


    def __del__(self):
        if self._filePtr:
            self._filePtr.close()
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
        if offset != 0:
            raise Exception("zfile.seek() does not support offsets != 0")
        self._filePtr.seek(0, whence)
        self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
        self._text = ""
        self._lines = list()
    #seek()

#zopen

def getNonOverlap(interval_in):
	
	curr_interval = None
	interval_out = []
	for interval in interval_in:
		if curr_interval is None:
			curr_interval = interval
		elif curr_interval[0] == interval[0]:
			if curr_interval[1] is None:
				pass
			elif curr_interval[1][1] >= interval[1][0]:
				curr_interval = (curr_interval[0], (curr_interval[1][0], interval[1][1]))
			else:
				interval_out.append(curr_interval)
				curr_interval = interval
		else:
			interval_out.append(curr_interval)
			curr_interval = interval
	
	interval_out.append(curr_interval)
	
	return interval_out

def readIntervals(interval_f):
	intervals = []
	for l in interval_f:
		curr_l = l.strip().split()
		try:
			intervals.append((curr_l[0], (int(curr_l[1]), int(curr_l[2])) ))
		except ValueError:
			print >> sys.stderr, "WARNING, cound not parse line '%s', ignoring" % l.strip()
	
	return intervals

def getIntervals(interval_str):
	intervals = []
	try:
		interval_f = file(interval_str, 'r')
		intervals = readIntervals(interval_f)
	except:
		try:
			interval_f = dxpy.bindings.dxfile.DXFile(interval_str, mode='r')
			intervals = readIntervals(interval_f)
		except:
			for interval in interval_str.split(","):
				try:
					curr_int = interval.split(':')
					if len(curr_int) > 1:
						curr_bound = curr_int[1].split("-")
						intervals.append((curr_int[0], (int(curr_bound[0]), int(curr_bound[1]))))
					else:
						intervals.append((curr_int[0], None))
				except:
					print >> sys.stderr, "WARNING: could not parse interval: '%s', ignoring" % interval
	
	return intervals
					
if __name__ =="__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("-L", "--intervals", help="Intervals of the file to download and print", action="append")
	parser.add_argument("-f", "--vcf", help="The VCF file to split.  Can be either a HTTP/FTP link or a DNANexus file string or a local filename")
	parser.add_argument("-i", "--index", help="The tabix index of the VCF file.  Can be HTTP/FTP link, DNANexus file string or local filename")
	
	args = parser.parse_args()
	
	intervals = []
	for iv in args.intervals:
		intervals.extend(getIntervals(iv))
	
	# make sure the intervals are sorted and non-overlapping
	intervals.sort()
	intervals = getNonOverlap(intervals)
	
	# OK, now I have vcf_f, idx_f, and a list of intervals
	# First, print the header of the file
	vcf_f = None
	try:
		vcf_f = file(args.vcf, 'r')
	except IOError:
		vcf_f = urllib.urlopen(args.vcf)
	
	# print the VCF header
	vcf_zfile = zopen(vcf_f)
	
	n_fields=None
	
	for l in vcf_zfile:
		if l.startswith("##"):
			print l.strip()
		elif l.startswith("#"):
			n_fields=len(l.strip().split())
			print l.strip()
		else:
			break
			
	# OK, the header is printed, now time to print the intervals
	# First, actually get the tabix object
	tb = tabix.open(args.vcf, args.index)
	
	for interval in intervals:
		
		tabix_itr = ()
		if interval[1] is None:
			tabix_itr = tb.querys(interval[0])
		else:
			tabix_itr = tb.query(interval[0], interval[1][0], interval[1][1])
		
		for line in tabix_itr:
			if n_fields is not None and len(line) != n_fields:
				raise IOError("Inconsistent number of fields")
				
			print "\t".join(line)
			#sys.stdout.write("\t".join(line) + "\n")	
	
