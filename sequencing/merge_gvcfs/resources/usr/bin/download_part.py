#!/usr/bin/env python

import sys
import tabix
import dxpy
import argparse
import zlib
import urllib

def zfile(filePtr, splitChar="\n", chunkSize=1*1024*1024):
        dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
        text = ""
        while dc:
                data = filePtr.read(chunkSize)
                if data:
                        text += dc.decompress(data)
                        data = None
                else:
                        text += dc.flush()
                        dc = None
                if text:
                        lines = text.split(splitChar)
                        i,x = 0,len(lines)-1
                        text = lines[x]
                        while i < x:
                                yield lines[i]
                                i += 1
                        lines = None
        #while data remains
        if text:
                yield text
#zfile()

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
	vcf_zfile = zfile(vcf_f)
	
	for l in vcf_zfile:
		if l.startswith("#"):
			print l.strip()
		else:
			break
			
	# OK, the header is printed, now time to print the intervals
	# First, actually get the tabix object
	tb = tabix.open(args.vcf, args.index)
	
	for interval in intervals:
		
		tabix_itr = ()
		if intervals[1] is None:
			tabix_itr = tb.querys(intervals[0])
		else:
			tabix_itr = tb.query(intervals[0], intervals[1][0], intervals[1][1])
		
		for line in tabix_itr:
			sys.stdout.write("\t".join(line) + "\n")

	
	
	
