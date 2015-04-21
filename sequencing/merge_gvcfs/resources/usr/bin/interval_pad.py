#!/usr/bin/env python

# This reads a bed file from standard input, pads the intervals, sorts and then
# outputs the intervals guranteed to be non-overlapping

import sys

if __name__ == "__main__":
	interval_pad = 0
	if len(sys.argv) > 1:
		try:
			interval_pad = int(sys.argv[1])
		except ValueError:
			print >> sys.stderr, "WARNING: argument", sys.argv[1], "not convertible to int, setting interval padding to 0."
	
	intervals = []
	for l in sys.stdin:
		items = l.strip().split()
		try:
			intervals.append((items[0], (int(items[1]) - interval_pad, int(items[2]) + interval_pad)))
		except ValueError:
			print >> sys.stderr, "Warning: unable to parse line '%s', ignoring" % l.strip()
	
	intervals.sort()
	
	curr_int = None
	for i in intervals:
		if curr_int is None:
			curr_int = i
		else:
			if curr_int[0] == i[0] and curr_int[1][1] >= i[1][0]:
				curr_int = (curr_int[0], (curr_int[1][0], i[1][1]))
			else:
				print "%s\t%d\t%d" % (curr_int[0], curr_int[1][0], curr_int[1][1])
				curr_int = i
		
	
	
	print "%s\t%d\t%d" % (curr_int[0], curr_int[1][0], curr_int[1][1])
