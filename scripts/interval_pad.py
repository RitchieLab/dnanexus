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
	
	n_splits=1
	if len(sys.argv) > 2:
		try:
			n_splits = int(sys.argv[2])
		except ValueError:
			print >> sys.stderr, "WARNING: argument", sys.argv[2], "not convertible to int, setting number of splits to 1."
	
	if n_splits < 1:
		print >> sys.stderr, "WARNING: Requested < 1 number of splits, setting number of splits to 1."
	
	
	intervals = []
	orig_intervals=[]
	for l in sys.stdin:
		items = l.strip().split()
		try:
			intervals.append((items[0], (int(items[1]) - interval_pad, int(items[2]) + interval_pad)))
			orig_intervals.append((items[0], (int(items[1]), int(items[2]))))
		except ValueError:
			print >> sys.stderr, "Warning: unable to parse line '%s', ignoring" % l.strip()
	
	intervals.sort()
	
	
	curr_int = None
	final_intervals = []
	final_orig = []
	curr_orig=[]
	for idx, i in enumerate(intervals):
		if curr_int is None:
			curr_int = i
			curr_orig=[orig_intervals[idx]]
		else:
			if curr_int[0] == i[0] and curr_int[1][1] >= i[1][0]:
				curr_int = (curr_int[0], (curr_int[1][0], i[1][1]))
				curr_orig.append(orig_intervals[idx])
			else:
				final_intervals.append((curr_int[0], curr_int[1][0], curr_int[1][1]))
				final_orig.append(curr_orig)
				curr_int = i
				curr_orig=[orig_intervals[idx]]
	
	
	final_intervals.append((curr_int[0], curr_int[1][0], curr_int[1][1]))
	final_orig.append(curr_orig)
	# OK, now, if we request the number of splits > 1
	if n_splits > 1:
		cs=1
		si=1
		n_int = len(final_intervals) / n_splits + (len(final_intervals) % n_splits > 0)
		for idx in range(len(final_intervals)):
			for ci in final_orig[idx]:
				print "%d:%s\t%d\t%d" % (cs, ci[0], ci[1][0], ci[1][1])

			si+=1
			if si > n_int:
				si=1
				cs+=1		
	else:
		for ci in final_intervals:
			print "%s\t%d\t%d" % ci
