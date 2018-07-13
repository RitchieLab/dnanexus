#!/usr/bin/env python

import collections
import sys


if (len(sys.argv) < 3) or ("%d" not in sys.argv[2]):
    print "usage: %s <baseline> <pathformat> [permutations]" % (sys.argv[0],)
    print "\t<baseline> is the unpermuted summary.tsv file"
    print "\t<pathformat> is the path to the permuted summary.tsv files,"
    print "\t\twhich must contain \"%%d\" in place of the permutation number"
    print "\t[permutations] is the number of permutations to scan; default=1000"
    sys.exit(2)


pMin = 1
pMax = 1000
if len(sys.argv) > 3:
    pRange = sys.argv[3].split('-',1)
    if len(pRange) > 1:
        pMin = int(pRange[0])
        pMax = int(pRange[1])
    else:
        pMax = int(pRange[0])

results = list()
resultTestPval = collections.defaultdict(dict)
resultTestBetter = collections.defaultdict(dict)
with open(sys.argv[1],'rU') as infile:
    header = infile.next().strip()
    tests = header.split("\t")[11:]
    l = 1
    for line in infile:
        l += 1
        cols = line.strip().split("\t")
        if len(cols) != len(tests) + 11:
            sys.exit("ERROR: found %d columns on line %d of %s, expected %d" % (len(cols),l,sys.argv[1],len(tests)+11))
        result = (cols[0],cols[1])
        results.append(result)
        for c in xrange(11,len(cols)):
            test = tests[c - 11]
            resultTestPval[result][test] = float(cols[c])
            resultTestBetter[result][test] = 0

for p in xrange(pMin,pMax+1):
    filename = sys.argv[2] % (p,)
    with open(filename,'rU') as infile:
        header = infile.next().strip()
        l = 1
        for line in infile:
            l += 1
            cols = line.strip().split("\t")
            if len(cols) != len(tests) + 11:
                sys.exit("ERROR: found %d columns on line %d of %s, expected %d" % (len(cols),l,filename,len(tests)+11))
            result = (cols[0],cols[1])
            for c in xrange(11,len(cols)):
                test = tests[c - 11]
                if float(cols[c]) < resultTestPval[result][test]:
                    resultTestBetter[result][test] += 1

header = "Outcome\tBin"
for test in tests:
    header += "\t" + test
print header
for result in results:
    line = "%s\t%s" % result
    for test in tests:
        b = resultTestBetter[result][test]
        line += "\t%s%g" % (("" if (b > 0) else "<"), (b if (b > 0) else 1) / float(pMax - pMin + 1))
    print line
