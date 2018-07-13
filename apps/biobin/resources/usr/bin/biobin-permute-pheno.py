#!/usr/bin/env python

import random
import sys


if len(sys.argv) < 2:
    print "usage: %s <filename.phe> [seed]" % (sys.argv[0],)
    sys.exit(2)

iids = list()
phenos = list()
with open(sys.argv[1],'rU') as infile:
    header = infile.next().strip()
    for line in infile:
        iid,pheno = line.split(None,1)
        iids.append(iid)
        phenos.append(pheno.strip())

if len(sys.argv) > 2:
    random.seed(long(sys.argv[2]))
random.shuffle(iids)

print header
for i in xrange(len(iids)):
    print iids[i] + "\t" + phenos[i]
