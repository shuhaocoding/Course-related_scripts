#!/usr/bin/env python

from __future__ import print_function

import sys
import os

def main(argc, argv):

    celltypes = [ "HEK293T", "HepG2", "K562", "MCF-7" ]
    d120 = os.path.join(os.path.dirname(os.path.realpath(__file__)), "120")
    with open(os.path.join(d120, "tmp0.out"), 'r') as f:
        tmatch = { line.strip(): {} for line in f }
    zeros = {}

    for celltype in celltypes:
        with open(os.path.join(d120, "%s.only.bed.top500.bed.500.txt" % celltype), 'r') as f:
            isheader = True
            for line in f:
                key = line.strip().split('\t')[0]
                tmatch[key][celltype] = ' ' + ' '.join(line.strip().split('\t')[1:])
                if isheader:
                    zeros[celltype] = ' ' + ' '.join([ '0' for i in xrange(len(line.strip().split('\t')[1:])) ])
                    isheader = False

    with open(os.path.join(d120, "tmp4.out"), 'w') as o:
        for key in sorted([ k for k, _ in tmatch.iteritems() ]):
            o.write(key)
            for celltype in celltypes:
                if celltype in tmatch[key]:
                    o.write("\t%s" % tmatch[key][celltype])
                else:
                    o.write("\t%s" % zeros[celltype])
            o.write('\n')

    return 0

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
