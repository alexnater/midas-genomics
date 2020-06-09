#!/usr/bin/env python3

"""
Takes a glob pattern of SFS files, sample sizes, output sfs name, and output format string and combines them into a single SFS.

"""

import sys
import os
import glob
import numpy as np

glob_pattern = sys.argv[1]
samplesizes = [ int(x) for x in sys.argv[2].split(',') ]
outfile = sys.argv[3]
outformat = sys.argv[4]


infiles = glob.glob(glob_pattern)
dimensions = [ 2 * x + 1 for x in samplesizes ]
sfs = np.zeros(np.prod(dimensions), dtype=float)

for infile in infiles:
    print("Reading in file {} ...".format(os.path.basename(infile)), file=sys.stderr)
    temp = np.fromfile(infile, dtype=float, count=-1, sep=" ")
    if not temp.size == np.prod(dimensions):
        raise Exception("ERROR: Incorrectly formated sfs file {}!".format(infile))
    sfs += temp

mask = np.zeros(sfs.size, dtype=int)
mask[0] = 1
mask[-1] = 1
sfs[0] += sfs[-1]
sfs[-1] = 0

print("Printing file {} ...".format(outfile), file=sys.stderr)
with open(outfile, 'w') as outhandle:
    if outformat == 'dadi':
        print("# {}D-SFS from ANGSD, total length of sequence: L={}".format(len(dimensions), sfs.sum()), file=outhandle)
        print("{}\tunfolded".format("\t".join(map(str, dimensions))), file=outhandle)
        sfs.tofile(outhandle, sep=" ", format="%.12f")
        print("\n", end="", file=outhandle)
        mask.tofile(outhandle, sep=" ", format="%d")
    elif outformat == 'fastsimcoal':
        print("1 observations. No. of demes and sample sizes are on next line", file=outhandle)
        print("{} {}".format(len(dimensions), " ".join((str(x - 1) for x in dimensions))), file=outhandle)
        sfs.tofile(outhandle, sep=" ", format="%.12f")
    else:
        raise Exception("ERROR: Output format {} not recognized!".format(outformat))



