## this script calculates the distribution of HB distance and Beta distance

import sys
import os
import numpy as np
import cPickle

if len(sys.argv) < 3:
	print 'Usage: AnalyzeBetaHBdistance.py *.HBBeta.pkl *.distcb'
	exit(-1)

## load parsedSST file

parsedSSTfile = sys.argv[1]
distcbFile = sys.argv[2]

f = open(parsedSSTfile, 'rb')
sst=cPickle.load(f)
f.close()

distcb = np.loadtxt(distcbFile)
assert (distcb.shape[0] == len(sst['sequence']) )

"""
HB = sst['HB']
HBrow = HB.row
HBcol = HB.col
for r, c in zip(HBrow, HBcol):
	print distcb[r, c]
"""

Beta = sst['BetaPairing']
for r, c in zip(Beta.row, Beta.col):
	print distcb[r, c]

