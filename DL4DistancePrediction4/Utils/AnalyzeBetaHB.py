## this script calculates the distribution of HB distance and Beta distance

import sys
import os
import numpy as np
import cPickle

if len(sys.argv) < 2:
	print 'Usage: AnalyzeBetaHBdistance.py *.HBBeta.pkl'
	exit(-1)

## load parsedSST file

parsedSSTfile = sys.argv[1]
#distcbFile = sys.argv[2]

f = open(parsedSSTfile, 'rb')
sst=cPickle.load(f)
f.close()

"""
distcb = np.loadtxt(distcbFile)
assert (distcb.shape[0] == len(sst['sequence']) )
"""

Beta = sst['BetaPairing']
LRcount = 0
MRcount = 0
SRcount = 0
NRcount = 0
count = len(Beta.row)
for r, c in zip(Beta.row, Beta.col):
	if abs(r-c) >= 24:
		LRcount += 1
	elif abs(r-c) >=12:
		MRcount += 1
	elif abs(r-c) >=6:
		SRcount += 1
	elif r!=c:
		NRcount += 1

HB = sst['HB']
HBLRcount = 0
HBMRcount = 0
HBSRcount = 0
HBNRcount = 0
HBcount = len(HB.row)

for r, c in zip(HB.row, HB.col):
	if abs(r-c)>=24:
		HBLRcount += 1
	elif abs(r-c)>=12:
		HBMRcount += 1
	elif abs(r-c)>=6:
		HBSRcount += 1
	elif r!=c:
		HBNRcount += 1

seqLen = Beta.shape[0]
assert ( Beta.shape == HB.shape)

print os.path.basename(parsedSSTfile).split('.')[0], seqLen, count, NRcount, SRcount, MRcount, LRcount, HBcount, HBNRcount, HBSRcount, HBMRcount, HBLRcount, 12*seqLen-30, 12*seqLen-108, 24*seqLen-432, (seqLen-23)*(seqLen-23)
