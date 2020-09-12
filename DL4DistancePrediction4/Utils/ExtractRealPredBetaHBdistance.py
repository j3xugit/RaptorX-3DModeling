import numpy as np
import cPickle
import os
import sys

if len(sys.argv) < 4:
	print 'Usage: ExtractRealPredBetaHBdistance.py HBBetaFile boundFile distcbFile '
	exit(-1)

betaFile = sys.argv[1]
boundFile = sys.argv[2]
distcbFile = sys.argv[3]

f=open(betaFile, 'rb')
HBBeta = cPickle.load(f)
f.close()

HB = HBBeta['HB']
Beta = HBBeta['BetaPairing']
missing = HBBeta['missing']
existing = [ 1-m for m in missing ]

## mapping from the full sequence (i.e. SEQRES) to DSSP sequence (i.e., SST sequence)
indexMap = np.cumsum(existing) - 1

f = open(boundFile, 'rb')
bounds = cPickle.load(f)
f.close()

predDist = bounds[0]['CbCb'][:,:,0]

distcb = np.loadtxt(distcbFile)

assert distcb.shape == predDist.shape
assert distcb.shape == Beta.shape

print os.path.basename(boundFile).split('.')[0], predDist.shape

for r, c in zip(Beta.row, Beta.col):
	print 'Beta:', r, c, predDist[r, c], distcb[r, c], indexMap[r], indexMap[c]

for r, c in zip(HB.row, HB.col):
	if abs(r-c)<5:
		continue
	print 'HB:', r, c, predDist[r, c], distcb[r, c], indexMap[r], indexMap[c]
