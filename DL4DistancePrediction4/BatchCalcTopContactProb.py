import numpy as np
import os
import sys
import glob
import cPickle

import ContactUtils

def Usage():
	print "python BatchCalcTopContactProb.py folder [fileSuffix]"
	print "	This script calculates the sum of top predicted contact probabilities for all the predicted contact/distance files in a folder"
	print "	folder: a folder containing a few files ending with fileSuffix (default .predictedDistMatrix.pkl) "
	
if len(sys.argv) < 2:
	Usage()
	exit(1)

inDir=sys.argv[1]

fileSuffix = ".predictedDistMatrix.pkl"
if len(sys.argv) >= 3:
	fileSuffix = sys.argv[2]

namePattern = inDir + "/*" + fileSuffix

files = glob.glob(namePattern)

topProbs = []
for f in files:
	with open(f, 'rb') as fh:
		predContactMatrix = cPickle.load(fh)[3]['CbCb']
		topProbSum = ContactUtils.TopContactProbSum(predContactMatrix)
		topProbs.append( (f, topProbSum) )

topProbs = sorted(topProbs, key=lambda x: np.sum(x[1]), reverse=True )

for f, p in topProbs:
	print np.sum(p), f, p

