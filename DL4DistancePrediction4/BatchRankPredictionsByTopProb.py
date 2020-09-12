import numpy as np
import os
import sys
import cPickle

import ContactUtils

##this script ranks all predictions of the same protein by the sum of top predicted contact probabilities 

def Usage():
	print "python BatchRankPredictionsByTopProb.py proteinList ResDir predFolder1 predFolder2 predFolder3 ..."
	print "		this script ranks the predictedDistMatrix of all proteins in the proteinList by top contact probability"
	print "		each folder contains predictions of all proteins generated under one setting"
	print "		each predFile ends with .predictedDistMatrix.pkl"
	
if len(sys.argv) < 4:
	Usage()
	exit(1)

proteinListFile = sys.argv[1]
resDir = sys.argv[2]
if not os.path.isdir(resDir):
	os.mkdir(resDir)

predFolders=sys.argv[3:]

with open(proteinListFile, 'r') as fh:
	proteins = [ line.strip() for line in list(fh) ]

for protein in proteins:
	topProbs = []
	for predFolder in predFolders:
		predFile = os.path.join(predFolder, protein + '.predictedDistMatrix.pkl')
		with open(predFile, 'rb') as fh:
			predContactMatrix = cPickle.load(fh)[3]['CbCb']
		topProbSum = ContactUtils.TopContactProbSum(predContactMatrix)
		topProbs.append( (predFile, topProbSum) )

	topProbs = sorted(topProbs, key=lambda x: np.sum(x[1]), reverse=True )
	topFile = topProbs[0][0]
	print topFile
	bname = os.path.basename(topFile)
	targetFile = os.path.join(resDir, bname)
	os.symlink(os.path.abspath(topFile), targetFile)

