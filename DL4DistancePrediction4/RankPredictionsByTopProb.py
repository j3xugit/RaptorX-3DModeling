import numpy as np
import os
import sys
import cPickle

#import config
from ContactUtils import TopContactProbSum
#import OrientationUtils

def Usage():
	print "python RankPredictionsByTopProb.py predFile1 predFile2 predFile3 ..."
	print "	This script ranks the predicted dist/orientation matrices of one protein by top predicted contact probabilities"
	print "	predFile: one predicted dist/orientation matrix file, ending with .predictedDistMatrix.pkl"
	print "	the result is a ranking list of all predictions and their associated prob sums"
	
if len(sys.argv) < 2:
	Usage()
	exit(1)

predFiles=sys.argv[1:]

topProbs = []
for f in predFiles:
	with open(f, 'rb') as fh:
		pred = cPickle.load(fh)

	CbCbProbSum = None

	distProbSumPool = []
	for apt in ['CbCb', 'CaCa', 'NO']:
		if not pred[3].has_key(apt):
			continue
		predContactMatrix = pred[3][apt]
		topProbSum = TopContactProbSum(predContactMatrix)
		distProbSumPool.append(topProbSum)
		if apt == 'CbCb':
			CbCbProbSum = np.average(topProbSum)

	distProbSum = np.average(distProbSumPool)

	topProbs.append( (f, distProbSum, CbCbProbSum) )

	"""
	oriContactMatrices = dict()
	for response, predMatrix in pred[2].iteritems():
		labelName, labelType, subType = config.ParseResponse(response)
		if labelName not in config.allOrientationNames:
			continue
		oriContactMatrices[labelName] = OrientationUtils.DeriveOriContactMatrix(predMatrix, response)

	oriProbSumPool = []
	for ori in ['Ca1Cb1Cb2Ca2','N1Ca1Cb1Cb2','Ca1Cb1Cb2']:
		if not oriContactMatrices.has_key(ori):
			continue
		topProbSum = ContactUtils.TopContactProbSum(oriContactMatrices[ori])
		oriProbSumPool.append(topProbSum)

	oriProbSum = np.average(oriProbSumPool)

	probSum = (distProbSum + oriProbSum)/2

	topProbs.append( (f, probSum, distProbSum, oriProbSum, CbCbProbSum) )
	"""

topProbs = sorted(topProbs, key=lambda x: x[1], reverse=True )

for p in topProbs:
	f, distProb, CbCbProb = p
	print f, distProb, CbCbProb

