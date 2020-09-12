import cPickle
import sys
import os
import numpy as np

import config
from config import ParseResponse

import LabelUtils
import DistanceUtils
import ContactUtils
import OrientationUtils

import getopt

def Usage():
    	print 'python AddRVMatrix2PredictedMatrix.py [-w weight | -s savefile ] predictedMatrixFile RVMatrixFile'
	print '	This script adds a real-valued (native) matrix to a predicted dist/orientation matrix'
	print '	-w: weight for the real-valued matrix'
	print ' -s: file for result saving'
    	print '	predictedMatrixFile: a PKL file for predicted matrix with name like XXX.predictedDistMatrix.pkl'
    	print '	    This file is a tuple of 6 items: name, primary sequence,  predicted distance prob, predicted contact prob, labelWeights and labelDistribution'
	print '	RVMatrixFile: a PKL file for real-valued matrix, e.g., XXX.native.pkl'

def AddRVMatrix(predFile, RVFile, w4RV=1.):
	assert predFile is not None
	assert RVFile is not None

        content = DistanceUtils.LoadRawDistProbFile(predFile)
        name, sequence, predProbMatrix, predContactMatrix, labelWeight, labelDistribution = content[:6]

	with open(RVFile, 'r') as fh:
		RV = cPickle.load(fh)

	RVseq = RV['sequence']
	assert sequence == RVseq

	RVMatrix = RV['atomDistMatrix']
	RVMatrix.update(RV['atomOrientationMatrix'])

	resProbMatrix = {}
        for response in predProbMatrix.keys():
		labelName, labelType, subType = ParseResponse(response)
		if not RVMatrix.has_key(labelName): 
			resProbMatrix[response] = predProbMatrix[response]
			continue

		if labelName in config.allDistLabelNames:
			distm = RVMatrix[labelName]
			labelMatrix, _, _  = DistanceUtils.DiscretizeDistMatrix(distm, config.distCutoffs[subType], invalidDistanceSeparated=(subType.endswith('Plus') or subType.endswith('Minus')) )

		elif labelName in config.allOrientationNames:
			orim = RVMatrix[labelName]
			if labelName in config.allDihedralNames:
                        	bins = config.dihedralCutoffs[subType]
                        else:
                                bins = config.angleCutoffs[subType]

			labelMatrix, _, _ = OrientationUtils.DiscretizeOrientationMatrix(orim, bins=bins, distThreshold4Orientation=20, distMatrix=RVMatrix['CbCb'], invalidEntrySeparated=(subType.endswith('Plus') or subType.endswith('Minus')) )
		else:
			print 'ERROR: unsupported label name: ', labelName
			exit(1)
		oneHotMatrix = LabelUtils.GetOneHotLabelMatrix(labelMatrix, numLabels=config.GetResponseProbDims(response) )
		resProbMatrix[response ] = (predProbMatrix[response] + w4RV * oneHotMatrix )/(1+w4RV)

	## here handle the situation when different label types of the same labelName appear in predFile
	contactMatrixProb = dict()
	contactMatrixCount = dict()
	for response in resProbMatrix.keys():
		labelName, labelType, subType = ParseResponse(response)
		if labelName not in config.allAtomPairNames:
			continue

		labelOf8 = DistanceUtils.LabelsOfOneDistance(config.ContactDefinition, config.distCutoffs[subType])
		if not contactMatrixProb.has_key(labelName):
			contactMatrixProb[labelName] = ContactUtils.Distance2Contact(resProbMatrix[response], labelOf8)
			contactMatrixCount[labelName] = 1
		else:
			contactMatrixProb[labelName] += ContactUtils.Distance2Contact(resProbMatrix[response], labelOf8)
			contactMatrixCount[labelName] += 1

	for labelName in contactMatrixProb.keys():
		contactMatrixProb[labelName] = contactMatrixProb[labelName]/contactMatrixCount[labelName]

        return (name, sequence, resProbMatrix, contactMatrixProb, labelWeight, labelDistribution)

def main(argv):
	w4RV = 1.
	savefile = None

    	try:
        	opts, args = getopt.getopt(argv,"w:s:",["weight=", "savefile="])
        	#print opts, args
    	except getopt.GetoptError:
        	Usage()
        	exit(1)

	if len(args) != 2:
		Usage()
		exit(1)

	predFile, RVFile = args

    	for opt, arg in opts:
		if opt in ("-w", "--weight"):
			w4RV = np.float32(arg)
			assert w4RV > 0
		elif opt in ("-s", "--savefile"):
			savefile = arg
		else:
	    		Usage()
	    		exit(1)

	content4save = AddRVMatrix(predFile, RVFile, w4RV)
	assert len(content4save) >= 6

	if savefile is None:
		savefile = content4save[0] + '-RV.predictedDistMatrix.pkl'

        with open(savefile, 'wb') as fh:
        	cPickle.dump(content4save, fh, protocol = cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    	main(sys.argv[1:])
