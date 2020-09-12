import cPickle
import sys
import os
import numpy as np

import getopt

import config
import DistanceUtils

def Usage():
    	print 'python ReducePredictedDistMatrix.py [-r targetResponse | -s savefolder] predictedDistMatrix_PKL_file '
	print '  This script converts a predicted distance prob matrix from a finer-grained discretization scheme to a coarser-grained one'
	print '       e.g., from 52 distance bins to 12 distance bins'
	print '	 -r: the resultant distance matrix is only for the target response, e.g., CbCb_Discrete14C (default), CbCb+CaCa_Discrete14C, AllAP_Discrete14C'
	print '		currently Plus or Minus is not supported, e.g., 12CPlus, 3CMinus '
	print '	 -s: the folder for result save, default current work directory'
    	print '  PKL_file: a distance prob matrix file named after XXX.predictedDistMatrix.pkl where XXX is the protein name'
	print '  the result is saved into a file XXX.predictedDistMatrix4YYY.pkl where YYY is the targetSubType, e.g., 14C'


def main(argv):

	targetResponse = 'CbCb_Discrete14C'
	savefolder = os.getcwd()

	if len(argv)<1:
		Usage()
		exit(1)

	try:
                opts, args = getopt.getopt(argv,"r:s:",["response=", "savefolder="])
                #print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

        if len(args) != 1:
                Usage()
                exit(1)

	for opt, arg in opts:
                if opt in ("-r", "--response"):
                        targetResponse = arg
		elif opt in ("-s", "--savefolder"):
			savefolder = arg
			if not os.path.isdir(savefolder):
				os.mkdir(savefolder)
		else:
			Usage()
			exit(1)
	
	predFile = args[0]
	if not os.path.isfile(predFile):
		print 'ERROR: the predicted distance prob matrix file does not exist: ', predFile
		exit(1)

	targetLabelName, targetLabelType, targetSubType = config.ParseResponse(targetResponse)

	if targetSubType.endswith('Plus') or targetSubType.endswith('Minus'):
		print 'ERROR: currently the target distance type cannot end with Plus or Minus'
		exit(1)

	if not config.distCutoffs.has_key(targetSubType):
		print 'ERROR: the key ', targetSubType, ' is not defined in config.distCutoffs'
		exit(1)

	dstDistCutoff = config.distCutoffs[targetSubType]

	with open(predFile, 'rb') as fh:
		pred = cPickle.load(fh)
	target, sequence, distProbMatrix, contProbMatrix, labelWeight, labelDistribution = pred[:6]

	newDistProbMatrix = dict()
	newLabelWeight = dict()
	newLabelDistribution = dict()

	for response in distProbMatrix.keys():
		labelName, labelType, subType = config.ParseResponse(response)
		if labelName not in config.allAtomPairNames:
			continue

		if labelName != targetLabelName:
			continue

		if not config.distCutoffs.has_key(subType):
			print 'ERROR: the dist prob matrix to be reduced has a discretization scheme undefined in config.distCutoffs: ', subType
			exit(1)

		srcDistCutoff = config.distCutoffs[subType]

		## convert distProbMatrix
		newDistProbMatrix[targetResponse] = DistanceUtils.MergeDistanceBinsBySum(distProbMatrix[response], srcDistCutoff, dstDistCutoff)		

		## convert labelDistribution
		newLabelDistribution[targetResponse] = DistanceUtils.MergeDistanceBinsBySum(labelDistribution[response], srcDistCutoff, dstDistCutoff)

		## convert labelWeight
		newLabelWeight[targetResponse] = DistanceUtils.MergeDistanceBinsByAverage(labelWeight[response], srcDistCutoff, dstDistCutoff, labelDistribution[response] )


	fields = os.path.basename(predFile).split('.')
        savefile = '.'.join( [fields[0], '.'.join(fields[1:-1]) + '4' + targetSubType, fields[-1] ] )
	savefile = os.path.join(savefolder, savefile)
	with open(savefile, 'wb') as fh:
		cPickle.dump((target, sequence, newDistProbMatrix, contProbMatrix, newLabelWeight, newLabelDistribution), fh, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    	main(sys.argv[1:])
