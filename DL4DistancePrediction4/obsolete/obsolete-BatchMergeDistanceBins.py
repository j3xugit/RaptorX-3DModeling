import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np

import config
import DistanceUtils
import ContactUtils
##from MergePredictedContactMatrix import MergeAndSaveOneProtein

import getopt

def Usage():

    	print 'python BatchMergeDistanceBins.py poteinList PKL_folder targetDistLabelType'
	print '  This script merges distance bins in a predicted distance prob matrix for a list of proteins'
	print '  For example, it can be used to merge 52 distance bins to 12 distance bins or even to 3 distance bins'
    	print '  PKL_folder: a folder containing distance prob matrix files with name XXX.correctedDistMatrix.pkl '
    	print '     A distance prob matrix file contains a dictionary, in which each key indicates one atom pair type such as CbCb_Discrete52C '
	print '  targetDistLabelType is the target distance bin information, e.g., it can be 12C, 3C, 25C'


def main(argv):


	if len(argv)<3:
		Usage()
		exit(1)

	proteinListFile = argv[0]
	predFolder = argv[1]
	targetType = argv[2]

	fileSuffix = '.correctedDistMatrix.pkl'

	if not os.path.isfile(proteinListFile):
		print 'the protein list file does not exist: ', proteinListFile
		exit(1)

	if not os.path.isdir(predFolder):
		print 'the folder for predicted matrix files does not exist: ', predFolder
		exit(1)

	fh = open(proteinListFile, 'r')
	proteins = [ line.strip() for line in list(fh) ]
	fh.close()

	for protein in proteins:
		predFile = os.path.join( predFolder, protein + fileSuffix ) 
		if not os.path.isfile(predFile):
			print 'the prediction file does not exist: ', predFile
			exit(1)

		fh = open(predFile, 'rb')
		pred = cPickle.load(fh)
		fh.close()

		newMatrices = dict()
		for response, m in pred.iteritems():
			apt = config.Response2LabelName(response)
			labelType = config.Response2LabelType(response)
			subType = labelType[len('Discrete'): ]

			srcDistCutoff = config.distCutoffs[subType]
			dstDistCutoff = config.distCutoffs[targetType]
			newDistMatrix = DistanceUtils.MergeDistanceBins(m, srcDistCutoff, dstDistCutoff)		

			newMatrices[apt] = newDistMatrix

		fields = predFile.split('.')
        	savefile = '.'.join( [fields[0], '.'.join(fields[1:-1]) + '4' + targetType, fields[-1] ] )
		fh = open(savefile, 'wb')
		cPickle.dump(newMatrices, fh, protocol=cPickle.HIGHEST_PROTOCOL)
		fh.close()



if __name__ == "__main__":
    	main(sys.argv[1:])
