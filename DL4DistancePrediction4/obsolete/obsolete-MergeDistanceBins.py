import cPickle
import sys
import os
import numpy as np

import config
import DistanceUtils
import ContactUtils

import getopt

def Usage():

    	print 'python MergeDistanceBins.py PKL_file targetDistLabelType'
	print '  This script merge distance bins in a predicted distance prob matrix.  For example, merge 52 distance bins to 12 or 3 distance bins'
    	print '  PKL_file: a file containing corrected distance prob matrix with name like XXX.correctedDistMatrix.pkl'
    	print '     This file contains a dictionary, in which each key indicates one atom pair type such as CbCb_Discrete52C '
	print '  targetDistLabelType is the target distance bin information, e.g., it can be 12C, 3C, 25C'


def main(argv):


	if len(argv)<2:
		Usage()
		exit(-1)

	predFile = argv[0]
	targetType = argv[1]

	if not os.path.isfile(predFile):
		print 'ERROR: the file for predicted matrix files does not exist: ', predFile
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
