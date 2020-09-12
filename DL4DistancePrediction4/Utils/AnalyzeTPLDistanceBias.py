import os
import sys
import numpy as np
import cPickle
from config import distCutoffs
from DistanceUtils import DiscretizeDistMatrix

def BinCount(queryDistMatrix, tplDistMatrix, labelType='25C'):
	if labelType == '25C':
		numLabels = 25
	elif labelType == '52C':
		numLabels = 52
	else:
		print 'ERROR: unsupported labelType in BinCount: ', labelType
		exit(-1)

	counts = dict()
	for k in queryDistMatrix.keys():
		if k not in set(['CbCb', 'CaCa', 'CgCg', 'CaCg', 'NO']):
			continue

		qMatrix = queryDistMatrix[k]
		tMatrix = tplDistMatrix[k]

		qLabelMatrix, _, _ = DiscretizeDistMatrix(qMatrix, bins=distCutoffs[labelType])
		tLabelMatrix, _, _ = DiscretizeDistMatrix(tMatrix, bins=distCutoffs[labelType])
		qtLabelMatrix = tLabelMatrix * numLabels + qLabelMatrix

		count = np.bincount(qtLabelMatrix.flatten(), minlength = numLabels * numLabels )
		counts[k] = count.reshape( (numLabels, numLabels) )

	return counts

def LoadTPLDistanceFeatureFiles(files=None):
	if files is None or len(files)==0:
                print 'Please provide valid TPL-based feature files'
                exit(-1)

        fhs = [ open(file, 'rb') for file in files ]
        data = sum([ cPickle.load(fh) for fh in fhs ], [])
        [ fh.close() for fh in fhs ]

        distances = []

	for d, counter in zip(data, xrange(len(data)) ):
		## read in the native atom dist matrix
		if not d.has_key('atomDistMatrix'):
			print 'the TPL-based feature file for pair ', d['name'], 'does not have native distance information for the query protein'
			continue

		if not d.has_key('tplDistMatrix'):
			print 'the TPL-based feature file for pair ', d['name'], 'does not have distance information for the template protein'
			continue

		distances.append( ( d['atomDistMatrix'], d['tplDistMatrix'] ) )

		if counter> 0 and counter%500 ==0 :
			print 'loaded information for %d protein pairs' % counter

	print 'In total loaded distance information for %d protein pairs ' % (len(distances) )

	return distances

def main(argv):
	files = argv
	distances = LoadTPLDistanceFeatureFiles(files)

	counts = dict()

	for d in distances:
		count = BinCount( d[0], d[1])
		for k, v in count.iteritems():
			if not counts.has_key(k):
				counts[k] = v
			else:
				counts[k] += v

	## normalize counts
	newCounts = dict()
	refProbs = dict()
	for k, v in counts.iteritems():
		rowSum = np.sum(v, axis=1, keepdims=True)
		newV = v*1. / rowSum
		newCounts[k] = newV

		oneSum = np.sum(v, axis=0)
		refProb = oneSum*1. / np.sum(oneSum)
		refProbs[k] =refProb

	print newCounts
	print refProbs



if __name__ == '__main__':
        main(sys.argv[1:])



	
