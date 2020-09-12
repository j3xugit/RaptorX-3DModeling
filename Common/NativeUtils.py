import os
import sys
import cPickle
import numpy as np

def LoadGroundTruth(name, location):
	filename = os.path.join(location, name + '.native.pkl')
	if not os.path.isfile(filename):
		print 'ERROR: ground truth file does not exist: ', filename
		exit(1)

	with open(filename, 'r') as fh:
		truth = cPickle.load(fh)

	return truth

def AddGroundTruth(rawData, truth):
	if rawData['sequence'] != truth['sequence']:
		print 'ERROR: mismatch between a query sequence and the sequence in its ground truth file: ', rawData['name']
		exit(1)

	seqLen = rawData['length']
	if truth.has_key('atomDistMatrix'):
		atomDistMatrix = truth['atomDistMatrix']
		if seqLen != atomDistMatrix['CbCb'].shape[0]:
			print 'ERROR: inconsistent query seq length with its dist matrix for ', rawData['name']
                	exit(1)
		rawData['atomDistMatrix'] = atomDistMatrix

	if truth.has_key('atomOrientationMatrix'):
		rawData['atomOrientationMatrix'] = truth['atomOrientationMatrix']

