import numpy as np
import sys
import os
import cPickle

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import config
from utils import PrettyPrint


def CalcDistAverageOfSingleMatrix(distM, minDist=None, maxDist=None, maxSeqSep=30):
	m=np.ones_like(distM, dtype=np.int32)

	if minDist is not None:
		np.putmask(m, distM<minDist, 0)

	if maxDist is not None:
		np.putmask(m, distM>maxDist, 0)

	np.putmask(m, distM<0.001, 0)

	nzIndices = m.nonzero()

	sumDist = np.zeros(maxSeqSep + 1, dtype=np.float32)
	count = np.zeros(maxSeqSep + 1, dtype= np.int32)
	
	for i, j in zip(nzIndices[0], nzIndices[1]):
		offset = min(abs(i-j), maxSeqSep)
		sumDist[offset] += distM[i, j]
		count[offset] +=1

	avg = [ s/(0.00001+c) for s, c in zip(sumDist, count) ]

	return avg, count
	
	

if __name__ == "__main__":

	maxSeqSep = 40
	minDist = None
	maxDist = None

	if len(sys.argv)<4:
		print 'Usage: python CalcMeanDistance.py proteinListFile folder4DistMatrix min-max'
		print '     This script calculates the average distance in a distance interval defined by [min, max] from a list of native distance matrices'
		print ' 		the average distance for each different sequence separation will be calculated'
		print '     folder4DistMatrix is the folder containing all the native distance matrix files. Each file has name like proteinName.atomDistMatrix.pkl'
		print '	The results will be written to screen as well as a PKL file'
		exit(1)

	targetFile = sys.argv[1]
	fh = open(targetFile, 'r')
	targets = [ line.strip() for line in list(fh) ]
	fh.close()

	distDir = sys.argv[2]

	distBoundStr= sys.argv[3]
	fields = distBoundStr.split('-')
	print fields
	if len(fields) == 2:
		if fields[0] !='' :
			minDist = np.float32(fields[0])

		if fields[1] !='' :
			maxDist = np.float32(fields[1])
	else:
		print 'ERROR: incorrect format for distance bound. Allowed examples are 5.0-10.0, -4.0, 16.0-'
		exit(1)

	print minDist, maxDist

	allDistMatrices = dict()
	for apt in config.allAtomPairTypes:
		allDistMatrices[apt] = []

	numValidTargets = 0
	for target in targets:
		distFile = os.path.join(distDir, target + '.atomDistMatrix.pkl')
		if not os.path.isfile(distFile):
			continue

		fh = open(distFile, 'rb')
		matrices = cPickle.load(fh)
		fh.close()

		for apt in config.allAtomPairTypes:
			allDistMatrices[apt].append(matrices[apt])

		numValidTargets += 1

	print 'In total loaded distance matrices for ', numValidTargets, ' targets from file: ', targetFile

	allcounts = dict()
	allcounts['notes'] = 'This is a dictionary for distance average calculated from proteins in ' + targetFile + '. You may access one entry by [atomPairType][1-30]. Each entry is a tuple in which the first element is the distribution and the 2nd is the distance boundaries.'
	allcounts['proteins'] = targets

	for apt, mlist in allDistMatrices.iteritems():
		print 'Calculating distance distribution for atom pairs ', apt

		accAvg = np.zeros(maxSeqSep + 1, dtype=np.float32)
		accCount = np.zeros(maxSeqSep+1, dtype=np.int32)
		for m in mlist:
			avg, count = CalcDistAverageOfSingleMatrix(m, minDist=minDist, maxDist=maxDist, maxSeqSep=maxSeqSep)
			accAvg += avg * count
			accCount += count
		accAvg = accAvg /(0.00001+accCount)

		allcounts[apt] = (accAvg, accCount)

		print apt
		print allcounts[apt]
		

	filename = os.path.basename(targetFile).split('.')[0] + '.MeanDistance.pkl'
	fh = open(filename, 'wb')
	cPickle.dump(allcounts, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	fh.close()

