import numpy as np
import sys
import os
import cPickle

from DL4DistancePrediction4.DistanceUtils import CalcDistProb
from DL4DistancePrediction4 import config
from DL4DistancePrediction4.utils import PrettyPrint

if __name__ == "__main__":

	if len(sys.argv)<3:
		print 'Usage: python CalcDistanceDistribution.py proteinListFile folder4DistMatrix'
		print '     This script calculates distance distribution from a list of native distance matrices'
		print '     folder4DistMatrix is the folder containing all the native distance matrix files with a name like proteinName.native.pkl or proteinName.atomDistMatrix.pkl'
		print '	    It is fine if folder4DistMatrix does not contain files for some proteins in proteinListFile'
		print '     The result is written to screen as well as a PKL file'
		exit(1)

	targetFile = sys.argv[1]
	fh = open(targetFile, 'r')
	targets = [ line.strip() for line in list(fh) ]
	fh.close()

	distDir = sys.argv[2]

	allDistMatrices = dict()
	for apt in config.allAtomPairNames:
		allDistMatrices[apt] = []

	numValidTargets = 0
	for target in targets:
		
		distFile = os.path.join(distDir, target + '.native.pkl')
		if not os.path.isfile(distFile):
			distFile = os.path.join(distDir, target + '.atomDistMatrix.pkl')
			if not os.path.isfile(distFile):
				continue

		fh = open(distFile, 'rb')
		matrices = cPickle.load(fh)
		fh.close()
		if distFile.endswith('native.pkl'):
			matrices = matrices['atomDistMatrix']

		for apt in config.allAtomPairNames:
			allDistMatrices[apt].append(matrices[apt])

		numValidTargets += 1

	print 'In total loaded distance matrices for ', numValidTargets, ' targets from file: ', targetFile

	allcounts = dict()
	allcounts['notes'] = 'This is a nested dictionary for distance distribution calculated from proteins in ' + targetFile + '. You may access one entry by [atomPairType][LabelType]. Each entry is a tuple in which the first element is the distribution and the 2nd is the distance boundaries.'
	allcounts['proteins'] = targets

	for apt, matrices in allDistMatrices.iteritems():
		print 'Calculating distance distribution for atom pairs ', apt
		allcounts[apt] = dict()
		for subType in ['47C', '47CPlus', '34C', '34CPlus', '12C', '12CPlus']:
			print 'Distance discretization type: ', subType
			bins = config.distCutoffs[subType]
			count = CalcDistProb(matrices, bins, invalidDistanceSeparated=subType.endswith('Plus'), numRanges=4)
			PrettyPrint( count )
			allcounts[apt][subType] = (count, bins)

	filename = os.path.basename(targetFile).split('.')[0] + '.DistDistribution.pkl'
	fh = open(filename, 'wb')
	cPickle.dump(allcounts, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	fh.close()
			
