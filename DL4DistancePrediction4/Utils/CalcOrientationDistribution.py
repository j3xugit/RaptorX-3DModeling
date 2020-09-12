import numpy as np
import sys
import os
import cPickle

from Common import NativeUtils
from DL4DistancePrediction4 import config
from DL4DistancePrediction4.OrientationUtils import CalcOrientationProb
from DL4DistancePrediction4.utils import PrettyPrint

if __name__ == "__main__":

	if len(sys.argv)<3:
		#print 'Usage: python CalcOrientationDistribution.py proteinListFile folder4OrientationMatrix folder4DistMatrix'
		print 'Usage: python CalcOrientationDistribution.py proteinListFile folder4Native'
		print '     This script calculates dihedral/angle distribution from a list of native orientation matrices'
		print '     folder4OrientationMatrix is the folder containing all the native orientation matrix files. Each file has name like proteinName.atomOrientationMatrix.pkl'
		print '     folder4DistMatrix is the folder containing all the native distance matrix files. Each file has name like proteinName.atomDistMatrix.pkl'
		print '	    The results are written to screen as well as a PKL file'
		exit(1)

	targetFile = sys.argv[1]
	fh = open(targetFile, 'r')
	targets = [ line.strip() for line in list(fh) ]
	fh.close()

	oriDir = sys.argv[2]

	allOriMatrices = dict()
	allCbCbMatrices = []

	numValidTargets = 0
	for target in targets:
		f = os.path.join(oriDir, target + '.native.pkl')
		#if not os.path.isfile(f) or not os.path.isfile(df):
		if not os.path.isfile(f):
			continue

		fh = open(f, 'rb')
		data = cPickle.load(fh)
		fh.close()
		matrices = data['atomOrientationMatrix']

		for k, v in matrices.iteritems():
			if allOriMatrices.has_key(k):
				allOriMatrices[k].append(v)
			else:
				allOriMatrices[k] = [v]

		matrices = data['atomDistMatrix']
		if not matrices.has_key('CbCb'):
			print 'ERROR: no Cb-Cb distance matrix found in the file: ', df
			exit(1)
		allCbCbMatrices.append(matrices['CbCb'])

		numValidTargets += 1

	print 'In total loaded orientation matrices for ', numValidTargets, ' targets from file: ', targetFile

	allcounts = dict()
	allcounts['notes'] = 'This is a nested dictionary for orientation distribution calculated from proteins in ' + targetFile + '. You may access one entry by [oirentationName][subType]. Each entry is a tuple in which the first element is the distribution and the 2nd is the discrete boundaries.'
	allcounts['proteins'] = targets

	for labelName, matrices in allOriMatrices.iteritems():
		print 'Calculating distribution for orientation ', labelName
		allcounts[labelName] = dict()
		if labelName in config.allDihedralNames:
			for subType in ['37C', '37CPlus', '25C', '25CPlus']:
				print 'dihedral discretization type: ', subType
				bins = config.dihedralCutoffs[subType]
				count = CalcOrientationProb(matrices, bins, distMatrices=allCbCbMatrices, invalidEntrySeparated=subType.endswith('Plus'), numRanges=4)
				PrettyPrint( count )
				allcounts[labelName][subType] = (count, bins)
		elif labelName in config.allAngleNames:
			for subType in ['19C', '19CPlus', '13C', '13CPlus']:
				print 'angle discretization type: ', subType
				bins = config.angleCutoffs[subType]
				count = CalcOrientationProb(matrices, bins, distMatrices=allCbCbMatrices, invalidEntrySeparated=subType.endswith('Plus'), numRanges=4)
				PrettyPrint( count )
				allcounts[labelName][subType] = (count, bins)

	filename = os.path.basename(targetFile).split('.')[0] + '.OriDistribution.pkl'
	fh = open(filename, 'wb')
	cPickle.dump(allcounts, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	fh.close()
			
