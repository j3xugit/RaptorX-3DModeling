import numpy as np
import sys
import os

if len(sys.argv) < 1:
	print 'Usage: AnalyzeDeepScore.py inputFile'
	exit(1)

dsFile = sys.argv[1]

name = os.path.basename(dsFile).split('.')[0]

dsfh = open(dsFile, 'r')
dsContent = [ line.strip() for line in list(dsfh) ]
dsfh.close()

if len(dsContent) < 1:
	print "Empty file: ", dsFile
	exit(1)

if dsContent[0].find(name) == -1:
	print 'the deepscore file has content inconsistent with target name: ', name
	exit(1)

def ParseDeepScoreFile(dsContent):

	method1_quality = []
	method1_first = []

	for line in dsContent:
		if line.find('model'):
			method1_quality.append( [ np.float32(v) for v in line.split()[3:8] ] )

		if line.find('model1.pdb') >=0:
			method1_first = [ np.float32(v) for v in line.split()[3:8] ]

	method1_quality.sort(key=lambda x: (x[1]+x[3])/2 )

	return method1_first, method1_quality[-1]


def ParseMaxClusterFile(msContent, dsContent):
	method1 = 'DistEC-t2-c15'
        method1_flag = 'Sig1.0'

	method2 = 'ContactEC'
        method2_flag = 'UP3.0'

	rankedModels = []
	numModels = 0
	for line in msContent[9:]:
		if line.find(method1) >=0 and line.find(method1_flag) >=0:
			rankedModels.append(line.split()[6])
			numModels += 1

		if line.find(method2) >=0 and line.find(method2_flag) >=0:
			rankedModels.append(line.split()[6])
			numModels += 1

		if numModels >=5:
			break

	##print rankedModels

	rankedQuality = []

	for model in rankedModels:
		for line in dsContent:
			if line.find(model) >=0 :
				quality = [ np.float32(v) for v in line.split()[4:8] ]
				rankedQuality.append(quality)

	quality_first = rankedQuality[0]

	rankedQuality.sort(key=lambda x: (x[0]+x[2])/2 )
	quality_best = rankedQuality[-1]

	return quality_first, quality_best
				

methodResults = ParseDeepScoreFile(dsContent)

finalResults = methodResults[0] + methodResults[1]

print '\t'.join([name] + [ str(v) for v in finalResults])
