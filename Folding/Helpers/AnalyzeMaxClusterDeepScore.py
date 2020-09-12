import numpy as np
import sys
import os

if len(sys.argv)<3:
	print 'Usage: python AnalyzeMaxClusterDeepScore.py MaxClusterResultFile DeepScoreResultFile'
	exit(-1)

mcFile = sys.argv[1]
dsFile = sys.argv[2]

name1 = os.path.basename(mcFile).split('.')[0]
name2 = os.path.basename(dsFile).split('.')[0]

assert( name1 == name2)

name = name1


mcfh = open(mcFile, 'r')
mcContent = [ line.strip() for line in list(mcfh) ]
mcfh.close()


dsfh = open(dsFile, 'r')
dsContent = [ line.strip() for line in list(dsfh) ]
dsfh.close()

if len(mcContent) < 10:
	print "Error in file: ", mcFile
	exit(-1)

if mcContent[0].find(name) == -1:
	print 'The maxcluster file may not be for this target'
	exit(-1)

if len(dsContent) < 20:
	print "Error in file: ", dsFile
	exit(-1)

if dsContent[0].find(name) == -1:
	print 'the deepscore file may not be for this target'
	exit(-1)

def ParseDeepScoreFile(dsContent):

	method1 = 'CbCb_off2'
	method1_flag = 'n5'

	method1_quality = []
	method1_first = []

	method2_quality = []
	method2_first = []

	method2 = 'ContactEC'
	method2_flag = 'UP3.0'

	for line in dsContent:
		if line.find(method1) >=0 and line.find(method1_flag) >=0:
			if line.find('model'):
				method1_quality.append( [ np.float32(v) for v in line.split()[4:8] ] )

			if line.find('model1.pdb') >=0:
				method1_first = [ np.float32(v) for v in line.split()[4:8] ]

		if line.find(method2) >=0 and line.find(method2_flag) >=0:
			if line.find('model'):
				method2_quality.append( [ np.float32(v) for v in line.split()[4:8] ] )

			if line.find('model1.pdb') >=0:
				method2_first = [ np.float32(v) for v in line.split()[4:8] ]

	method1_quality.sort(key=lambda x: (x[0]+x[2])/2 )
	method2_quality.sort(key=lambda x: (x[0]+x[2])/2 )

	return method1_first, method1_quality[-1], method2_first, method2_quality[-1]


def ParseMaxClusterFile(msContent, dsContent):
	method1 = 'CbCb_off2'
        method1_flag = 'n5'

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
mcResults = ParseMaxClusterFile(mcContent, dsContent)

finalResults = mcResults[0] + methodResults[0] + methodResults[2] + mcResults[1] + methodResults[1] + methodResults[3]

print '\t'.join([name] + [ str(v) for v in finalResults])
