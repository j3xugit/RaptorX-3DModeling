import numpy as np
import sys
import os

if len(sys.argv)<3:
	print 'Usage: python AnalyzeMaxCluster.py MaxClusterResultFile DeepScoreResultFile'
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

if len(dsContent) < 10:
	print "Error in file: ", dsFile
	exit(-1)

if dsContent[0].find(name) == -1:
	print 'the deepscore file may not be for this target'
	exit(-1)

def ParseDeepScoreFile(dsContent):

	modelScores = dict()

	for line in dsContent:
		columns = line.split()
		modelScores[columns[-1]] =  [ np.float32(v) for v in columns[3:8] ]
		
	return modelScores


def ParseMaxClusterFile(msContent, modelScores):

	rankedModels = []
	numModels = 0
	for line in msContent[9:]:
		rankedModels.append(line.split()[6])
		numModels += 1

		if numModels >=5:
			break

	##print rankedModels

	rankedQuality = []

	for model in rankedModels:
		if modelScores.has_key(model):
			rankedQuality.append(modelScores[model])
		else:
			print 'Error: cannot find quality score for model: ', model
			exit(-1)

	quality_first = rankedQuality[0]

	rankedQuality.sort(key=lambda x: (x[1]+x[3])/2 )
	quality_best = rankedQuality[-1]

	return quality_first, quality_best
				

modelScores = ParseDeepScoreFile(dsContent)
first, best = ParseMaxClusterFile(mcContent, modelScores)

finalResults = first + best

print '\t'.join([name] + [ str(v) for v in finalResults])
