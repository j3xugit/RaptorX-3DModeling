import numpy as np
import cPickle
import os
import sys

"""
The protein sequence embedding parameters are derived from the following paper
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0141287#pone.0141287.ref010
"""

def ReadEmbeddingParamsInText(paramfile):
	f = open(paramfile, 'r')
	content = [ line.strip() for line in list(f) ]
	f.close()

	params = dict()

	for c in content:
		fields = c.split()
		assert len(fields) == 2
		assert len(fields[0]) == 3
		features = [ np.float32(x) for x in fields[1].split(',') ]
		assert len(features) == 100
		params[fields[0].upper()] = features

	assert len(params.keys()) == 9048

	return params

def SaveEmbeddingParamsInPKL(params, savefile):
	fh = open(savefile, 'wb')
	cPickle.dump(params, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	fh.close()

def LoadEmbeddingParamsInPKL(savefile):
	fh = open(savefile, 'rb')
	params = cPickle.load(fh)

	return params

def EmbedOneSequence(sequence, params):
	featureList = []
	featureList.append(params['UNK'])

	for i in range(1, len(sequence)-1 ):
		word = sequence[i-1:i+2].upper()
		if not params.has_key(word):
			featureList.append(params['UNK'])
		else:
			featureList.append(params[word])

	featureList.append(params['UNK'])

	return np.array(featureList)

if __name__ == "__main__":

	if len(sys.argv) < 2:
		print 'please specify the embedding param file in text format!'
		exit(-1)

	print 'the embedding param file is: ', sys.argv[1]
	params = ReadEmbeddingParamsInText(sys.argv[1])

	savefile = 'Mofrad-PLoSOne-2015Nov.3GramEmbeddingParams.pkl'

	SaveEmbeddingParamsInPKL(params, savefile)
