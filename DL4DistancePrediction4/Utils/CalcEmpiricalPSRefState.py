import numpy as np
import cPickle
import os
import sys
import random

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import config
from config import Response2LabelName, Response2LabelType
import ParseSimpleCommandLine as ParseCommandLine

import DataProcessor

eps=np.finfo(float).eps

def str_display(ls):
	if not isinstance(ls, (list, tuple, np.ndarray)):
		str_ls = '{0:.7f}'.format(ls)
		return str_ls
		
	str_ls = ['{0:.7f}'.format(v) for v in ls ]
	str_ls2 = '[' + ' '.join(str_ls) + ']'
	return str_ls2

## calculate the probability of labels from a set of label matrices
def CountLabels(labelMatrices, numLabels, freq=None, maxSeqSep=600):

	if freq is None:
		freq = np.zeros((maxSeqSep+1, numLabels), dtype=np.int)

       	for m in labelMatrices:
		for seqSep in range(1, min(m.shape[0], maxSeqSep+1)):
         		values = np.diagonal(m, offset=seqSep)
                	res = np.bincount(values, minlength=numLabels )
			freq[seqSep] += res

        return freq


## batch is a list of original data. Each data is a dict() and data['atomLabelMatrix'] is a dict and data['atomLabelMatrix'][response] is a label matrix
## data['seqLen'] has the length information
def CountLabels4OneBatch(batch, modelSpecs, count=None, maxSeqSep=600):
	 ## collect all discrete label matrices
        allLabelMatrices = dict()
        for response in modelSpecs['responses']:
                name = Response2LabelName(response)
                labelType = Response2LabelType(response)
                if labelType.startswith('LogNormal') or labelType.startswith('Normal'):
                        continue
                allLabelMatrices[response] = [ d['atomLabelMatrix'][response] for d in batch ]

        ## count
	if count is None:
        	count = dict()
        for response in modelSpecs['responses']:
                name = Response2LabelName(response)
                labelType = Response2LabelType(response)
                if labelType.startswith('LogNormal') or labelType.startswith('Normal'):
                        continue
                else:
			if count.has_key(response):
                        	count[response] = CountLabels(labelMatrices = [ m.astype(np.int32) for m in allLabelMatrices[response] ], numLabels = config.responseProbDims[labelType], freq=count[response], maxSeqSep=maxSeqSep )
			else:
                        	count[response] = CountLabels(labelMatrices = [ m.astype(np.int32) for m in allLabelMatrices[response] ], numLabels = config.responseProbDims[labelType], maxSeqSep=maxSeqSep )

	return count


## batches is a list of minibatches and each minibatch is a list of original data
def CalcRefState(batches, modelSpecs, maxSeqSep=600):

	count = None
	for minibatch in batches:
		count = CountLabels4OneBatch(minibatch, modelSpecs, count=count, maxSeqSep=maxSeqSep)

	refState = dict()
	for response, c in count.iteritems():
		seqLenCutoff = 400
		assert c.shape[0] > seqLenCutoff
		c2 = np.zeros((seqLenCutoff+1, c.shape[1]), dtype=np.int)
		c2[:seqLenCutoff, : ] = c[:seqLenCutoff, :]
		c2[seqLenCutoff, :] = np.sum( c[seqLenCutoff :, :], axis=0)

		s = np.sum(c2, axis=1, keepdims=True)
		refState[response] = c2*1./(s + eps)

	return refState, count	

def main(argv):

    	modelSpecs = config.InitializeModelSpecs()
	modelSpecs = ParseCommandLine.ParseArguments(argv, modelSpecs)

	## load the datasets. Data is a list of proteins and each protein is represented as a dict()
	Data = DataProcessor.LoadDistanceLabelMatrices(modelSpecs['dataset'], modelSpecs=modelSpecs )
        print '#proteins loaded from the dataset: ', len(Data)
	allProteins = [ d['name'] for d in Data ]
	maxSeqLen = max( [ d['seqLen'] for d in Data ] )

        print 'Preparing batch data for training...'
        groupSize = modelSpecs['minibatchSize']
        batches = DataProcessor.SplitData2Batches(data=Data, numDataPoints=groupSize, modelSpecs=modelSpecs)
        print "#batches:", len(batches)

	## add code here to calculate empirical reference state
	## RefState is a dict, RefState[response] 
	refState, count = CalcRefState(batches=batches, modelSpecs=modelSpecs, maxSeqSep=maxSeqLen-1 )
	refState['dataset']=modelSpecs['dataset']
	refState['proteins'] = allProteins

	## save RefState
	responseStr = '-'.join(modelSpecs['responses'])
	file4save = 'EmpiricalPSRefState-' + responseStr + '-' + str(os.getpid()) + '.pkl'
	fh = open(file4save, 'wb')
	cPickle.dump(refState, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	fh.close()

	dist = np.arange(4.25, 20, 0.5)
	## print the reference state
	for response in modelSpecs['responses']:
		for row in refState[response]:
			print str_display(row)
		for row in refState[response][6:, 1:33]:
			print str_display( np.log(row/row[0])/np.log(dist/dist[0]) )


if __name__ == "__main__":

	## generate a random seed
	a = list( str(os.getpid())  +  os.urandom(8) )
	random.shuffle(a)
	seed = ''.join(a)
	random.seed(a=seed)

   	main(sys.argv[1:])
