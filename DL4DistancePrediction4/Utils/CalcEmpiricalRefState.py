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
from utils import SampleBoundingBox

def str_display(ls):
	if not isinstance(ls, (list, tuple, np.ndarray)):
		str_ls = '{0:.4f}'.format(ls)
		return str_ls
		
	str_ls = ['{0:.4f}'.format(v) for v in ls ]
	str_ls2 = '[' + ' '.join(str_ls) + ']'
	return str_ls2

## calculate the probability of labels from a set of label matrices
def CalcLabelProb(labelMatrices, numLabels, minSeqSep=3):

	freq = []
       	for m in labelMatrices:
      		index = np.triu_indices(m.shape[0], minSeqSep)
         	values = m[index]
                res = np.bincount(values, minlength=numLabels )
                freq.append(res)

	## in case that the matrix is not symmetric
       	for m in labelMatrices:
      		index = np.tril_indices(m.shape[0], minSeqSep)
         	values = m[index]
                res = np.bincount(values, minlength=numLabels )
                freq.append(res)

        count = np.sum(freq, axis=0)
        frequency = count * 1./(np.finfo(float).eps + np.sum(count) )

	lengths = [ m.shape[0] for m in labelMatrices ]
	avgLen = np.average(lengths)

        return np.array(frequency).astype(np.float32), avgLen


## batch is a list of original data. Each data is a dict() and data['atomLabelMatrix'] is a dict and data['atomLabelMatrix'][response] is a label matrix
## data['seqLen'] has the length information
def CalcRefState4OneBatch(batch, modelSpecs, minSeqSep=3):
	 ## collect all discrete label matrices
        allLabelMatrices = dict()
        for response in modelSpecs['responses']:
                name = Response2LabelName(response)
                labelType = Response2LabelType(response)
                if labelType.startswith('LogNormal') or labelType.startswith('Normal'):
                        continue
                allLabelMatrices[response] = [ d['atomLabelMatrix'][response] for d in batch ]

        ## calculate the discrete label distribution
        allRefProbs = dict()
        for response in modelSpecs['responses']:
                name = Response2LabelName(response)
                labelType = Response2LabelType(response)
                if labelType.startswith('LogNormal') or labelType.startswith('Normal'):
                        allRefProbs[response] = np.array([1.]).astype(np.float32)
                        continue

                if modelSpecs.has_key('UseBoundingBox4RefProbs') and (modelSpecs['UseBoundingBox4RefProbs'] is True):
                        ## here we sample a sub label matrix using BoundingBox to account for the real training scenario
                        newLabelMatrices = []
                        for lMatrix in allLabelMatrices[response]:
                                bounds = SampleBoundingBox( (lMatrix.shape[0], lMatrix.shape[1]),  modelSpecs['maxbatchSize'] )
                                new_lMatrix = lMatrix[ bounds[0]:bounds[2], bounds[1]:bounds[3] ].astype(np.int32)
                                newLabelMatrices.append(new_lMatrix)
                        allRefProbs[response], avgLen = CalcLabelProb(labelMatrices = newLabelMatrices, numLabels = config.responseProbDims[labelType], minSeqSep=minSeqSep)
                else:
                        allRefProbs[response], avgLen = CalcLabelProb(labelMatrices = [ m.astype(np.int32) for m in allLabelMatrices[response] ], numLabels = config.responseProbDims[labelType], minSeqSep=minSeqSep)

	return allRefProbs, avgLen


## batches is a list of minibatches and each minibatch is a list of original data
def CalcRefState(batches, modelSpecs, minSeqSep=3):
	RefStateList = dict()
	for response in modelSpecs['responses']:
		RefStateList[response] = []

	for minibatch in batches:
		ref, length = CalcRefState4OneBatch(minibatch, modelSpecs, minSeqSep)
		for response in modelSpecs['responses']:
			RefStateList[response].append( (length, len(minibatch), ref[response]) )

	RefState = dict()
	for response in modelSpecs['responses']:
		refs = [ ref for length, count, ref in RefStateList[response] ]
		## calculate the length-ind ref
		avgRef = np.average(refs, axis=0)
		RefState[response] = (avgRef, RefStateList[response])

	return RefState	
		
	

def main(argv):

    	modelSpecs = config.InitializeModelSpecs()
	modelSpecs = ParseCommandLine.ParseArguments(argv, modelSpecs)

	## load the datasets. Data is a list of proteins and each protein is represented as a dict()
	Data = DataProcessor.LoadDistanceLabelMatrices(modelSpecs['dataset'], modelSpecs=modelSpecs )
        print '#proteins loaded from the dataset: ', len(Data)
	allProteins = [ d['name'] for d in Data ]

        print 'Preparing batch data for training...'
        groupSize = modelSpecs['minibatchSize']
        batches = DataProcessor.SplitData2Batches(data=Data, numDataPoints=groupSize, modelSpecs=modelSpecs)
        print "#batches:", len(batches)

	## add code here to calculate empirical reference state
	## RefState is a dict, RefState[response] = (length-independent ref, length-dependent ref)
	## length-independent ref is an 1d array, length-dependent ref is a list with each element being an tuple (length, 1d array)
	RefState = CalcRefState(batches=batches, modelSpecs=modelSpecs)
	RefState['dataset']=modelSpecs['dataset']
	RefState['proteins'] = allProteins

	## save RefState
	responseStr = '-'.join(modelSpecs['responses'])
	file4save = 'EmpRefState-' + responseStr + '-' + str(os.getpid()) + '.pkl'
	fh = open(file4save, 'wb')
	cPickle.dump(RefState, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	fh.close()

	## print the length-ind reference state
	for response in modelSpecs['responses']:
		print RefState[response][0]


if __name__ == "__main__":

	## generate a random seed
	a = list( str(os.getpid())  +  os.urandom(8) )
	random.shuffle(a)
	seed = ''.join(a)
	random.seed(a=seed)

   	main(sys.argv[1:])
