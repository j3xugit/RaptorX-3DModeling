import numpy as np
import cPickle
import os
import os.path
import sys
import time
import datetime
import random
import gzip
import theano

import config
from config import Response2LabelType, Response2LabelName
import PropertyUtils

from Common import SequenceEmbedding

##this file contains some functions for data processing
def LoadNativeLabelsFromFile(filename):
	if not os.path.isfile(filename):
		print 'ERROR: cannot find the file for native property labels: ', filename
		exit(1)

	with open(filename, 'rb') as fh:
		property = cPickle.load(fh)
	return property

##this file contains some functions for data processing
def LoadNativeLabels(name, location='pdb25-7952-properties/', responses=None):
	filename = os.path.join(location, name+'.nativeProperties.pkl')
	if not os.path.isfile(filename):
		print 'WARNING: cannot find the file for native property labels: ', filename
		return None

	with open(filename, 'rb') as fh:
		property = cPickle.load(fh)
	return property

def LoadPropertyFeatures(files=None, modelSpecs=None, forTrainValidation=True):
    	if files is None or len(files)==0:
       		print 'ERROR: the feature files is empty'
		exit(1)

	data = []
	for infile in files:
		with open(infile, 'rb') as fh:
			data.extend(cPickle.load(fh) )

	EmbeddingModel = None
	if modelSpecs.has_key('UseSequenceEmbedding') and modelSpecs['UseSequenceEmbedding']:
		EmbeddingModelFile = os.path.join( os.environ['DL4PropertyPredHome'], 'data', 'Mofrad-PLoSOne-2015Nov.3GramEmbeddingParams.pkl')
		EmbeddingModel = SequenceEmbedding.LoadEmbeddingParamsInPKL( EmbeddingModelFile )

    	## each protein has sequential features as input 
    	proteinFeatures=[]
    	counter = 0

   	for d in data:
        	oneprotein = dict()
        	oneprotein['name'] = d['name']

		##collecting sequential features...
		seqMatrices = []

        	seqMatrices.append( d['PSSM'] )
        	##seqMatrices.append( d['PSFM'] )

		##Load sequence embedding features here
		if EmbeddingModel is not None:
			seqMatrices.append ( SequenceEmbedding.EmbedOneSequence( d['sequence'], EmbeddingModel ) )

		if modelSpecs.has_key('UsePSFM') and modelSpecs['UsePSFM']:
        		seqMatrices.append( d['PSFM'] )

		if modelSpecs.has_key('UseOneHotEncoding') and modelSpecs['UseOneHotEncoding']:
        		seqMatrices.append( config.SeqOneHotEncoding(d['sequence']) )

		## add template similarity score here
		if modelSpecs.has_key('UseTemplate') and modelSpecs['UseTemplate']:
	    		#print 'Using template similarity score...'
			if not d.has_key('tplSimScore'):
				print 'ERROR: no tplSimScore for target', d['name'], 'which is needed since you specify to use template information'
				exit(1)
			if d['tplSimScore'].shape[1] != 10:
				print 'ERROR: the number of query-template similarity features is not 10 in data for', d['name']
				exit(1)

			if not d.has_key('tplProperties'):
				print 'ERROR: no tplProperties for target', d['name'], 'which is needed since you specify to use template information'
				exit(1)

			if d['tplProperties'].shape[1] < 15:
				print 'ERROR: #template local structure properties shall be at least 15 for target', d['name']
				exit(1)

			## the query-template similarity score shall be arranged in the order of: AA identity (binary), blosum80, blosum62, blosum45, spScore, spScore_ST, ppScore, pmScore, cc, hdsm
	    		seqMatrices.append( d['tplSimScore'] )

			##we do not use omg information from the template, the first 8 features shall be the 8-state secondary structure, then followed by pACC, CNa, CNb, Phi, Psi, Theta and Tau
	    		#seqMatrices.append( d['tplProperties'][:,:15] )
	    		seqMatrices.append( d['tplProperties'][:,:8] )
			for r in modelSpecs['responses']:
				if r.startswith('ACC'):
	    				seqMatrices.append( d['tplProperties'][:,8:9] )
				elif r.startswith('Phi') or r.startswith('Psi') or r.startswith('CLE'):
	    				seqMatrices.append( d['tplProperties'][:,11:13] )
				elif r.startswith('Theta') or r.startswith('Tau'):
	    				seqMatrices.append( d['tplProperties'][:,13:15] )
				elif r.startswith('CNa') or r.startswith('CNb'):
	    				seqMatrices.append( d['tplProperties'][:,9:11] )
				else:
					print 'ERROR: unsupported response', r
					exit(1)

		if d.has_key('otherSeqFeatures'):
	    		seqMatrices.append(d['otherSeqFeatures'])

		## all the features shall have shape (seqLen, nFeatures) where nFeatures is variable, but seqLen is the sequence length of one protein
		seqFeature = np.concatenate( seqMatrices, axis=1).astype(np.float32)

		oneprotein['sequence'] = d['sequence']
		oneprotein['seqLen'] = seqFeature.shape[0]
		oneprotein['seqFeatures'] = seqFeature

		if not d.has_key('DISO') and d.has_key('Missing'):
			d['DISO'] = d['Missing']

		##collecting labels...
		for r in modelSpecs['responses']:
			labelName = Response2LabelName(r)
			labelType = Response2LabelType(r)

			if not d.has_key(labelName) and forTrainValidation:
				print 'ERROR: missing label information for protein ', d['name'], ' and response ', r
				exit(1)
			elif not d.has_key(labelName):
				continue

			labels = d[labelName]

			## need some special handling of discrete labels
			if labelType.startswith('Discrete'):
				if r.startswith('SS3'):
					labels = np.array( [ PropertyUtils.SS3Letter2Code[c] for c in labels ] ).reshape( (-1, 1) )
				elif r.startswith('SS8'):
					labels = np.array( [ PropertyUtils.SS8Letter2Code[c] for c in labels ] ).reshape( (-1, 1) )
				elif r.startswith('ACC') or r.startswith('DISO'):
					labels = labels.reshape( (-1, 1) )
				elif r.startswith('CLE'):
					labels = np.array( [ PropertyUtils.CLELetter2Code[c] for c in labels ] ).reshape( (-1, 1) )
				else:
					print 'ERROR: please specify how to convert your discrete labels to numbers for response ', r
					exit(1)

			oneprotein[labelName] = labels
				
		##at this point, finish collecting features and labels for one protein
		if d.has_key('Missing'):
			oneprotein['missing'] = d['Missing']
		elif forTrainValidation:
			print 'ERROR: for training data, we need information to specify which residues have no 3D coordinates'
			exit(1)

        	proteinFeatures.append(oneprotein)

		counter += 1
		if (counter %500 ==1):
            		print 'assembled features and labels for ', counter, ' proteins.'

    	"""
    	tmpfile = open(files[0] + '.contactInput.pkl', 'wb')
    	cPickle.dump(proteinFeatures, tmpfile, protocol = cPickle.HIGHEST_PROTOCOL)
    	tmpfile.close()
    	"""

    	return proteinFeatures


## this function prepares one batch of data for training, validation and test
## data is a list of protein features and possibly labels, generated by LoadDistanceFeatures
def AssembleOneBatch( data, modelSpecs, forTrainValidation=True ):
	if not data:
		print 'WARNING: the list of data is empty'
		return None

	numSeqs = len(data)
	seqLens = [ d['seqLen'] for d in data ]
	maxSeqLen = max( seqLens )
	minSeqLen = min( seqLens )
	#print 'maxSeqLen= ', maxSeqLen, 'minSeqLen= ', minSeqLen

        X1d = np.zeros(shape=(numSeqs, maxSeqLen, data[0]['seqFeatures'].shape[1] ), dtype = theano.config.floatX)
	## for mask
	M1d = np.zeros(shape=(numSeqs, maxSeqLen - minSeqLen ), dtype=np.int8 )

	## Y shall be a list of labels, each for one type
	##we always need a weight vector to deal with residues without 3D coordinates in training and validation, if modelSpecs['UseSampleWeight']:
        Y = []
	W = []
	for res in modelSpecs['responses']:
		labelType = Response2LabelType(res)
		labelName = Response2LabelName(res)

		dataType = ( np.int32 if labelType.startswith('Discrete') else theano.config.floatX )
		if forTrainValidation:
			if not data[0].has_key(labelName):
				print 'ERROR: label information is needed for training protein ', data['name'], ' and response ', res
				exit(1)
            		Y.append(np.zeros( shape=(numSeqs, maxSeqLen, config.responseValueDims[labelType] ), dtype = dataType )  )

			if not data[0].has_key('missing'):
				print 'ERROR: missing information is needed for training protein ', data['name']
				exit(1)

            		W.append(np.zeros( shape=(numSeqs, maxSeqLen, 1), dtype = theano.config.floatX )  )

        for j in range(len(data) ):
            	seqLen = data[j]['seqLen']
            	X1d[j, maxSeqLen - seqLen :, : ] = data[j]['seqFeatures']
            	M1d[j, maxSeqLen - seqLen : ].fill(1)

		for y, w, res in zip(Y, W, modelSpecs['responses']):
                	y[j, maxSeqLen - seqLen :, ] = data[j][Response2LabelName(res) ]

			if res.startswith('DISO'):
				## for disorder prediction, all the residues shall be considered since those residues without 3D coordinates are positive examples
				## we may assign a larger weight to positive examples since they are only 6% of the whole data set
                		w[j, maxSeqLen - seqLen :, ] = np.reshape(data[j]['missing'], (-1, 1) ) * (modelSpecs['w4diso'] -1.) + 1.
			else:
				## assign weight 0 to those residues without coordinates, otherwise 1
                		w[j, maxSeqLen - seqLen :, ] = 1.0 - np.reshape(data[j]['missing'], (-1, 1) )

	onebatch = [X1d, M1d]
	onebatch.extend(Y)
	onebatch.extend(W)

	return onebatch

##split data into minibatch, each minibatch numDataPoints data points
def SplitData2Batches(data=None, numDataPoints=1000000, modelSpecs=None, forTrainValidation=True):

    	if data is None:
        	print 'ERROR: Please provide data for process!'
		exit(1)

    	if numDataPoints < 10:
		print 'ERROR: Please specify the number of data points in a minibatch'
		exit(1)

    	## sort proteins by length from large to small
    	data.sort(key=lambda x: x['seqLen'], reverse=True)

	##seqDataset stores the resultant data
    	batches = []
	names = []

    	i = 0
    	while i < len(data):

        	currentSeqLen = data[i]['seqLen']
		numSeqs = min( len(data) - i, max(1, numDataPoints/currentSeqLen ) )
		#print 'This batch contains ', numSeqs, ' sequences'

		names4onebatch = [ d['name'] for d in data[i: i+numSeqs] ]
		oneBatch = AssembleOneBatch( data[i : i+numSeqs], modelSpecs, forTrainValidation )
		batches.append(oneBatch)
		names.append(names4onebatch)
	
		i += numSeqs

    	return batches, names

