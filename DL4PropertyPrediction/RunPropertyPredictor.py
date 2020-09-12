import numpy as np
import sys
import theano
import theano.tensor as T
from numpy import random as rng
import os
import os.path
import time
import datetime
import gzip
import cPickle
import gc

import config
from config import Response2LabelType
import DataProcessor
import PropertyUtils
import Model4PropertyPrediction


import getopt

def Usage():
    print 'python RunPropertyPredictor.py -m modelfiles -p predfiles [-d save_folder] [-g ground_truth_folder]'
    print '-p: specify one or multiple files containing data to be predicted in PKL format, separated by semicolon'
    print '	 Each file is a list of proteins and each protein stored as a dictionary'
    print '-m: specify one or multiple model files in PKL format, separated by semicolon'
    print '-d: specify where to save the result files (default current work directory). The result file is named after proteinName.predictedProperties.pkl'
    print '	Each file saves a tuple of 4 items: Name, Sequence, predictedProbability (for discrete label) or mean(variance) for continuous labels, predictedValue'
    print '-g: specify the ground truth folder containing native property in PKL format'
    print '	when this option is provided, prediction accuracy will be calculated'
    print '     Each protein is stored as a dictionary of native labels or values'

def LoadModels(modelFiles):
    	## load all the models from the files. Each file contains specification for one model.
	models = []
	for mFile in modelFiles:
		if not os.path.isfile(mFile):
			print 'ERROR: the model file does not exist: ', mFile
			exit(1)

    		fh = open(mFile, 'rb')
    		model = cPickle.load(fh)
    		fh.close()
		models.append(model)

	return models

def BuildPredictors(models):
	predictors = []
	for model in models:
		if not ( model['network'].startswith('ResNet1d') or model['network'].startswith('ResNet1D') ):
			print 'ERROR: unsupproted network architecture: ', model['network']
			exit(1)

		propertyPredictor, x, xmask = Model4PropertyPrediction.BuildModel(model, forTrain=False)
		inputVariables = [ x, xmask]

	  	pred4prob = propertyPredictor.output4prob
	  	pred = propertyPredictor.y_pred
        	predict = theano.function(inputVariables, outputs=[pred4prob, pred], on_unused_input='warn' )

		## set model parameter values
		[ p.set_value(v) for p, v in zip(propertyPredictor.params, model['paramValues']) ]

		predictors.append( (predict, inputVariables) )

	return predictors

def PredictProperty(models, predictors, predFiles):

	allsequences = dict()

	##allresults shall be a nested dictionary, e.g, allresults[proteinName][response] = predicted_property_list
	allresults4prob = dict()
	allresults = dict()

	for model, predictor in zip(models, predictors):

		predict, inputVariables = predictor

		## We shall load these files for each model separately since each model may use a different set of features
		predData = DataProcessor.LoadPropertyFeatures(predFiles, modelSpecs=model, forTrainValidation=False)

		##make sure the input has the same number of features as the model
		rindex = np.random.randint(0, high=len(predData) )
		assert model['n_in_seq'] == predData[rindex]['seqFeatures'].shape[1]

		## collecting sequences
		for d in predData:
			if not allsequences.has_key(d['name']):
				allsequences[d['name']] = d['sequence']
			elif allsequences[d['name']] != d['sequence']:
				print 'ERROR: inconsistent primary sequence for the same protein in the protein feature files'
				exit(1)
				
		predSeqData, names = DataProcessor.SplitData2Batches(data=predData, numDataPoints=30, modelSpecs=model, forTrainValidation=False)
		print '#predData: ', len(predData), '#batches: ', len(predSeqData)

		for onebatch, names4onebatch in zip(predSeqData, names):
			input = onebatch[ : len(inputVariables) ]
			result4prob, result = predict(*input)

			## x1d has shape (batchSize, maxSeqLen, numFeatures) and x1dmask has shape (batchSize, #cols_to_be_masked)
			x1d, x1dmask = input[0:2]
			seqLens = x1d.shape[1] - x1dmask.shape[1] + np.sum(x1dmask, axis=1)
			maxSeqLen = x1d.shape[1]

			##result4prob has shape (batchSize, maxSeqLen, sum( responseProbDims{res] for res in modelSpecs['responses'])  )
			assert result4prob.shape[2] == sum( [ config.responseProbDims[ Response2LabelType(res) ] for res in model['responses'] ] )

			##result has shape (batchSize, maxSeqLen, sum( responseValueDims{res] for res in modelSpecs['responses'])  )
			assert result.shape[2] == sum( [ config.responseValueDims[ Response2LabelType(res) ] for res in model['responses'] ] )

			nameGenerator = ( name for name in names4onebatch if not allresults.has_key(name) )
			for name in nameGenerator:
				allresults[name]=dict() 
				allresults4prob[name]=dict() 

			dims = [ config.responseProbDims[ Response2LabelType(res) ] for res in model['responses'] ]
			endPositions = np.cumsum(dims)
			startPositions =  endPositions - dims

			for res, start, end in zip(model['responses'], startPositions, endPositions):
				nameGenerator = ( name for name in names4onebatch if not allresults4prob[name].has_key(res) )
				for name in nameGenerator:
					allresults4prob[name][res] = []	

				## remove masked positions
				revised_batchres = [ tmp[ maxSeqLen-seqLen:, : ] for tmp, seqLen in zip(result4prob[:,:,start:end], seqLens) ]

				[ allresults4prob[name][res].append( res4one ) for res4one, name in zip(revised_batchres, names4onebatch) ]

			dims = [ config.responseValueDims[ Response2LabelType(res) ] for res in model['responses'] ]
			endPositions = np.cumsum(dims)
			startPositions =  endPositions - dims

			for res, start, end in zip(model['responses'], startPositions, endPositions):
				nameGenerator = ( name for name in names4onebatch if not allresults[name].has_key(res) )
				for name in nameGenerator:
					allresults[name][res] = []	

				## remove masked positions
				revised_batchres = [ tmp[ maxSeqLen-seqLen:, : ] for tmp, seqLen in zip(result[:,:,start:end], seqLens) ]
				[ allresults[name][res].append( res4one ) for res4one, name in zip(revised_batchres, names4onebatch) ]


	## calculate the final result, which is the average of all the predictd properties for the same protein and response name
	finalresults = dict()
	for name, results in allresults.iteritems():
		if not finalresults.has_key(name):
			finalresults[name] = dict()
		for response in results.keys():
			tmpresult = np.average(allresults[name][response], axis=0)

			##convert coding of discrete labels to more meaningful representation
			labelType = Response2LabelType(response)
			if not labelType.startswith('Discrete'):
				finalresults[name][response] = tmpresult

	finalresults4prob = dict()
	for name, results in allresults4prob.iteritems():
		if not finalresults4prob.has_key(name):
			finalresults4prob[name] = dict()
		for response in results.keys():
			finalresults4prob[name][response] = np.average(allresults4prob[name][response], axis=0)

			labelType = Response2LabelType(response)
			if labelType.startswith('Discrete'):
				tmpresult = np.argmax(finalresults4prob[name][response], axis=1)
				finalresults[name][response] = PropertyUtils.Coding2String(tmpresult, response)

	"""
	## collect the average label distributions and weight matrix. We collect all the matrices and then calculate their average.
	labelDistributions = dict()
	labelWeights = dict()
	for model in models:
		for apt in model['responseNames']:
			if not labelDistributions.has_key(apt):
				labelDistributions[apt] = []
			if not labelWeights.has_key(apt):
				labelWeights[apt] = []

			labelDistributions[apt].append(model['labelRefProbs'][apt])
			labelWeights[apt].append(model['weight4' + model['labelType'] ][apt])

	finalLabelDistributions = dict()
	finalLabelWeights = dict()

	for apt in labelDistributions.keys():
		finalLabelDistributions[apt] = np.average(labelDistributions[apt], axis=0)
	for apt in labelWeights.keys():
		finalLabelWeights[apt] = np.average(labelWeights[apt], axis=0)
	"""

	return finalresults4prob, finalresults, allsequences

def PredictProperty2(modelFiles, predFiles):
	models = LoadModels(modelFiles)
	predictors = BuildPredictors(models)
	return PredictProperty(models, predictors, predFiles)

def SaveResults(finalresults4prob, finalresults, allsequences, savefolder=None):
	##write all the results here
	for name, results in finalresults.iteritems():
		savefilename = name + '.predictedProperties.pkl'
		if savefolder is not None:
			savefilename = os.path.join(savefolder, savefilename)

		fh = open(savefilename, 'wb')
		cPickle.dump( (name, allsequences[name], finalresults4prob[name], finalresults[name] ), fh, protocol=cPickle.HIGHEST_PROTOCOL)
		fh.close()

	#return finalresults4prob, finalresults, allsequences

def main(argv):

    	modelFiles = None
    	predFiles = None
	nativefolder = None
	savefolder = None

	if len(argv) < 4:
		Usage()
		exit(1)

    	try:
       		opts, args = getopt.getopt(argv, "p:m:g:d:", ["predictfile=", "model=", "nativefolder=", "savefolder="])
        	#print opts, args
    	except getopt.GetoptError as err:
		print err
       		Usage()
        	exit(1)

    	if len(opts) < 2:
        	Usage()
        	exit(1)

    	for opt, arg in opts:
		if opt in ("-p", "--predictfile"):
           		 predFiles = arg.split(';')
		elif opt in ("-m", "--model"):
            		modelFiles = arg.split(';')
		elif opt in ("-d", "--savefolder"):
			savefolder = arg
			if not os.path.isdir(savefolder):
				print 'ERROR: the folder for result saving does not exist:', savefolder
				exit(1)

		elif opt in ("-g", "--nativefolder"):
			nativefolder = arg
			if not os.path.isdir(nativefolder):
				print 'ERROR: the folder for native properties does not exist or is not accessible:', nativefolder
				exit(1)
		else:
            		Usage()
            		exit(1)

    	print 'modelFiles=', modelFiles
    	print 'predFiles=', predFiles
	print 'savefolder=', savefolder
	print 'nativefolder=', nativefolder

	assert len(modelFiles) > 0
	assert len(predFiles) > 0

	finalresults4prob, finalresults, allsequences = PredictProperty2(modelFiles, predFiles)

	SaveResults(finalresults4prob, finalresults, allsequences, savefolder)

	if nativefolder is not None:
		#print 'nativeFolder=', nativefolder
		propertyPrediction = finalResults
		avgErrPerTarget, avgErrPerResidue, allerrors = PropertyUtils.EvaluatePropertyPrediction(propertyPrediction, nativefolder)

		for response, err in avgErrPerTarget.iteritems():
			print 'average prediction error per target for ', response
			print err
			print 'average prediction error per residue for ', response
			print avgErrPerResidue[response]

	
if __name__ == "__main__":
    main(sys.argv[1:])

