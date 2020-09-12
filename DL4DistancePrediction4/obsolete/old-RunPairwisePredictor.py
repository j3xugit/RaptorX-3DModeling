import os
import sys
import numpy as np
import cPickle
import getopt
import glob
from scipy.stats import norm

import theano
import theano.tensor as T

import config
from config import ParseResponse, Response2LabelName, GetResponseProbDims, GetResponseValueDims

import FeatureUtils
import ContactUtils
import DistanceUtils
from Alignment import AlignmentUtils

import DataProcessor
#import Model4DistancePrediction
from Model4PairwisePrediction import BuildModel
from utils import Compatible

import gc
import resource

def Usage():
    	print 'python RunPairwisePredictor.py -m modelfiles -p proteinName or nameFile -i featureFolders [-a aliFolders] [-t tplFolder] [-d save_folder] [-g ground_truth_folder]'
	print '	This script predicts inter-atom and inter-residue relationship  (e.g., distance and orientation) for one or multiple proteins'
    	print '  -m: one or multiple deep model files in PKL format, separated by semicolon'
	print ' '
    	print '  -p: protein name or a file for protein names. A string ending with .list is interpreted as a file name; otherwise a protein name.'
	print '		when a file is specified, both featureFolders and aliFolders (if exists) shall be a single folder'
	print '		when a name is specified, both featureFolders and aliFolders (if exists) can be multiple folders'
    	print '  -i: one or multiple folders containing input features, separated by semicolon'
    	print '		Each folder has one set of input features. Three typical input feature files are: XXX.inputFeatures.pkl, XXX.extraCCM.pkl and XXX.a2m'
    	print '		Some deep models only need XXX.inputFeatures.pkl while others may need the other two types of files'
    	print '  -a: one or multiple folders containing alignment files (optional), which is needed if you want to predict from a query-template alignment instead of a sequence'
    	print '  -t: one folder for template files ending with tpl.pkl (optional), which is needed when the -a option is used'
    	print ' '
    	print '  -d: the folder saving the result files (default current work directory). multiple result files may be generated depending on input options'
    	print '	 	Each result file is a tuple of 6 items: proteinName, proteinSequence, predictedDistMatrixProb, predictedContactMatrix, labelWeightMatrix, and labelDistributionMatrix'
    	print '  -g: the folder for groundtruth files saved in PKL format. When provided, contact prediction accuracy will be calculated'


## load all the models from the files. Each file contains specification for one model.
def LoadModels(modelFiles, CheckConsistency=True):
	models = []
	for mFile in modelFiles:
    		with open(mFile, 'rb') as fh:
    			model = cPickle.load(fh)

		if not model.has_key('Version') or model['Version'] < 3.0:
			print 'ERROR: the deep model version in the file is lower than 3.0: ', mFile
			exit(1)

		if not model['network'] in config.allNetworks:
			print 'ERROR: unsupported network architecture: ', model['network']
			exit(1)

		model['modelFile'] = mFile
		models.append(model)

	if not CheckConsistency:
		return models

	## check consistency among models. All the models shall have the same labelType for the same labelName
	labelTypes = dict()
	for model in models:
		for response in model['responses']:
			labelName, labelType, _ = ParseResponse(response)
			if not labelTypes.has_key(labelName):
				labelTypes[labelName] = labelType
			elif labelTypes[labelName] != labelType:
				print 'ERROR: inconsistent response in the deep model file: ', model['modelFile']
				exit(1)

	return models
		
## calculate the average label weight and distribution matrix
def CollectLabelWeightNDistribution(models):
	labelDistributions = dict()
	labelWeights = dict()
	for model in models:
		for response in model['responses']:
			if not labelDistributions.has_key(response):
				labelDistributions[response] = []
			if not labelWeights.has_key(response):
				labelWeights[response] = []

			if model.has_key('labelDistributions'):
				labelDistributions[response].append(model['labelDistributions'][response])
			else:
				labelDistributions[response].append(model['labelRefProbs'][response])
			labelWeights[response].append(model['weight4labels'][response])

	finalLabelDistributions = dict()
	finalLabelWeights = dict()

	for response, distributions in labelDistributions.iteritems():
		finalLabelDistributions[response] = (np.average(distributions, axis=0)).astype(np.float32)
	for response, ws in labelWeights.iteritems():
		finalLabelWeights[response] = (np.average(ws, axis=0) ).astype(np.float32)

	return finalLabelWeights, finalLabelDistributions

## convert predicted distance probability matrix into contact matrix 
## A predicted distance prob matrix has 3 dimensions while a predicted contact matrix has 2 dimensions
def DeriveContactMatrix(finalresults):
	predictedContactMatrices = dict()
	for name, results in finalresults.iteritems():
		predictedContactMatrices[name] = dict()
		for response in results.keys():
			labelName, labelType, subType = ParseResponse(response)

			if labelName in config.allAtomPairNames:
				if labelType.startswith('Discrete'):
					#subType = labelType[len('Discrete'): ]
					labelOf8 = DistanceUtils.LabelsOfOneDistance(config.ContactDefinition, config.distCutoffs[subType])
					predictedContactMatrices[name][labelName] =  np.sum( finalresults[name][response][:, :, :labelOf8], axis=2)
				elif labelType.startswith('Normal'):
					assert labelType.startswith('Normal1d2')
					normDistribution =  norm( loc=finalresults[name][response][:, :, 0], scale=finalresults[name][response][:,:,1])
					predictedContactMatrices[name][labelName] =  (normDistribution.cdf(config.ContactDefinition)).astype(np.float32)
				elif labelType.startswith('LogNormal'):
					assert labelType.startswith('LogNormal1d2')
					normDistribution =  norm( loc=finalresults[name][response][:, :, 0], scale=finalresults[name][response][:,:,1])
					predictedContactMatrices[name][labelName] =  (normDistribution.cdf(np.log(config.ContactDefinition) )).astype(np.float32)
				else:
					print 'unsupported label type in response: ', response
					exit(1)

			elif labelName in ['HB', 'Beta']:
				predictedContactMatrices[name][labelName] =  finalresults[name][response][:, :, 0]

			elif labelName in config.allOrientationNames:
				continue
			else:
				print 'ERROR: unknown label name in response: ', response
				exit(1)

	return predictedContactMatrices

def FindAllTemplates(query, aliDir):
	filenamePattern = aliDir + '/*-' + query + '.fasta'
	allfiles = glob.glob(filenamePattern)
	templates = [ os.path.basename(aliFile).split('-')[0] for aliFile in allfiles ]
	return templates

def LoadProteinData4OneModel(model, names, inputFolder, aliFolder=None, tplFolder=None):

	if config.UseTemplate(model):
		if aliFolder is None or tplFolder is None:
			print 'ERROR: the deep model needs template information, but aliFolder or tplFolder is None'
			exit(1)

	data = []
	for name in names:
		rawData = dict()
		feature = FeatureUtils.LoadFeaturePKL(name, location=inputFolder, modelSpecs=model)
                rawData.update(feature)
		rawData['length'] = len(rawData['sequence'])
		rawData['name'] = name

		if not config.UseTemplate(model):
			data.append(rawData)
			continue

		templates = FindAllTemplates(query=name, aliDir=aliLocation)
		for template in templates:
			feature = AlignmentUtils.GenerateAlignmentFeatures2( (name, template), queryData=rawData, aliDir=aliFolder, tplDir=tplFolder, modelSpecs=model)
			rawData.update(feature)
			rawData['name'] = name + '-' + template
			rawData['template'] = template

			data.append(rawData)

	return data
	
def LoadOneProteinData4OneModel(model, name, inputFolders, aliFolders=None, tplFolder=None):

	if config.UseTemplate(model):
		if aliFolders is None or tplFolder is None:
			print 'ERROR: the deep model needs template information, but aliFolders or tplFolder is None'
			exit(1)

	data = []
	for featureLocation, findex in zip(inputFolders, range(len(inputFolders)) ):
		rawData = dict()
		tmpName = name + '-F' + str(findex)
		rawData['featureDir'] = featureLocation

		feature = FeatureUtils.LoadFeaturePKL(name, location=featureLocation, modelSpecs=model)
                rawData.update(feature)
		rawData['length'] = len(rawData['sequence'])
		rawData['name'] = tmpName

		if not config.UseTemplate(model):
			data.append(rawData)
			continue

		for aliLocation, aindex in zip(aliFolders, range(len(aliFolders)) ):
			templates = FindAllTemplates(query=name, aliDir=aliLocation)
			for template in templates:
				feature = AlignmentUtils.GenerateAlignmentFeatures2( (name, template), queryData=rawData, aliDir=aliLocation, tplDir=tplFolder, modelSpecs=model)
				rawData.update(feature)
				tmpName2 = tmpName + '-A' + str(aindex) + '-' + template
				rawData['name'] = tmpName2
				rawData['aliDir'] = aliLocation
				rawData['template'] = template

				data.append(rawData)

	return data
		
def BuildPredictor(model):
	#distancePredictor, x, y, xmask, ymask, xem = Model4DistancePrediction.BuildModel(model, forTrain=False)
	distancePredictor, x, y, xmask, ymask, xem = BuildModel(model, forTrain=False)

	inputVariables = [ x, y, xmask, ymask]
	if xem is not None:
		inputVariables.append(xem)

	pred_prob = distancePredictor.output_prob
        predict = theano.function(inputVariables, pred_prob, on_unused_input='warn' )

	## set model parameter values
	if not Compatible(distancePredictor.params, model['paramValues']):
		print 'FATAL ERROR: the deep model definition is not compatible with the model parameters loaded from: ', model['modelFile']
		exit(1)
	[ p.set_value(v) for p, v in zip(distancePredictor.params, model['paramValues']) ]

	return predict, inputVariables

## predict labels (distance and orientation) for a set of proteins
def PredictLabels4Proteins(models, proteinNames, inputFolder, aliFolder=None, tplFolder=None, saveFolder=None):

	##allresults is a nested dictionary, i.e., allresults[proteinName][response] = sum of predicted_prob_matrices
	##We predict one prob_matrix by each model for each protein and each response and then average them by the number of models
	##two different models may share common responses

	allsequences = dict()
	allresults = dict()  ## the predicted results
	numModels = dict() ## count the number of models that may predict each response

	for model in models:
		predict, inputVariables = BuildPredictor(model)

		## load raw data for each model separately since each model may have different requirement of the data
		rawData = LoadProteinData4OneModel(model, proteinNames, inputFolder, aliFolder, tplFolder)
		predData = DataProcessor.ExtractFeaturesNLabels(rawData, modelSpecs=model, forTrainValidation=False)

		##make sure the input has the same number of features as the model
		FeatureUtils.CheckModelNDataConsistency(model, predData)

		for d in predData:
			name = d['name']
			if not allresults.has_key(name):
				allresults[name]=dict() 
				numModels[name] =dict()

                        if not allsequences.has_key(name):
                                allsequences[name] = d['sequence']
                        elif allsequences[name] != d['sequence']:
                                print 'ERROR: inconsistent primary sequence for the same protein in the protein feature files'
                                exit(1)

		predSeqData = DataProcessor.SplitData2Batches(data=predData, numDataPoints=624, modelSpecs=model)
		print '#predData: ', len(predData), '#batches: ', len(predSeqData)

		##for onebatch, names4onebatch in zip(predSeqData, names):
		for minibatch in predSeqData:
			onebatch, names4onebatch = DataProcessor.AssembleOneBatch(minibatch, model)
			input = onebatch[ : len(inputVariables) ]
			result = predict(*input)
			##result is a 4-d tensor. The last dimension is the concatenation of the predicted prob parameters for all responses in this model
			assert result.shape[3] == sum( [ GetResponseProbDims(response) for response in model['responses'] ] )

			## calculate the start and end positions of each response in the last dimension of result
			dims = [ GetResponseProbDims(response) for response in model['responses'] ]
                        endPositions = np.cumsum(dims)
                        startPositions =  endPositions - dims

			x1d, x2d, x1dmask, x2dmask = input[0:4]
			seqLens = x1d.shape[1] - x1dmask.shape[1] + np.sum(x1dmask, axis=1)
			maxSeqLen = x1d.shape[1]

			for response, start, end in zip(model['responses'], startPositions, endPositions):

				## batchres is a batch of result, its ndim=4
				## the 1st dimension of batchres is batchSize, the 2nd and 3rd dimensions are contact/distance matrix sizes and the 4th is for the predicted probability parameters
				batchres = result[:, :, :, start:end ]
				## remove masked positions
				revised_batchres = [ probMatrix[ maxSeqLen-seqLen:, maxSeqLen-seqLen:, : ] for probMatrix, seqLen in zip(batchres, seqLens) ]

				for res4one, name in zip(revised_batchres, names4onebatch):
                                        if not allresults[name].has_key(response):
                                                allresults[name][response] = res4one
                                                numModels[name][response] = np.int32(1)
                                        else:
                                                ## here we save sum to reduce memory consumption, which could be huge when many deep models are used to predict a large set of proteins
                                                allresults[name][response] +=  res4one
                                                numModels[name][response] += np.int32(1)


	## calculate the final result, which is the average of all the predictd prob matrices for the same protein and response
	finalresults = dict()
	for name, results in allresults.iteritems():
		if not finalresults.has_key(name):
			finalresults[name] = dict()

		## finalresults has 3 dimensions. 
		for response in results.keys():
			finalresults[name][response] = (allresults[name][response]/numModels[name][response]).astype(np.float32)

			##make the predicted distance prob matrices symmetric for some reponses. This also slightly improves accuracy.
			labelName = Response2LabelName(response)
			if config.IsSymmetricLabel( labelName ):
				finalresults[name][response] = (finalresults[name][response] + np.transpose(finalresults[name][response], (1, 0, 2) ) )/2.

	## convert predicted distance probability matrix into contact matrix. 
	predictedContactMatrices = DeriveContactMatrix(finalresults)

	## collect the average label distributions and weight matrix
	finalLabelWeights, finalLabelDistributions = CollectLabelWeightNDistribution(models)

	##write all the results here
	## for each protein, we have a output file saving a tuple (name, sequence, predicted distance matrix, predicted contact matrix, labelWeight, labelDistribution)
	for name, results in finalresults.iteritems():
		savefilename = name + '.predictedDistMatrix.pkl'
		if saveFolder is not None:
			savefilename = os.path.join(saveFolder, savefilename)

		originalName = name
		for n in proteinNames:
			if name.startswith(n):
				originalName = n
				break

		with open(savefilename, 'wb') as fh:
			cPickle.dump( (originalName, allsequences[name], results, predictedContactMatrices[name], finalLabelWeights, finalLabelDistributions), fh, protocol=cPickle.HIGHEST_PROTOCOL)

	res4return = (predictedContactMatrices, finalresults, allsequences, finalLabelWeights, finalLabelDistributions )

	return res4return

## predict labels for a single proteins under different conditions, e.g., features generated from different MSAs and alignments to different templates
def PredictLabels4OneProtein(models, targetName, inputFolders, aliFolders=None, tplFolder=None, saveFolder=None):

	##allresults is a nested dictionary, i.e., allresults[proteinName][response] = sum of predicted_prob_matrices
	##We predict one prob_matrix by each model for each protein and each response and then average them per protein and response to get the final results
	##two different models may share common responses

	allsequences = dict()
	allresults = dict()  ## the results predicted from the real input
	numModels = dict() ## count the number of models that may predict each response

	for model in models:
		predict, inputVariables = BuildPredictor(model)

		## load raw data for each model separately since each model may have different requirement of the data
		rawData = LoadOneProteinData4OneModel(model, targetName, inputFolders, aliFolders, tplFolder)
		predData = DataProcessor.ExtractFeaturesNLabels(rawData, modelSpecs=model, forTrainValidation=False, returnMode='list')

		##make sure the input has the same number of features as the model 
		FeatureUtils.CheckModelNDataConsistency(model, predData)

		for d in predData:
			name = d['name']
			if not allresults.has_key(name):
				allresults[name]=dict() 
				numModels[name] =dict()

                        if not allsequences.has_key(name):
                                allsequences[name] = d['sequence']
                        elif allsequences[name] != d['sequence']:
                                print 'ERROR: inconsistent primary sequence for the same protein in the protein feature files'
                                exit(1)


		predSeqData = DataProcessor.SplitData2Batches(data=predData, numDataPoints=1000000, modelSpecs=model)
		print '#predData: ', len(predData), '#batches: ', len(predSeqData)

		##for onebatch, names4onebatch in zip(predSeqData, names):
		for minibatch in predSeqData:
			onebatch, names4onebatch = DataProcessor.AssembleOneBatch(minibatch, model)
			input = onebatch[ : len(inputVariables) ]
			result = predict(*input)
			##result is a 4-d tensor. The last dimension is the concatenation of the predicted prob parameters for all responses in this model
			assert result.shape[3] == sum( [ GetResponseProbDims(response) for response in model['responses'] ] )

			## calculate the start and end positions of each response in the last dimension of result
			dims = [ GetResponseProbDims(response) for response in model['responses'] ]
                        endPositions = np.cumsum(dims)
                        startPositions =  endPositions - dims

			x1d, x2d, x1dmask, x2dmask = input[0:4]
			seqLens = x1d.shape[1] - x1dmask.shape[1] + np.sum(x1dmask, axis=1)
			maxSeqLen = x1d.shape[1]

			for response, start, end in zip(model['responses'], startPositions, endPositions):

				## batchres is a batch of result, its ndim=4
				## the 1st dimension of batchres is batchSize, the 2nd and 3rd dimensions are distance/orientation matrix sizes and the 4th is for the predicted probability parameters
				batchres = result[:, :, :, start:end ]
				## remove masked positions
				revised_batchres = [ probMatrix[ maxSeqLen-seqLen:, maxSeqLen-seqLen:, : ] for probMatrix, seqLen in zip(batchres, seqLens) ]

				for res4one, name in zip(revised_batchres, names4onebatch):
                                        if not allresults[name].has_key(response):
                                                allresults[name][response] = res4one
                                                numModels[name][response] = np.int32(1)
                                        else:
                                                ## here we save sum to reduce memory consumption, which could be huge when many deep models are used to predict a large set of proteins
                                                allresults[name][response] +=  res4one
                                                numModels[name][response] += np.int32(1)


	## calculate the final result, which is the average of predictd prob matrices by all models for the same protein and the same response
	finalresults = dict()
	for name, results in allresults.iteritems():
		if not finalresults.has_key(name):
			finalresults[name] = dict()

		## finalresults has 3 dimensions. 
		for response in results.keys():
			finalresults[name][response] = (allresults[name][response]/numModels[name][response]).astype(np.float32)

			##make the predicted distance prob matrices symmetric for some reponses. This also slightly improves accuracy.
			labelName = Response2LabelName(response)
			if config.IsSymmetricLabel( labelName ):
				finalresults[name][response] = (finalresults[name][response] + np.transpose(finalresults[name][response], (1, 0, 2) ) )/2.

	## convert predicted distance probability matrix into contact matrix 
	predictedContactMatrices = DeriveContactMatrix(finalresults)

	## collect the average label distributions and weight matrix
	finalLabelWeights, finalLabelDistributions = CollectLabelWeightNDistribution(models)

	##write all the results here
	## for each protein, we have a output file saving a tuple (name, sequence, predicted distance matrix, predicted contact matrix, labelWeight, labelDistribution)
	for name, results in finalresults.iteritems():

		savefilename = name + '.predictedDistMatrix.pkl'
		if saveFolder is not None:
			savefilename = os.path.join(saveFolder, savefilename)

		with open(savefilename, 'wb') as fh:
			#cPickle.dump( (name, allsequences[name], results, predictedContactMatrices[name], finalLabelWeights, finalLabelDistributions), fh, protocol=cPickle.HIGHEST_PROTOCOL)
			cPickle.dump( (targetName, allsequences[name], results, predictedContactMatrices[name], finalLabelWeights, finalLabelDistributions), fh, protocol=cPickle.HIGHEST_PROTOCOL)

	res4return = (predictedContactMatrices, finalresults, allsequences, finalLabelWeights, finalLabelDistributions )

	return res4return


def main(argv):

	if len(argv)<6:
		Usage()
		exit(1)

    	modelFiles = None
    	inputFolders = None
	aliFolders = None
	tplFolder = None

	nativefolder = None
	savefolder = None

	name = None
	nameFile = None

    	try:
       		opts, args = getopt.getopt(argv, "p:i:a:t:m:g:d:", ["name=", "inputFolders=", "aliFolders=", "tplFolder=", "model=", "nativefolder=", "savefolder="])
        	#print opts, args
    	except getopt.GetoptError as err:
		print err
       		Usage()
        	exit(1)


    	if len(opts) < 3:
        	Usage()
        	exit(1)

    	for opt, arg in opts:
		if opt in ("-p", "--name"):
			if arg.endswith('.list'):
				if not os.path.isfile(arg):
					print 'ERROR: file for protein names does not exist: ', arg
					exit(1)
				nameFile = arg
			else:
				name = arg

		elif opt in ("-i", "--inputFolders"):
			inputFolders = arg.split(';')
			for f in inputFolders:
				if not os.path.isdir(f):
					print "ERROR: one input feature folder does not exist: ", f
					exit(1)

		elif opt in ("-a", "--aliFolders"):
			aliFolders = arg.split(';')
			for f in aliFolders:
				if not os.path.isdir(f):
					print "ERROR: one alignment folder does not exist: ", f
					exit(1)

		elif opt in ("-t", "--tplFolder"):
			tplFolder = arg
			if not os.path.isdir(tplFolder):
				print "ERROR: folder for tpl.pkl files does not exist: ", tplFolder
				exit(1)

		elif opt in ("-m", "--model"):
            		modelFiles = arg.split(';')
			for m in modelFiles:
				if not os.path.isfile(m):
					print "model file does not exist: ", m
					exit(1)

		elif opt in ("-d", "--savefolder"):
			savefolder = arg
			if not os.path.isdir(savefolder):
				print 'The specified folder for result save does not exist:', savefolder
				exit(1)

		elif opt in ("-g", "--nativefolder"):
			nativefolder = arg
			if not os.path.isdir(nativefolder):
				print 'The specified folder does not exist or is not accessible:', nativefolder
				exit(1)
		else:
            		Usage()
            		exit(1)

	print 'protein nameFile=', nameFile
	print 'protein name=', name
    	print 'modelFiles=', modelFiles
    	print 'inputFolders=', inputFolders
	print 'aliFolders=', aliFolders
	print 'tplFolder=', tplFolder

	print 'savefolder=', savefolder
	print 'nativefolder=', nativefolder

	assert len(modelFiles) > 0
	assert len(inputFolders) > 0

	models = LoadModels(modelFiles)

	if name is not None:
		contPredictions = PredictLabels4OneProtein(models, name, inputFolders, aliFolders, tplFolder, savefolder)[0]

	elif nameFile is not None:
		assert len(inputFolders)==1
		inputFolder = inputFolders[0]
		aliFolder = None

		if aliFolders is not None:
			assert len(aliFolders)==1
			aliFolder = aliFolders[0]

		with open(nameFile, 'r') as fh:
			names = [ n.strip() for n in list(fh) ]
		contPredictions = PredictLabels4Proteins(models, names, inputFolder, aliFolder, tplFolder, savefolder)[0]
	else:
		print 'ERROR: both name and nameFile are empty'
		exit(1)

	if nativefolder is not None:
		#print 'nativeFolder=', nativefolder
		avgacc, allacc = ContactUtils.EvaluateContactPredictions(contPredictions, nativefolder)
		ContactUtils.PrintAllContactAccuracy(avgacc, allacc)

	
if __name__ == "__main__":
    	main(sys.argv[1:])
	#sys.setrecursionlimit(40000)
	resource.setrlimit(resource.RLIMIT_STACK, [0x10000000, resource.RLIM_INFINITY])
	sys.setrecursionlimit(10000000)

	"""
	recursionlimit = sys.getrecursionlimit()
        print 'recursionlimit = ', recursionlimit
        sys.setrecursionlimit(3*recursionlimit)
	"""

