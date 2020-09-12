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
    	print 'python RunPairwisePredictor.py -m modelfiles -p proteinName or nameFile or inputFeature_PKL -i featureFolders [-a aliFile or aliFolders] [-t tplFile or tplFolder] [-d save_folder] [-g ground_truth_folder]'
	print '	This script predicts inter-atom and inter-residue relationship  (e.g., distance and orientation) for one or multiple proteins or one or multiple query-template alignments'
    	print '  -m: one or multiple deep model files in PKL format, separated by semicolon'
	print ' '
    	print '  -p: a string ending with .list and .inputFeatures.pkl is interpreted as a name file and a feature file, respectively; otherwise a protein name'
	print '		when a protein name is specified, featureFolders shall be one or multiple folders and the -a option may specify an alignment file or one/multiple folders for query-template alignment'
	print '		when a feature file is specified, featureFolders is ignored and the -a option may specify an alignment file or one/multiple folders for query-template alignment'
	print '		when a name file is specified, featureFolders, aliFolders (if exists) and tplFolder (if exists) shall be a single folder'
    	print '  -i: one or multiple folders containing input features, separated by semicolon'
    	print '		Each folder has one set of input features. Three typical input feature files are: XXX.inputFeatures.pkl, XXX.extraCCM.pkl and XXX.a2m'
    	print '		Some deep models use only XXX.inputFeatures.pkl while others may use the other two types of files'
	print ' '
    	print '  -a: optional, a string ending with .fasta is interpreted as an alignment file; otherwise as one or multiple folders for alignment files (optional)'
    	print '  -t: optional, a string ending with .tpl.pkl is interpreted as a template file; otherwise a template folder'
	print '		when an alignment file is specified by the -a option, this string is interpreted as a template file'
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

## load the data of multiple proteins. 
def LoadProteinData4OneModel(model, names, inputFolder, aliFolder=None, tplFolder=None):

	if not os.path.isdir(inputFolder):
		print 'ERROR: folder for protein features does not exist: ', inputFolder
		exit(1)

	if config.UseTemplate(model):
		if not os.path.isdir(aliFolder) or not os.path.isdir(tplFolder):
			print 'ERROR: the deep model needs alignment and template information, but either aliFolder or tplFolder is an invalid directory'
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

		templates = FindAllTemplates(query=name, aliDir=aliFolder)
		for template in templates:
			feature = AlignmentUtils.GenerateAlignmentFeatures2( (name, template), queryData=rawData, aliDir=aliFolder, tplDir=tplFolder, modelSpecs=model)
			rawData.update(feature)
			rawData['name'] = name + '-' + template
			rawData['template'] = template
			data.append(rawData)

	return data

## load the data of a single protein, which may have multiple sets of features when multiple methods are used to generate multiple sequence alignments	
def LoadOneProteinData4OneModel(model, name, inputFolders, aliFolders=None, tplFolder=None):

	if config.UseTemplate(model):
		if aliFolders is None:
			print 'ERROR: the deep model needs query-template alignment information, but aliFolders is None'
			exit(1)

		if not os.path.isdir(tplFolder):
			print 'ERROR: the deep model needs template information, but tplFolder is an invalid directory'
			exit(1)

	data = []
	for featureLocation, findex in zip(inputFolders, range(len(inputFolders)) ):
		rawData = dict()

		bname = os.path.basename(featureLocation)
		if bname.startswith("feat_" + name + "_"):
			bname = bname[ len("feat_" + name + "_"): ]

		"""
		if len(inputFolders) == 1:
			tmpName = name
		else:
			tmpName = name + '-' + bname
		"""

		tmpName = name
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

				if len(aliFolders) == 1:
					tmpName2 = tmpName + '-' + template
				else:
					tmpName2 = tmpName + '-A' + str(aindex) + '-' + template
				rawData['name'] = tmpName2
				rawData['aliDir'] = aliLocation
				rawData['template'] = template

				data.append(rawData)

	return data

## load the data of a single query-template alignment. it is possible that there are multiple sets of input features for the query protein.
def LoadOneAlignment4OneModel(model, name, inputFolders, aliFile, tplFile):

	data = []
	for featureLocation, findex in zip(inputFolders, range(len(inputFolders)) ):
		rawData = dict()
		rawData['featureDir'] = featureLocation
		if len(inputFolders) == 1:
			tmpName = name
		else:
			tmpName = name + '-F' + str(findex)

		feature = FeatureUtils.LoadFeaturePKL(name, location=featureLocation, modelSpecs=model)
                rawData.update(feature)
		rawData['length'] = len(rawData['sequence'])
		rawData['name'] = tmpName

		if not config.UseTemplate(model):
			data.append(rawData)
			continue

		if aliFile is None or tplFile is None:
			print 'ERROR: the deep model is trained by query-template information, but aliFile or tplFile is not provided'
			exit(1)

		feature = AlignmentUtils.GenerateAlignmentFeatures3(queryData=rawData, aliFile=aliFile, tplFile=tplFile, queryName=name, modelSpecs=model)
		rawData.update(feature)
		template = feature['template']
		rawData['name'] = tmpName + '-' + template

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

def BuildPredictors(models):
	predictors = [ ]
	for model in models:
		predict, inputVariables = BuildPredictor(model)
		predictors.append( (predict, inputVariables) )
	return predictors

## predict labels for one protein, one alignment or multiple proteins under different conditions, e.g., features generated from different MSAs and alignments to different templates
def PredictMatrixLabels(models, predictors, names, inputFolders, aliFolders=None, tplFolder=None, aliFile=None, tplFile=None, saveFolder=None):

	if not isinstance(names, (list, tuple)):
		targetName = names
	else:
		targetName = None

	##allresults is a nested dictionary, i.e., allresults[proteinName][response] = sum of predicted_prob_matrices
	##We predict one prob_matrix by each model for each protein and each response and then average them per protein and response to get the final results
	##two different models may share common responses

	allsequences = dict()
	allresults = dict()  ## the results predicted from the real input
	numModels = dict() ## count the number of models that may predict each response

	for model, predictor in zip(models, predictors):
		#predict, inputVariables = BuildPredictor(model)
		predict, inputVariables = predictor

		## load data for each model separately since each model may have a different specification
		if targetName is None:
			## in this case, inputFolders and aliFolders shall be a single folder instead of a list of folders
			rawData = LoadProteinData4OneModel(model, names, inputFolders, aliFolders, tplFolder)

		elif aliFile is not None and tplFile is not None:
			rawData = LoadOneAlignment4OneModel(model, targetName, inputFolders, aliFile, tplFile)
		else:
			rawData = LoadOneProteinData4OneModel(model, targetName, inputFolders, aliFolders, tplFolder)

		predData = DataProcessor.ExtractFeaturesNLabels(rawData, modelSpecs=model, forTrainValidation=False, returnMode='list')

		##make sure the input has the same number of features as the model 
		FeatureUtils.CheckModelNDataConsistency(model, predData)

		## check sequence consistency
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

		if targetName is not None:
			originalName = targetName
		else:
			for n in names:
				if name.startswith(n):
					originalName = n
					break

		with open(savefilename, 'wb') as fh:
			#cPickle.dump( (name, allsequences[name], results, predictedContactMatrices[name], finalLabelWeights, finalLabelDistributions), fh, protocol=cPickle.HIGHEST_PROTOCOL)
			cPickle.dump( (originalName, allsequences[name], results, predictedContactMatrices[name], finalLabelWeights, finalLabelDistributions), fh, protocol=cPickle.HIGHEST_PROTOCOL)

	return (predictedContactMatrices, allsequences)
	"""
	res4return = (predictedContactMatrices, finalresults, allsequences, finalLabelWeights, finalLabelDistributions )
	return res4return
	"""

def main(argv):

	if len(argv)<4:
		Usage()
		exit(1)

    	modelFiles = None

    	inputFolders = None

	aliFolders = None
	aliFile = None

	tplFolder = None
	tplFile = None

	nativefolder = None
	savefolder = None

	name = None
	nameFile = None
	inputFeature = None

	nameStr = None
	inputStr = None
	aliStr = None
	tplStr = None

    	try:
       		opts, args = getopt.getopt(argv, "m:p:i:a:t:g:d:", ["model=", "name=", "input=", "alignment=", "template=", "nativefolder=", "savefolder="])
        	#print opts, args
    	except getopt.GetoptError as err:
		print err
       		Usage()
        	exit(1)


    	if len(opts) < 2:
        	Usage()
        	exit(1)

    	for opt, arg in opts:
		if opt in ("-m", "--model"):
            		modelFiles = arg.split(';')
			for m in modelFiles:
				if not os.path.isfile(m):
					print "ERROR: model file does not exist: ", m
					exit(1)

		elif opt in ("-p", "--name"):
			nameStr = arg

		elif opt in ("-i", "--inputFolders"):
			inputStr = arg

		elif opt in ("-a", "--aliFolders"):
			aliStr = arg

		elif opt in ("-t", "--tplFolder"):
			tplStr = arg

		elif opt in ("-d", "--savefolder"):
			savefolder = arg
			if not os.path.isdir(savefolder):
				print 'ERROR: the specified folder for results does not exist: ', savefolder
				exit(1)

		elif opt in ("-g", "--nativefolder"):
			nativefolder = arg
			if not os.path.isdir(nativefolder):
				print 'ERROR: the specified folder for ground truth does not exist: ', nativefolder
				exit(1)
		else:
            		Usage()
            		exit(1)

	if nameStr.endswith('.inputFeatures.pkl'):
		if not os.path.isfile(nameStr):
			print 'ERROR: input feature files does not exist: ', nameStr
			exit(1)
		inputFeature = nameStr

	elif nameStr.endswith('.list'):
		if not os.path.isfile(nameStr):
			print 'ERROR: file for protein names does not exist: ', nameStr
			exit(1)
		nameFile = nameStr
	else:
		name = nameStr

	if name is not None or nameFile is not None:
		assert inputStr is not None

		inputFolders = inputStr.split(';')
		for f in inputFolders:
			if not os.path.isdir(f):
				print "ERROR: one input feature folder does not exist: ", f
				exit(1)

	if tplStr is not None and tplStr.endswith('.tpl.pkl'):
		tplFile = tplStr
		if not os.path.isfile(tplFile):
			print "ERROR: the template file does not exist: ", tplFile
			exit(1)
	elif tplStr is not None:
		tplFolder = tplStr
		if not os.path.isdir(tplFolder):
			print "ERROR: the template folder does not exist: ", tplFolder
			exit(1)

	if aliStr is not None and aliStr.endswith('.fasta'):
		aliFile = aliStr
		if not os.path.isfile(aliFile):
			print "ERROR: the alignment file does not exist: ", aliFile
			exit(1)
		if tplFile is None:
			print "ERROR: the template file is None"
			exit(1)

	elif aliStr is not None:
		aliFolders = aliStr.split(';')
		for f in aliFolders:
			if not os.path.isdir(f):
				print "ERROR: one alignment folder does not exist: ", f
				exit(1)
		if tplFolder is None:
			print 'ERROR: the template folder is None'
			exit(1)

	assert len(modelFiles) > 0
    	print 'modelFiles=', modelFiles

	print 'protein nameFile=', nameFile
	"""
	print 'protein name=', name
	print 'input feature file=', inputFeature
	"""

    	print 'inputFolders=', inputFolders

	"""
	print 'aliFolders=', aliFolders
	print 'aliFile=', aliFile

	print 'tplFolder=', tplFolder
	print 'tplFile=', tplFile
	"""

	print 'savefolder=', savefolder
	"""
	print 'nativefolder=', nativefolder
	"""

	models = LoadModels(modelFiles)
	predictors = BuildPredictors(models)

	if name is not None:
		if aliFile is not None:
			contPredictions = PredictMatrixLabels(models, predictors, name, inputFolders, aliFile=aliFile, tplFile=tplFile, saveFolder=savefolder)[0]
		else:
			contPredictions = PredictMatrixLabels(models, predictors, name, inputFolders, aliFolders=aliFolders, tplFolder=tplFolder, saveFolder=savefolder)[0]

		if nativefolder is not None:
			avgacc, allacc = ContactUtils.EvaluateContactPredictions(contPredictions, nativefolder)
			ContactUtils.PrintAllContactAccuracy(avgacc, allacc)

	elif nameFile is not None:
		assert len(inputFolders) == 1

		if aliFolders is not None:
			assert len(aliFolders)==1
			aliFolder = aliFolders[0]
		else:
			aliFolder = None

		with open(nameFile, 'r') as fh:
			names = [ n.strip() for n in list(fh) ]

		## for a batch of proteins, we predict 100 proteins every time to save CPU memory consumption
		if nativefolder is not None:
			allaccuracy = dict()
			avgaccuracy = dict()

		for i in range(0, len(names), 100):
			group = names[i : min(i+100, len(names))]
			contPredictions = PredictMatrixLabels(models, predictors, group, inputFolders[0], aliFolders=aliFolder, tplFolder=tplFolder, saveFolder=savefolder)[0]
			if nativefolder is not None:
				avgacc, allacc = ContactUtils.EvaluateContactPredictions(contPredictions, nativefolder)
				allaccuracy.update(allacc)
				for k, v in avgacc.iteritems():
					if not avgaccuracy.has_key(k):
						avgaccuracy[k] = v * len(group)
					else:
						avgaccuracy[k] += v * len(group)

		if nativefolder is not None:
			for k, v in avgaccuracy:
				avgaccuracy[k] = v/len(names)
			ContactUtils.PrintAllContactAccuracy(avgaccuracy, allacc)

	elif inputFeature is not None:
		inputFolder = os.path.dirname(inputFeature)
		name = os.path.basename(inputFeature)[ : -len('.inputFeatures.pkl') ]
		if aliFile is not None:
			contPredictions = PredictMatrixLabels(models, predictors, name, [inputFolder], aliFile=aliFile, tplFile=tplFile, saveFolder=savefolder)[0]
		else:
			contPredictions = PredictMatrixLabels(models, predictors, name, [inputFolder], aliFolders=aliFolders, tplFolder=tplFolder, saveFolder=savefolder)[0]

		if nativefolder is not None:
			avgacc, allacc = ContactUtils.EvaluateContactPredictions(contPredictions, nativefolder)
			ContactUtils.PrintAllContactAccuracy(avgacc, allacc)
	else:
		print 'ERROR: at least one of name, nameFile and inputFeature shall not be None'
		exit(1)

	"""
	if nativefolder is not None:
		avgacc, allacc = ContactUtils.EvaluateContactPredictions(contPredictions, nativefolder)
		ContactUtils.PrintAllContactAccuracy(avgacc, allacc)
	"""
	
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
