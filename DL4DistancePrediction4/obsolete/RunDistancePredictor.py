import numpy as np
import os
import sys
import resource
import cPickle
import gc
import getopt

import theano
import theano.tensor as T

import config
from config import ParseResponse, Response2LabelName, GetResponseProbDims, GetResponseValueDims

import FeatureUtils
import DistanceUtils
import DataProcessor
import Model4DistancePrediction
import ContactUtils
from utils import Compatible


def Usage():
    print 'python RunDistancePredictor.py -m modelfiles -p predfiles [-d save_folder] [-g ground_truth_folder]'
    print '  -m: one or multiple deep model files in PKL format, separated by semicolon'
    print '  -p: one or multiple files containing input features in PKL format, separated by semicolon'
    print '  -d: the folder saving the result file (default current work directory), which is named after proteinName.predictedDistMatrix.pkl'
    print '	 Each file is a tuple of 6 or 7 items: proteinName, proteinSequence, predictedDistMatrixProb, predictedContactMatrix, labelWeightMatrix, labelDistributionMatrix and refDistMatrixProb'
    print '  -g: the folder for the true distance matrix saved in PKL format. When provided, contact prediction accuracy will be calculated'


def PredictDistMatrix(modelFiles, predFiles, savefolder=None):
    	## load all the models from the files. Each file contains specification for one model.
	models = []
	for mFile in modelFiles:
    		fh = open(mFile, 'rb')
    		model = cPickle.load(fh)
    		fh.close()
		if not model.has_key('Version') or model['Version'] < 3.0:
			print 'ERROR: The deep model version in the file is lower than 3.0: ', mFile
			exit(1)
		models.append(model)

	## check consistency among models. All the models shall have the same labelType for the same labelName
	labelTypes = dict()
	for model, mFile in zip(models, modelFiles):
		for response in model['responses']:
			labelName, labelType, _ = ParseResponse(response)
			if not labelTypes.has_key(labelName):
				labelTypes[labelName] = labelType
			elif labelTypes[labelName] != labelType:
				print 'ERROR: inconsistent response in the deep model: ', mFile
				exit(1)
					

	##allresults is a nested dictionary, i.e., allresults[proteinName][response] = list of predicted_prob_matrices
	##We predict one prob_matrix by each model for each protein and each response and then average them per protein and response to get the final results
	##two different models may share some common responses.

	allsequences = dict()
	allresults = dict()  ## the results predicted from the real input
	allresults_ref = dict() ## for the results predicted from reference input
	numModels = dict() ## count the number of models that may predict each response

	for model, mfile in zip(models, modelFiles):
		if not model['network'] in config.allNetworks:
			print 'ERROR: unsupported network architecture: ', model['network']
			exit(1)

		distancePredictor, x, y, xmask, ymask, xem = Model4DistancePrediction.BuildModel(model, forTrain=False)

		inputVariables = [ x, y, xmask, ymask]
		if xem is not None:
			inputVariables.append(xem)

	  	pred_prob = distancePredictor.output_prob
        	predict = theano.function(inputVariables, pred_prob, on_unused_input='warn' )

		## set model parameter values
		if not Compatible(distancePredictor.params, model['paramValues']):
			print 'FATAL ERROR: the deep model definition is not compatible with the loaded model parameters in the model file: ', mfile
			exit(1)

		[ p.set_value(v) for p, v in zip(distancePredictor.params, model['paramValues']) ]

		## We shall load these files for each model separately since each model may have different requirement of the data
		predData = DataProcessor.LoadDistanceFeatures(predFiles, modelSpecs = model, forTrainValidation=False)

		##make sure the input has the same number of features as the model. We do random check here to speed up
		FeatureUtils.CheckModelNDataConsistency(model, predData)

		## check if all the proteins of the same name have exactly the same sequence
		for d in predData:
			if not allsequences.has_key(d['name']):
				allsequences[d['name']] = d['sequence']
			elif allsequences[d['name']] != d['sequence']:
				print 'ERROR: inconsistent primary sequence for the same protein in the protein feature files'
				exit(1)
			
		## predSeqData and names are in the exactly the same order, so we know which data is for which protein	
		predSeqData = DataProcessor.SplitData2Batches(data=predData, numDataPoints=624, modelSpecs=model)
		print '#predData: ', len(predData), '#batches: ', len(predSeqData)

		##for onebatch, names4onebatch in zip(predSeqData, names):
		for minibatch in predSeqData:
			onebatch, names4onebatch = DataProcessor.AssembleOneBatch(minibatch, model)
			input = onebatch[ : len(inputVariables) ]
			result = predict(*input)
			##result is a 4-d tensor. The last dimension is the concatenation of the predicted prob parameters for all responses in this model
			assert result.shape[3] == sum( [ GetResponseProbDims(response) for response in model['responses'] ] )

			if config.UseRefState(model):
				onebatch_ref,_ = DataProcessor.AssembleOneBatch(minibatch, model, forRefState=True)
				input_ref = onebatch_ref[ : len(inputVariables) ]
				result_ref = predict(*input_ref)
				assert result_ref.shape[3] == sum( [ GetResponseProbDims(response) for response in model['responses'] ] )

			## calculate the start and end positions of each response in the last dimension of result
			dims = [ GetResponseProbDims(response) for response in model['responses'] ]
                        endPositions = np.cumsum(dims)
                        startPositions =  endPositions - dims

			for name in names4onebatch:
				if not allresults.has_key(name):
					allresults[name]=dict() 
					if config.UseRefState(model):
						allresults_ref[name]=dict() 
					numModels[name] =dict()

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

				if config.UseRefState(model):
					batchres_ref = result_ref[:, :, :, start:end ]
					revised_batchres_ref = [ probMatrix[ maxSeqLen-seqLen:, maxSeqLen-seqLen:, : ] for probMatrix, seqLen in zip(batchres_ref, seqLens) ]

					for res4one_ref, name in zip(revised_batchres_ref, names4onebatch):
                                        	if not allresults_ref[name].has_key(response):
                                                	allresults_ref[name][response] = res4one_ref
                                        	else:
                                                	## here we save sum to reduce memory consumption, which could be huge when many deep models are used to predict a large set of proteins
                                                	allresults_ref[name][response] +=  res4one_ref


		del predict
		del predData
		del predSeqData
		gc.collect()


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

	finalresults_ref = dict()
	for name, results in allresults_ref.iteritems():
		if not finalresults_ref.has_key(name):
			finalresults_ref[name] = dict()

		## finalresults has 3 dimensions. 
		for response in results.keys():
			finalresults_ref[name][response] = (allresults_ref[name][response]/numModels[name][response]).astype(np.float32)

			##make the predicted distance prob matrices symmetric for some reponses. This also slightly improves accuracy.
			labelName = Response2LabelName(response)
			if config.IsSymmetricLabel( labelName ):
				finalresults_ref[name][response] = (finalresults_ref[name][response] + np.transpose(finalresults_ref[name][response], (1, 0, 2) ) )/2.

	## collect the average label distributions and weight matrix. We collect all the matrices and then calculate their average.
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

	## convert the predicted distance probability matrix into a predicted contact matrix. 
	## Each predicted prob matrix has 3 dimensions while Each predicted contact matrix has 2 dimensions
	predictedContactMatrices = dict()
	from scipy.stats import norm
	for name, results in finalresults.iteritems():
		predictedContactMatrices[name] = dict()
		for response in results.keys():
			labelName, labelType, subType = ParseResponse(response)

			if labelName in config.allAtomPairNames:
				if labelType.startswith('Discrete'):
					subType = labelType[len('Discrete'): ]
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


	##write all the results here
	## for each protein, we have a output file saving a tuple (name, sequence, predicted distance matrix, predicted contact matrix, labelWeight, labelDistribution, ref distance prob)
	for name, results in finalresults.iteritems():

		savefilename = name + '.predictedDistMatrix.pkl'
		if savefolder is not None:
			savefilename = os.path.join(savefolder, savefilename)

		fh = open(savefilename, 'wb')

		if not finalresults_ref:
			cPickle.dump( (name, allsequences[name], results, predictedContactMatrices[name], finalLabelWeights, finalLabelDistributions), fh, protocol=cPickle.HIGHEST_PROTOCOL)
		else:
			cPickle.dump( (name, allsequences[name], results, predictedContactMatrices[name], finalLabelWeights, finalLabelDistributions, finalresults_ref[name]), fh, protocol=cPickle.HIGHEST_PROTOCOL)
		fh.close()

	#return finalresults, predictedContactMatrices, finalresults_ref, allsequences
	if not finalresults_ref:
		res4return = (predictedContactMatrices, finalresults, allsequences, finalLabelWeights, finalLabelDistributions )
	else:
		res4return = (predictedContactMatrices, finalresults, allsequences, finalLabelWeights, finalLabelDistributions, finalresults_ref)

	return res4return


def main(argv):

    	modelFiles = None
    	predFiles = None
	nativefolder = None
	savefolder = None

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
			for p in predFiles:
				if not os.path.isfile(p):
					print "input feature file does not exist: ", p
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

    	print 'modelFiles=', modelFiles
    	print 'predFiles=', predFiles
	print 'savefolder=', savefolder
	print 'nativefolder=', nativefolder

	assert len(modelFiles) > 0
	assert len(predFiles) > 0

	contPredictions = PredictDistMatrix(modelFiles, predFiles, savefolder)[0]

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

