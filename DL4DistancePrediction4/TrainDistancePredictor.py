import os
import sys
import time
import datetime
import random
import numpy as np

import cPickle
import json

import multiprocessing
import threading

import theano
import theano.tensor as T

from shared_ndarray import SharedNDArray

import config
from config import ParseResponse, GetResponseProbDims, GetResponseValueDims
import ParseCommandLine

#from Model4DistancePrediction import BuildModel
from Model4PairwisePrediction import BuildModel

from Optimizers import SGDM, SGDM2, Nesterov
from Adams import Adam, AdamW, AdamWAMS, AMSGrad

import TrainUtils
import FeatureUtils
import LabelUtils
import DataProcessor

from utilsNoT import SampleBoundingBox, Compatible, str_display, SplitList
from utils import ToFloatX
from Initialize import InitializeModelSpecs

##a small function to convert shared ndarray to np.ndarray
## d is a list, in which some elements are shared ndarray
def ToNonSharedArray(d):
	newD = []
	for e in d:
		if isinstance(e, SharedNDArray):
			newD.append(e.array)
		else:
			newD.append(e)
	return newD

## data is a list of arrays, some of which are SharedNDArray
def FreeSharedArrayInOneBatch(data):
	if data is None:
		return
	for e in data:
		if isinstance(e, SharedNDArray):
			e.unlink()

## data is a list of batches. Each batch is a list of arrays, some of which are SharedNDArray
def FreeSharedArrayInBatches(data):
	if data is None:
		return
	for batch in data:
		FreeSharedArrayInOneBatch(batch)

##calculate the weight for each minibatch where weightList is the weight matrix and base is the normalization constant
##weightList is a list of tensors, each having shape(batchSize, seqLen, seqLen)
##base is an 1d array, with the same length as weightList
def CalcBatchWeight(weightList, base):
	if weightList is None or len(weightList)<=0:
		return None
	assert len(weightList) == base.shape[0]
	weight = [ T.sum(w)*1./b for w, b in zip(weightList, base) ]

	return T.stack(weight)

def ScaleLossByBatchWeight(loss, weightList, modelSpecs):
	weightedLoss = loss
        if weightList is not None and len(weightList)>0:
		if modelSpecs.has_key('batchWeightBase'):
			batchWeight = CalcBatchWeight(weightList, modelSpecs['batchWeightBase'])
			weightedLoss = T.mul(loss, batchWeight)

	return weightedLoss

## pdecay is the decay of params. It is only used for AdamW and AdamWAMS
## if pdecay is None, then in AdamW and AdamWAMS, pdecay is set to params
## Only SGDM, Adam and AdamW are fully tested
def UpdateAlgorithm(alg, params, gparams, pdecay, lr=None, l2reg=None):

	other_params = []

    	if alg == 'SGD':
        	updates = GD(params, gparams, lr)

    	elif alg == 'SGDM':
        	updates, other_params = SGDM(params, gparams, np.float32(0.90), lr)

    	elif alg == 'SGDM2':
        	updates, other_params = SGDM2(params, gparams, np.float32(0.90), lr)

	elif alg == 'SGNA':
		updates, other_params = Nesterov(params, gparams, np.float32(0.90), lr)

    	elif alg == 'Adam':
	    	updates, other_params = Adam(params, gparams, lr = lr)

	elif alg == 'AMSGrad':
	    	updates, other_params = AMSGrad(params, gparams, lr = lr)

	elif alg == 'AdamW':
		assert (l2reg is not None)
	    	updates, other_params = AdamW(params, gparams, pdecay=pdecay, l2reg=l2reg, lr = lr)

	elif alg == 'AdamWAMS':
		assert (l2reg is not None)
	    	updates, other_params = AdamWAMS(params, gparams, pdecay=pdecay, l2reg=l2reg, lr = lr)

    	else:
        	raise NotImplementedError

	return updates, other_params

## Initialize the check point. We do not set the values of other_params since when this function is called, other_params may not be defined
def InitializeChkpoint(params, modelSpecs):
        chkpoint = dict()
        chkpoint['best_validation_loss'] = np.array([np.inf ] * len(modelSpecs['responses']) )
        chkpoint['bestParamValues'] = None
        chkpoint['bestOtherParamValues'] = None
        chkpoint['epoch'] = 0
        chkpoint['best_epoch'] = 0

        checkPointFile = None
        if modelSpecs.has_key('checkpointFile'):
                checkPointFile = modelSpecs['checkpointFile']

	restart = False
        if checkPointFile is not None and os.path.isfile(checkPointFile):
                print 'restarting from the check point file: ', checkPointFile
                saved_fh = open(checkPointFile, 'rb')
                chkpoint = cPickle.load( saved_fh )
                saved_fh.close()

                if not chkpoint.has_key('param_values') or not chkpoint.has_key('other_param_values'):
                        print 'ERROR: the check point file does not have param_values or other_param_values!'
                        exit(1)

                if Compatible(params, chkpoint['param_values']):
                        [ p.set_value(v) for p, v in zip(params, chkpoint['param_values']) ]
                else:
                        print 'ERROR: param_values in the chkpoint file is incompatible with params!'
                        exit(1)

		restart = True

        else:
                print 'Could not find the check point file: ', checkPointFile
                print 'Training the model from scratch ...'

        return chkpoint, restart

## calculate loss, error and accuracy for a set of data
## validate is a theano function that executes the validation procedure
## validate has two possible implementations: fullvalidate or quickvalidate

## batches is a list of minibatch, each containing a list of protein information
def ValidateAllData(validData, validate, modelSpecs, forRefState=False):
	accs = []
	losses = []
	errs = []
        numSamples = []

        if config.UseSampleWeight(modelSpecs):
                w4losses = []
                w4errors = []
        else:
                w4losses = None
                w4errors = None


        for batch in validData:
		
		input = ToFloatX( ToNonSharedArray(batch) )
		onebatch = input[:-1]

	       	onebatch_res = validate( *input )
		los = onebatch_res[0]
		err = onebatch_res[1]
		losses.append(los)
	    	errs.append(err)

		if len(onebatch_res) > 2:
			acc = onebatch_res[2]
	    		accs.append( acc )
			##numSamples is the number of proteins in one batch
	    		numSamples.append(onebatch[0].shape[0])

                if config.UseSampleWeight(modelSpecs):
                        #weights = onebatch[ len(onebatch) - len(modelSpecs['responses']) : ]
                        weights = onebatch[-len(modelSpecs['responses']) : ]
                        w4loss = []
                        w4error = []
                        for res, w in zip(modelSpecs['responses'], weights):
                                wSum = np.sum(w)
                                w4loss.append(wSum)
                                w4error.extend( [ wSum] * GetResponseValueDims(res) )
                        w4losses.append(w4loss)
                        w4errors.append(w4error)

	## The loss and err is normalized by the weight of each minibatch. This is equivalent to minimize loss and err per residue pair
	## The top accuracy is not normalized by the weight of a minibatch, i.e.,  we want to maximize per-protein accuracy.
	if len(accs)>0 and len(numSamples)>0 :
        	return np.average(losses, axis=0, weights=w4losses), np.average(errs, axis=0, weights=w4errors), np.average(accs, axis=0, weights=numSamples)
	else:
        	return np.average(losses, axis=0, weights=w4losses), np.average(errs, axis=0, weights=w4errors)

def TrainByOneBatch(batch, train, modelSpecs, forRefState=False):

	## batch is a list of protein locations, so we need to load the real data here
	minibatch = DataProcessor.LoadRealData(batch, modelSpecs)

	## add code here to make sure that the data has the same input dimension as the model specification
	FeatureUtils.CheckModelNDataConsistency(modelSpecs, minibatch)

	onebatch, names4onebatch = DataProcessor.AssembleOneBatch(minibatch, modelSpecs, forRefState=forRefState)
        x1d, x2d, x1dmask, x2dmask = onebatch[0:4]

	## crop a large protein to deal with limited GPU memory. For sequential and embedding features, the theano model itself will crop based upon bounding box
        bounds = SampleBoundingBox( (x2d.shape[1], x2d.shape[2]), modelSpecs['maxbatchSize'] )

        #x1d_new = x1d[:, bounds[1]:bounds[3], :]
        x1d_new = x1d
        x2d_new = x2d[:, bounds[0]:bounds[2], bounds[1]:bounds[3], :]
        #x1dmask_new = x1dmask[:, bounds[1]:x1dmask.shape[1] ]
        x1dmask_new = x1dmask
        x2dmask_new = x2dmask[:, bounds[0]:x2dmask.shape[1], bounds[1]:bounds[3] ]

	input = [x1d_new, x2d_new, x1dmask_new, x2dmask_new]

	## if embedding is used
	##if any( k in modelSpecs['seq2matrixMode'] for k in ('SeqOnly', 'Seq+SS') ):
	if config.EmbeddingUsed(modelSpecs):
              	embed = onebatch[4]
                #embed_new = embed[:, bounds[1]:bounds[3], : ]
                embed_new = embed
                input.append(embed_new)

	      	remainings = onebatch[5:]
        else:
		remainings = onebatch[4:]


        ##crop the ground truth and weight matrices
	for x2d0 in remainings:
		if len( x2d0.shape ) == 3:
			input.append( x2d0[:, bounds[0]:bounds[2], bounds[1]:bounds[3] ] )
		else:
			input.append( x2d0[:, bounds[0]:bounds[2], bounds[1]:bounds[3], : ] )

	## add bounding box to the input list
	input.append(bounds)

	if config.TrainByRefLoss(modelSpecs):
		if forRefState:
			input.append(np.int32(-1) )
		else:
			input.append(np.int32(1) )

        train_loss, train_errors, param_L2 = train(*input)

	return train_loss, train_errors, param_L2

## train the model for one epoch. Each epoch scan the whole list of sampled training proteins once.
def RunOneEpoch(epoch, data, chkpoint, params, other_params, train, quickValidate, fullValidate, modelSpecs):
	print '\nstart time of epoch ', epoch, ': ', datetime.datetime.now()

	## trainSeqData and validSeqData are a list of batches. Each batch contain file paths for feature and labels of a protein
	## the real content of a training/validation protein will be loaded right before running training and validation
	trainSeqData, validSeqData = data

        #random.shuffle(trainSeqData), 
	## we just use trainSeqData to approximately estimate the number of batches. the real train data is loaded by another process
	n_train_batches = len(trainSeqData)

        validation_frequency = min( max( len(validSeqData), len(trainSeqData)/12 ), len(trainSeqData)/6)
	validation_frequency *= max(0.8, min(3.5, 16.0/(epoch+0.0001)))
	validation_frequency = np.int32(validation_frequency)

	chkpt_save_frequency = max(1, min(validation_frequency, np.int32( len(trainSeqData) / 30 ) ) )

        t_loss = []
        t_errors = []
        results = []

        for minibatch_index in xrange( n_train_batches ):
		#input = TrainUtils.FetchOneBatch(sharedQ)
		#input = TrainUtils.FetchOneBatch(sharedQ, modelSpecs)
		global trainSharedQ, trainSharedLabelPool, trainSharedLabelWeightPool
		d = trainSharedQ.get()

		#print 'convert shared array to non-shared array...'
		input = ToFloatX(ToNonSharedArray(d) )

		#print 'feed non-shared array to train...'
		train_loss, train_errors, param_L2 = train(*input)

		#print 'free shared array after train...'
		FreeSharedArrayInOneBatch(d)

                t_loss.append(train_loss)
                t_errors.append(train_errors)

		""" need some revision here if TrainByRefLoss is enabled
		if config.TrainByRefLoss(modelSpecs):
			## add code here to train by a reference input
			_, _, param_L2 = TrainByOneBatch(minibatch, train, modelSpecs, forRefState=True)
		"""

		## validate the model	
                if ( minibatch_index + 1 ) % validation_frequency  == 0:
			global validData, validSharedQ, validDataLoader, stopValidDataLoader
			if validData is None:
				print 'Fetching all validation data at ', datetime.datetime.now()
				validData = [ validSharedQ.get() for i in range(len(validSeqData) ) ]
				stopValidDataLoader.set()
				validDataLoader.terminate()
				validDataLoader.join()	
				print 'Finish fetching validation data at ', datetime.datetime.now()

			valid_loss, valid_errors = ValidateAllData(validData, quickValidate, modelSpecs)

			if config.UseRefState(modelSpecs):
				ref_loss, ref_errors = ValidateAllData(validData, quickValidate, modelSpecs, forRefState=True)
                        	print( 'epoch %2d, minibatch %4d/%4d, train loss %s, train error %s, paramL2 %.2f, valid loss %s, valid error %s, ref loss %s, ref error %s' % (epoch, minibatch_index + 1, n_train_batches, str_display(train_loss), str_display(train_errors), param_L2, str_display(valid_loss), str_display(valid_errors), str_display(ref_loss), str_display(ref_errors) ) )
			else:
                        	print( 'epoch %2d, minibatch %4d/%4d, train loss %s, train error %s, paramL2 %.2f, valid loss %s, valid error %s' % (epoch, minibatch_index + 1, n_train_batches, str_display(train_loss), str_display(train_errors), param_L2, str_display(valid_loss), str_display(valid_errors) ) )

                        result = dict()
                        result['valid_loss'] = valid_loss
                        result['valid_errors'] = valid_errors

			if config.UseRefState(modelSpecs):
                        	result['ref_loss'] = ref_loss
                        	result['ref_errors'] = ref_errors

                        results.append(result)

                        if np.mean(valid_loss) < np.mean(chkpoint['best_validation_loss'] ):
                                chkpoint['best_validation_loss'] = valid_loss
                                chkpoint['best_validation_err'] = valid_errors
                                chkpoint['best_epoch'] = epoch
                                chkpoint['bestParamValues'] = [param.get_value () for param in params ]
                                chkpoint['bestOtherParamValues'] = [param.get_value () for param in other_params ]
				chkpoint['train_loss4best_validation_loss'] = train_loss

		## save check point for restart 
		if ( minibatch_index + 1 ) % chkpt_save_frequency  == 0:

                        chkpoint['param_values'] = [ p.get_value(borrow=True) for p in params ]
                        chkpoint['other_param_values'] = [ u.get_value(borrow=True) for u in other_params ]
                        chkpoint['epoch'] = epoch

                        ## save checkpoint for restart
                        if modelSpecs.has_key('checkpointFile'):
                                checkPointFile = modelSpecs['checkpointFile']
                                if checkPointFile is not None:
                                        fh = open(checkPointFile, 'wb')
                                        cPickle.dump(chkpoint, fh, protocol=cPickle.HIGHEST_PROTOCOL)
                                        fh.flush()
                                        fh.close()

		## check to see if we shall skip a portion of one epoch. Read a file SkipOneEpoch.pid to obtain which epoch to skip. 
		## For example, if the file has a float 16.4, then the 17th epoch will be skipped after training 40% of the epoch
		## To skip one whole epoch, please see AdjustEpochs()
		chk_skip_frequency = chkpt_save_frequency
                if ( minibatch_index + 1 ) % chk_skip_frequency ==0 and TrainUtils.SkipOneEpoch(epoch, minibatch_index, n_train_batches):
			print 'Skip the rest of epoch ... ', epoch 
			break

        ## statistics
        avg_train_loss = np.average( t_loss, axis=0  )
        avg_train_errors = np.average( t_errors, axis=0 )
        avg_valid_errors = np.average( [ res['valid_errors'] for res in results ], axis=0 )
        avg_valid_loss = np.average( [ res['valid_loss'] for res in results ], axis=0 )
        print( 'average result at epoch %i: train loss %s, train error %s, valid loss %s, valid error %s' %
                 ( epoch, str_display(avg_train_loss), str_display(avg_train_errors), str_display(avg_valid_loss), str_display(avg_valid_errors) ) )

        ## test on the validation data
	assert validData is not None
        validLoss, validErr, validAcc = ValidateAllData(validData, fullValidate, modelSpecs)

	if config.UseRefState(modelSpecs):
        	refLoss, refErr = ValidateAllData(validData, quickValidate, modelSpecs, forRefState=True)
        	print 'valid loss: ', validLoss, 'valid err: ', validErr, 'ref loss: ', refLoss, 'ref err: ', refErr
	else:
        	print 'valid loss: ', validLoss, 'valid err: ', validErr

	#print validAcc
	print "validAcc: ", [ str_display(vAcc[:, 0]) for vAcc in validAcc ], ' for top seqLen *', modelSpecs['topRatios'], ' long-, medium-, long+medium- and short-range predictions'

        if np.mean(validLoss) < np.mean(chkpoint['best_validation_loss'] ):
                chkpoint['best_validation_loss'] = validLoss
                chkpoint['best_validation_err'] = validErr
                chkpoint['best_validation_acc'] = validAcc
		if config.UseRefState(modelSpecs):
                	chkpoint['best_validation_ref_loss'] = refLoss
                chkpoint['best_epoch'] = epoch
                chkpoint['bestParamValues'] = [param.get_value () for param in params ]
                chkpoint['bestOtherParamValues'] = [param.get_value () for param in other_params ]
		chkpoint['train_loss4best_validation_loss'] = train_loss

	print 'best epoch: ', chkpoint['best_epoch'], 'best valid loss: ', chkpoint['best_validation_loss'], 'train loss: ', chkpoint['train_loss4best_validation_loss']
	print 'end time of epoch ', epoch, ': ', datetime.datetime.now()

## train the model for one stage, which consists of a few epochs with the same learning rate
## if startFromBest is True, then start from the previously best model parameters
## data is a tuple (trainMetaData, validMetaData)
def RunOneStage(epoch_start, epoch_end, data, chkpoint, loss4train, loss4validate, pgrads, pdecay, modelSpecs, lr=np.float32(0.0001), startFromBest=(False, False) ):

	params = modelSpecs['params']
	errors = modelSpecs['errors']
	topAcc = modelSpecs['topAcc']
	algorithm = modelSpecs['algorithm']
	variable4train = modelSpecs['variable4train']
	variable4validate = modelSpecs['variable4validate']
	paramL2 = modelSpecs['paramL2']

        updates, other_params = UpdateAlgorithm(algorithm, params, pgrads, pdecay, lr=np.float32(lr), l2reg=modelSpecs['L2reg'] )

        if startFromBest[0] and chkpoint.has_key('bestParamValues') and Compatible(params, chkpoint['bestParamValues']):
        	[ u.set_value(v) for u, v in zip(params, chkpoint['bestParamValues']) ]

        if startFromBest[1] and chkpoint.has_key('bestOtherParamValues') and Compatible(other_params, chkpoint['bestOtherParamValues']):
                [ u.set_value(v) for u, v in zip(other_params, chkpoint['bestOtherParamValues']) ]

        train = theano.function(variable4train, [loss4train, errors, paramL2], updates=updates, on_unused_input='warn')
        quickValidate = theano.function(variable4validate, [loss4validate, errors], on_unused_input='warn')
        fullValidate = theano.function(variable4validate, [loss4validate, errors, topAcc], on_unused_input='warn')

        epoch = epoch_start
        while (epoch < epoch_end):
        	epoch += 1
                RunOneEpoch(epoch, data, chkpoint, params, other_params, train, quickValidate, fullValidate, modelSpecs)
                if epoch>=14 and (epoch - chkpoint['best_epoch'] > modelSpecs['patience']):
                	break

		## check to see if we want to run more or fewer epochs 
		change = TrainUtils.AdjustEpochs()
		if change != 0:
			epoch = max(0, epoch - change )
			chkpoint['best_epoch'] = max(0, chkpoint['best_epoch'] - change )
			if change >0:
				print 'Running ', change, ' more epochs...'
			elif change <0:
				print 'Running ', -change, ' fewer epochs...'

	## make a copy of check point file at the end of each stage
	if modelSpecs.has_key('checkpointFile') and (modelSpecs['checkpointFile'] is not None):
		dst = modelSpecs['checkpointFile'] + '-epoch' + str(epoch) 
		from shutil import copyfile
		copyfile(modelSpecs['checkpointFile'], dst)

        return epoch

def PrepareModel(modelSpecs=None):

	distancePredictor, x, y, xmask, ymask, xem, labelList, weightList, box, trainByRefLoss = BuildModel(modelSpecs)

        variables = [x, y, xmask, ymask]
	if xem is not None:
		variables.append(xem)

	assert (labelList is not None and len(labelList)>0 )
	variables.extend(labelList)
	variables.extend(weightList)

	if box is not None:
		variables.append(box)

	variable4validate = variables

	if trainByRefLoss is not None:
		variables.append(trainByRefLoss)
	variable4train = variables

	print '#all variables (including response) used in the model: ', len(variables)
	print 'all variables (including response) used in the model: ', variables

    	params = distancePredictor.params
	params4var = distancePredictor.params4var
        params4mean = list( set(params) - set(params4var) )

    	param_shapes = [ p.get_value().shape for p in params ]
    	param_sizes = map( np.prod, param_shapes )
    	modelSpecs['numParams'] = sum(param_sizes)

    	param4mean_shapes = [ p.get_value().shape for p in params4mean ]
    	param4mean_sizes = map( np.prod, param4mean_shapes )
    	modelSpecs['numParams4mean'] = sum(param4mean_sizes)

    	param4var_shapes = [ p.get_value().shape for p in params4var ]
    	param4var_sizes = map( np.prod, param4var_shapes )
    	modelSpecs['numParams4var'] = sum(param4var_sizes)

    	print 'The model has ', modelSpecs['numParams'], ' parameters, including ', modelSpecs['numParams4mean'], ' for mean and ', modelSpecs['numParams4var'], ' for variance.'

	regularizer = 0
	"""
        paramL1 = distancePredictor.paramL1
        if modelSpecs.has_key('L1reg') and modelSpecs['L1reg']>0:
                regularizer += (modelSpecs['L1reg'] * paramL1)
	"""

	## we use only L2 regularization
        paramL2 = distancePredictor.paramL2
        if modelSpecs.has_key('L2reg') and modelSpecs['L2reg']>0:
                regularizer += (modelSpecs['L2reg'] * paramL2)

    	if weightList is not None and len(weightList)>0:
		errors = distancePredictor.errors(labelList, weightList)
    	else:
        	errors = distancePredictor.errors(labelList)
    	topAcc = distancePredictor.TopAccuracyByRange(labelList)

	modelSpecs['variable4train'] = variable4train
	modelSpecs['variable4validate'] = variable4validate
	modelSpecs['params'] = params
	modelSpecs['params4mean'] = params4mean
	modelSpecs['params4var'] = params4var
	modelSpecs['paramL2'] = paramL2
	modelSpecs['regularizer'] = regularizer
	modelSpecs['topAcc'] = topAcc
	modelSpecs['errors'] = errors
	modelSpecs['labelList'] = labelList
	modelSpecs['weightList'] = weightList
	modelSpecs['trainByRefLoss'] = trainByRefLoss

	return distancePredictor, variable4train, variable4validate, params, params4mean, params4var, paramL2, regularizer, topAcc, errors, labelList, weightList, trainByRefLoss

## clean up after training
def Cleanup():

	global validDataLoader, stopValidDataLoader, validData
	#print 'In Cleanup: ', stopValidDataLoader
	if validDataLoader.is_alive():	
		#print 'notify the validation data loader to stop'
		stopValidDataLoader.set()

	global stopTrainDataLoader, trainDataLoaders
	#print 'In Cleanup: ', stopTrainDataLoader
	for trainDataLoader in trainDataLoaders:
		if trainDataLoader.is_alive():
			#print 'notify the train data loader to stop'
			stopTrainDataLoader.set()

	if validData is not None:
		print 'free validation data...'
		FreeSharedArrayInBatches(validData)
		del validData[:]
		validData = None

	time.sleep(10)

	##free shared memory occupied by those data still in trainSharedQ
	print 'free shared array in trainSharedQ, which currently has about ', trainSharedQ.qsize(), ' objects'
	#print 'trainSharedQ.empty(): ', trainSharedQ.empty()
	#while not trainSharedQ.empty():
	while trainSharedQ.qsize()>0:
		d = trainSharedQ.get()
		FreeSharedArrayInOneBatch(d)
	print 'Please double check files (starting with RaptorX-' + str(os.getpid()), ') at /dev/shm/ to make sure all shared memory has been freed'

	if validDataLoader.is_alive():
		print 'stop the validation data loader...'
		validDataLoader.terminate()
	validDataLoader.join()

def TrainModel(modelSpecs, trainValidData=None, predDataFile=None):
    	if (not trainValidData):
        	print 'Please provide train and validation data for model training'
		exit(1)

    	if modelSpecs is None:
        	print 'Please provide a model specification for training'
        	exit(1)

	distancePredictor, variable4train, variable4validate, params, params4mean, params4var, paramL2, regularizer, topAcc, errors, labelList, weightList, trainByRefLoss = PrepareModel(modelSpecs)

	chkpoint, restart = InitializeChkpoint(params, modelSpecs)

	assert ( len(modelSpecs['numEpochs'] ) > 0 )
	numEpochs4stages = np.cumsum(modelSpecs['numEpochs'] )
	## train parameters not related to variance and correlation
        epoch = chkpoint['epoch']

	if epoch < numEpochs4stages[-1]:

                if weightList is not None and len(weightList)>0:
                        loss4train = distancePredictor.loss(labelList, useMeanOnly=True, weightList=weightList, trainByRefLoss=trainByRefLoss)
                        loss4validate = distancePredictor.loss(labelList, useMeanOnly=True, weightList=weightList)
                else:
                        loss4train = distancePredictor.loss(labelList, useMeanOnly=True, trainByRefLoss=trainByRefLoss)
                        loss4validate = distancePredictor.loss(labelList, useMeanOnly=True)

		"""
		## weightedLoss is only used for cost, i.e., gradient calculation
		if modelSpecs.has_key('ScaleLoss4Cost') and (modelSpecs['ScaleLoss4Cost'] is True):
			weightedLoss = ScaleLossByBatchWeight(loss, weightList, modelSpecs)
		else:
			weightedLoss = loss
		"""
		if modelSpecs['algorithm'] in set(['AdamW', 'AdamWAMS']):
                	cost = T.sum( T.mul(loss4train, modelSpecs['w4responses']) ) / np.sum(modelSpecs['w4responses'])
		else:
                	cost = T.sum( T.mul(loss4train, modelSpecs['w4responses']) ) / np.sum(modelSpecs['w4responses']) + regularizer

                params4var_set = set(params4var)
                pgrads = [ T.grad(cost, p, consider_constant=weightList,disconnected_inputs='warn') if p not in params4var_set else T.zeros_like(p) for p in params ]
                pdecay = [ p if p not in params4var_set else T.zeros_like(p) for p in params ]

	for stage, lr, epoch_end in zip(xrange(len(numEpochs4stages) ), modelSpecs['lrs'], numEpochs4stages):
		if epoch >= epoch_end:
			continue

                print 'training for mean using a learning rate ', lr, ' ...'
		startFromBest = (stage>0 and epoch == numEpochs4stages[stage-1] )
                epoch_start = epoch 
                epoch = RunOneStage(epoch_start, epoch_end, trainValidData, chkpoint, loss4train, loss4validate, pgrads, pdecay, modelSpecs, lr=lr, startFromBest=(startFromBest, startFromBest) )
		

        ## train parameters only specific to variance and correlation
        numEpochs4var = modelSpecs['numEpochs4var']
	lrs = modelSpecs['lrs4var']

        if len(params4var) > 0:
		assert ( len(numEpochs4var) > 0 )
		assert ( len(lrs) > 0 )
	
		previousEpochs4Stages = numEpochs4stages	
		numEpochs4stages = np.cumsum( numEpochs4var ) + numEpochs4stages[-1]

		if epoch < numEpochs4stages[-1]:
                	print 'Training the parameters specific to correlation and variance ...'

                	if weightList is not None and len(weightList)>0:
                        	loss4train = distancePredictor.loss(labelList, weightList = weightList, trainByRefLoss=trainByRefLoss)
                        	loss4validate = distancePredictor.loss(labelList, weightList = weightList)
                	else:
                        	loss4train = distancePredictor.loss(labelList)
                        	loss4validate = distancePredictor.loss(labelList)

			"""
			## weightedLoss is only used for cost, i.e., gradient calculation
			if modelSpecs.has_key('ScaleLoss4Cost') and (modelSpecs['ScaleLoss4Cost'] is True):
				weightedLoss = ScaleLossByBatchWeight(loss, weightList, modelSpecs)
			else:
				weightedLoss = loss
			"""

			if modelSpecs['algorithm'] in set(['AdamW', 'AdamWAMS']):
                		cost = T.sum( T.mul(loss4train, modelSpecs['w4responses']) ) / np.sum(modelSpecs['w4responses'])
			else:
                		cost = T.sum( T.mul(loss4train, modelSpecs['w4responses']) ) / np.sum(modelSpecs['w4responses']) + regularizer

                	params4var_set = set(params4var)
                	pgrads = [ T.grad(cost, p, consider_constant=weightList, disconnected_inputs='raise') if p in params4var_set else T.zeros_like(p) for p in params ]
                	pdecay = [ p if p in params4var_set else T.zeros_like(p) for p in params ]

		for stage, lr, epoch_end in zip(xrange(len(lrs)), lrs, numEpochs4stages):
			if epoch >= epoch_end:
				continue

			print 'training for variance using a learning rate ', lr, ' ...'
			startFromBest = ( (stage==0 and epoch == previousEpochs4Stages[-1])  or (stage>0 and epoch == numEpochs4stages[stage-1] ) )
			epoch_start = epoch
                	epoch = RunOneStage(epoch_start, epoch_end, trainValidData, chkpoint, loss4train, loss4validate, pgrads, pdecay, modelSpecs, lr=lr, startFromBest=(startFromBest, startFromBest and (stage>0) ) )

    	resultModel = {}
    	resultModel['dateTrained']=datetime.datetime.now()
    	#resultModel['validLoss'] = validLoss
    	resultModel['validLoss'] = chkpoint['best_validation_loss']
    	#resultModel['validErr'] = validErr
	if chkpoint.has_key('best_validation_err'):
    		resultModel['validErr'] = chkpoint['best_validation_err']

	resultModel['trainLoss'] = chkpoint['train_loss4best_validation_loss']
    	#resultModel['validAcc']= validAcc
	if chkpoint.has_key('best_validation_acc'):
    		resultModel['validAcc']= chkpoint['best_validation_acc']

    	resultModel['paramValues']=chkpoint['bestParamValues']

    	bestParamL2norm = np.sum([(v**2).sum() for v in chkpoint['bestParamValues'] ])
    	resultModel['bestParamL2norm']=bestParamL2norm

    	bestParamL1norm = np.sum([abs(v).sum() for v in chkpoint['bestParamValues'] ])
    	resultModel['bestParamL1norm']=bestParamL1norm

    	print 'best param L1 norm: ', bestParamL1norm, 'L2 norm: ', bestParamL2norm

	Cleanup()

    	#test on prediction data if it is given. Here the prediction data shall be small to save memory and contain ground truth.
     	if modelSpecs['predFile'] is not None:
		predMetaData = DataProcessor.LoadMetaData(modelSpecs['predFile'])
		predDataLocation = DataProcessor.SampleProteinInfo(predMetaData)
                predBatches = DataProcessor.SplitData2Batches(predDataLocation, numDataPoints=624, modelSpecs=modelSpecs)
		print '\nLoading prediction data...'
                print "#predData minibatches:", len(predBatches)

		predData = []
		for batch in predBatches:
			data = DataProcessor.LoadRealData(batch, modelSpecs, returnMode='list')
			FeatureUtils.CheckModelNDataConsistency(modelSpecs, data)
			#input = TrainUtils.PrepareInput4Prediction(data, modelSpecs, floatType=np.float16)
			input = TrainUtils.PrepareInput4Prediction(data, modelSpecs, floatType=theano.config.floatX)
			predData.append(input)

        	if weightList is not None and len(weightList)>0:
                	loss4validate = distancePredictor.loss(labelList, weightList = weightList)
        	else:
                	loss4validate = distancePredictor.loss(labelList)

		fullValidate = theano.function(variable4validate, [loss4validate, errors, topAcc], on_unused_input='warn')
		if config.UseRefState(modelSpecs):
			quickValidate = theano.function(variable4validate, [loss4validate, errors], on_unused_input='warn')

    		## set model parameters for valiation and possibly prediction
    		for param, value in zip(params, chkpoint['bestParamValues']):
       			param.set_value(value)

        	predLoss, predErr, predAcc = ValidateAllData(predData, fullValidate, modelSpecs)
		if config.UseRefState(modelSpecs):
        		refLoss, refErr = ValidateAllData(predData, quickValidate, modelSpecs, forRefState=True)
        		print 'pred loss: ', predLoss, 'pred err: ', predErr, 'ref loss: ', refLoss, 'ref err: ', refErr
		else:
        		print 'pred loss: ', predLoss, 'pred err: ', predErr
        	resultModel['predLoss'] = predLoss
        	resultModel['predErr'] = predErr

        	print "predAcc: ", [ str_display(pAcc[:, 0]) for pAcc in predAcc ], 'for top ', modelSpecs['topRatios']
        	resultModel['predAcc'] = predAcc

		del predData[:]

	## training is done, remove the checkpoint file since it has been copied at the end of each stage
	if modelSpecs.has_key('checkpointFile') and (modelSpecs['checkpointFile'] is not None):
		try:
			os.remove(modelSpecs['checkpointFile'])
		except IOError:
			print 'WARNING: error in deleting the check point file: ', modelSpecs['checkpointFile']
	
	## remove theano variables from modelSpecs
	keys4removal = ['variable4train', 'variable4validate', 'params', 'params4mean', 'params4var', 'paramL2', 'regularizer', 'topAcc', 'errors', 'labelList', 'weightList', 'trainByRefLoss']
	for k in keys4removal:
		if modelSpecs.has_key(k):
			del modelSpecs[k]
	
    	return resultModel

def main(argv):

    	modelSpecs = InitializeModelSpecs()
	modelSpecs = ParseCommandLine.ParseArguments(argv, modelSpecs)

	startTime = datetime.datetime.now()

	trainMetaData = DataProcessor.LoadMetaData(modelSpecs['trainFile'] )
	FeatureUtils.DetermineFeatureDimensionBySampling(trainMetaData, modelSpecs)
	## calculate label distribution and weight at the very beginning
	print 'Calculating label distribution...'
	LabelUtils.CalcLabelDistributionNWeightBySampling(trainMetaData, modelSpecs)

	if config.TrainByRefLoss(modelSpecs) or config.UseRefState(modelSpecs):
		print 'Calculating feature expection by sampling...'
		FeatureUtils.CalcFeatureExpectBySampling(trainMetaData, modelSpecs)

	## trainMetaData is a list of groups. Each group contains a set of related proteins (seq-template alignments) and files for their features
	trainDataLocation = DataProcessor.SampleProteinInfo(trainMetaData)
        trainSeqData = DataProcessor.SplitData2Batches(trainDataLocation, numDataPoints=modelSpecs['minibatchSize'], modelSpecs=modelSpecs)
	print 'approximate #batches for train data: ', len(trainSeqData)

	#global trainSharedQ, stopTrainDataLoader, trainDataLoaders, trainSharedLabelPool, trainSharedLabelWeightPool
	global trainSharedQ, stopTrainDataLoader, trainDataLoaders
	trainSharedQ = multiprocessing.Queue(config.QSize(modelSpecs) )
	stopTrainDataLoader = multiprocessing.Event()
	#trainSharedLabelPool = multiprocessing.Manager().dict()
	#trainSharedLabelWeightPool = multiprocessing.Manager().dict()
	#print stopTrainDataLoader

	numTrainDataLoaders = config.NumTrainDataLoaders(modelSpecs)
	metaDatas = DataProcessor.SplitMetaData(trainMetaData, numTrainDataLoaders)

	trainDataLoaders = []
	for i, metaData in zip(xrange(numTrainDataLoaders), metaDatas):
        	#trainDataLoader = multiprocessing.Process(name='TrainDataLoader ' + str(i) + ' for ' + str(os.getpid()), target=TrainUtils.TrainDataLoader, args=(trainSharedQ, metaData, modelSpecs, True, True))
        	trainDataLoader = multiprocessing.Process(name='TrainDataLoader ' + str(i) + ' for ' + str(os.getpid()), target=TrainUtils.TrainDataLoader2, args=(trainSharedQ, stopTrainDataLoader, metaData, modelSpecs, True, True))
        	#trainDataLoader = multiprocessing.Process(name='TrainDataLoader ' + str(i) + ' for ' + str(os.getpid()), target=TrainUtils.TrainDataLoader3, args=(trainSharedQ, trainSharedLabelPool, trainSharedLabelWeightPool, stopTrainDataLoader, metaData, modelSpecs, True, True))
		trainDataLoader.daemon=True
		trainDataLoaders.append(trainDataLoader)

        print 'start the train data loaders...'
	for trainDataLoader in trainDataLoaders:
        	trainDataLoader.start()

	validMetaData = DataProcessor.LoadMetaData(modelSpecs['validFile'])
	validDataLocation = DataProcessor.SampleProteinInfo(validMetaData)

	## split data into batches, but do not load the real data from disk
	#validSeqData = DataProcessor.SplitData2Batches(validDataLocation, numDataPoints=modelSpecs['minibatchSize'], modelSpecs=modelSpecs)
	validSeqData = DataProcessor.SplitData2Batches(validDataLocation, numDataPoints=500*500, modelSpecs=modelSpecs)
	print '#batches for validation data: ', len(validSeqData)

	global validSharedQ, validDataLoader, stopValidDataLoader
	validSharedQ = multiprocessing.Queue(len(validSeqData) )
	stopValidDataLoader = multiprocessing.Event()
	#print stopValidDataLoader
	## shared memory is a limited resource, so avoid using it as much as possible
	## here we do not use shared array for validation data since we only need to load it once
        #validDataLoader = multiprocessing.Process(name='ValidDataLoader for '+str(os.getpid()), target=TrainUtils.ValidDataLoader, args=(validSharedQ, validSeqData, modelSpecs, True, False))
        validDataLoader = multiprocessing.Process(name='ValidDataLoader for '+str(os.getpid()), target=TrainUtils.ValidDataLoader2, args=(validSharedQ, stopValidDataLoader, validSeqData, modelSpecs, True, False))
        print 'start the validation data loader...'
        validDataLoader.start()

	"""
	if modelSpecs.has_key('ScaleLoss4Cost') and (modelSpecs['ScaleLoss4Cost'] is True):
		##calculate the average weight per minibatch
		maxDeviation = DataProcessor.CalcAvgWeightPerBatch(trainSeqDataset, modelSpecs)
		print 'maxWeightDeviation=', maxDeviation
	"""

	beforeTrainTime = datetime.datetime.now()
	print 'time spent before training :', beforeTrainTime - startTime

        result = TrainModel(modelSpecs=modelSpecs, trainValidData=(trainSeqData, validSeqData))


        ##merge ModelSpecs and result
        resultModel = modelSpecs.copy()
        resultModel.update(result)

	modelFile = TrainUtils.GenerateModelFileName(resultModel)
	print 'Writing the resultant model to ', modelFile
        cPickle.dump(resultModel, file(modelFile, 'wb'), cPickle.HIGHEST_PROTOCOL)

	afterTrainTime = datetime.datetime.now()
	print 'time spent on training:', afterTrainTime - beforeTrainTime

	## clean up again
	print 'Cleaning up again...'
	Cleanup()

if __name__ == "__main__":

	print 'Start time: ',  str(datetime.datetime.now() )

   	recursionlimit = sys.getrecursionlimit()
   	#print 'recursionlimit = ', recursionlimit
   	sys.setrecursionlimit(3*recursionlimit)

	## generate a random seed
	a = list( str(os.getpid())  +  os.urandom(8) )
	random.shuffle(a)
	seed = ''.join(a)
	#print 'setting random seed: ', seed
	random.seed(a=seed)

	validSharedQ = None
	validData = None 
	validDataLoader = None
	stopValidDataLoader = None

	trainSharedQ = None
	#trainSharedLabelPool = None
	#trainSharedWeightPool = None
	#trainData = None
	trainDataLoaders = []
	stopTrainDataLoader = None

   	main(sys.argv[1:])
