import numpy as np
import cPickle
import os
import sys
import time
import datetime
import random
import gc
import copy

import theano.tensor as T
import theano

import config
from config import Response2LabelName, Response2LabelType, ParseResponse, GetResponseProbDims, GetResponseValueDims
import ParseCommandLine as ParseCommandLine

from Model4DistancePrediction import BuildModel

from Optimizers import SGDM, SGDM2, Nesterov
from Adams import Adam, AdamW, AdamWAMS, AMSGrad

import DistanceUtils
import FeatureUtils
import LabelUtils
import DataProcessor

from utilsNoT import SampleBoundingBox, Compatible, str_display, IsNumber

from Initialize import InitializeModelSpecs


## pdecay is the decay of params. It is only used for AdamW and AdamWAMS
## if pdecay is None, then in AdamW and AdamWAMS, pdecay is set to params
def UpdateAlgorithm(alg, params, gparams, pdecay, lr=None, l2reg=None):
    	## currently the code mainly work for Adam. Other algorithms are not fully tested.

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

## Initialize check point. We do not set the values of other_params since when this function is called, other_params may not be defined
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
                print 'could not find the check point file: ', checkPointFile
                print 'training the model from scratch ...'

        return chkpoint, restart

## calculate loss, error and accuracy for a set of data
## validate is the function that executes the validation procedure
## validate has two options: fullvalidate or quickvalidate
def ValidateAllData(SeqDataset, validate, modelSpecs, forRefState=False):
	accs = []
	losses = []
	errs = []
        numSamples = []

        if modelSpecs['UseSampleWeight']:
                w4losses = []
                w4errors = []
        else:
                w4losses = None
                w4errors = None

        for minibatch in SeqDataset:
		onebatch, _= DataProcessor.AssembleOneBatch(minibatch, modelSpecs, forRefState=forRefState)

		## add a bounding box, which shall be equivalent to the shape of pairwise features
		x2d = onebatch[1]
		boundingbox = np.array([0, 0, x2d.shape[1], x2d.shape[2] ]).astype(np.int32)

	       	onebatch_res = validate( *(onebatch + [boundingbox]) )
		los = onebatch_res[0]
		err = onebatch_res[1]
		losses.append(los)
	    	errs.append(err)

		if len(onebatch_res) > 2:
			acc = onebatch_res[2]
	    		accs.append( acc )
			##numSamples is the number of proteins
	    		numSamples.append(onebatch[0].shape[0])

                if modelSpecs['UseSampleWeight']:
                        weights = onebatch[ len(onebatch) - len(modelSpecs['responses']) : ]
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

## this function is used to skip a portion of one epoch. It reads a file SkipOneEpoch.pid to obtain which epoch to skip where pid is the process id of the training program. 
## For example, if the file has a float 16.4, then the 17th epoch will be skipped after training 40% of the 16th epoch
## To skip one whole epoch, please see AdjustEpochs()
def SkipOneEpoch(epoch, minibatch_index, n_train_batches):

	file4skip = 'SkipOneEpoch.' + str(os.getpid() )
	if not os.path.isfile(file4skip):
		return False

	first_line='0'
	try:
		with open(file4skip, 'r') as f:
    			first_line = f.readline()
	except IOError:
		print 'WARNING: could not read file ', file4skip
		return False

	try:
		epoch4stop = np.float32( first_line.strip() )
	except ValueError:
		print 'WARNING: value in file is not a valid float: ', file4skip
		return False

	if epoch4stop <= (epoch - 1) or epoch4stop >= epoch:
		## epoch4stop is not for this epoch
		return False

	n_cutoff_batches = (epoch4stop + 1 - epoch) * n_train_batches

	if n_cutoff_batches > minibatch_index:
		return False

	## now we shall skip the rest of this epoch
	try:
		os.remove(file4skip)
	except IOError:
		print 'WARNING: could not remove file ', file4skip

	return True

## this function can be used to add or remove a few epochs in training
def AdjustEpochs():
	## this function reads a number from a file named AdjustEpoch.pid where pid is the process id returned by os.getpid()
	## When the number is positive, then run a few more epochs. When it is negative, then reduce a few epochs.

	file4EpochChange = 'AdjustEpoch.' + str(os.getpid() )
	if not os.path.isfile(file4EpochChange):
		return 0

	try:
		with open(file4EpochChange, 'r') as f:
    			first_line = f.readline()
		os.remove(file4EpochChange)

	except IOError:
		print 'ERROR: could not read or write file ', file4EpochChange
		return 0

	try:
		change = np.int32(first_line.strip())
	except ValueError:
		print 'WARNING: value in file is not a valid integer: ', file4EpochChange
		return 0

	return change

def TrainByOneBatch(minibatch, train, modelSpecs, forRefState=False):

	## we have to crop input for a very large protein to avoid crash due to limited GPU memory
	onebatch, names4onebatch = DataProcessor.AssembleOneBatch(minibatch, modelSpecs, forRefState=forRefState)
        x1d, x2d, x1dmask, x2dmask = onebatch[0:4]
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


        ##crop the remaining input including the ground truth and weight matrices
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

## train the model for one epoch. That is, scanning all the training data once
def RunOneEpoch(epoch, chkpoint, params, other_params, trainSeqData, validSeqData, train, quickValidate, fullValidate, modelSpecs):
	print 'start time of epoch ', epoch, ': ', datetime.datetime.now()

        random.shuffle(trainSeqData)
	n_train_batches = len(trainSeqData)

	if epoch > 16:
		validation_frequency = np.int32(0.8 * modelSpecs['validation_frequency'])
	elif epoch > 13:
		validation_frequency = np.int32(1.0 * modelSpecs['validation_frequency'])
	elif epoch > 11:
		validation_frequency = np.int32(1.3 * modelSpecs['validation_frequency'])
	elif epoch > 9:
		validation_frequency = np.int32(1.8 * modelSpecs['validation_frequency'])
	elif epoch > 5:
		validation_frequency = np.int32(2.2 * modelSpecs['validation_frequency'])
	else:
		validation_frequency = np.int32(2.5 * modelSpecs['validation_frequency'])

	chkpt_save_frequency = max(1, min(validation_frequency, np.int32( len(trainSeqData) / 30 ) ) )

        t_loss = []
        t_errors = []
        results = []

        for minibatch_index, minibatch in zip(xrange( n_train_batches ), trainSeqData):
		train_loss, train_errors, param_L2 = TrainByOneBatch(minibatch, train, modelSpecs)
                t_loss.append(train_loss)
                t_errors.append(train_errors)

		if config.TrainByRefLoss(modelSpecs):
			## add code here to train by a reference input
			_, _, param_L2 = TrainByOneBatch(minibatch, train, modelSpecs, forRefState=True)

		## validate the model	
                if ( minibatch_index + 1 ) % validation_frequency  == 0:

			valid_loss, valid_errors = ValidateAllData(validSeqData, quickValidate, modelSpecs)

			if config.UseRefState(modelSpecs):
				ref_loss, ref_errors = ValidateAllData(validSeqData, quickValidate, modelSpecs, forRefState=True)
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
                if ( minibatch_index + 1 ) % chk_skip_frequency ==0 and SkipOneEpoch(epoch, minibatch_index, n_train_batches):
			print 'Skip the rest of epoch ', epoch 
			break
				

        ## statistics
        avg_train_loss = np.average( t_loss, axis=0  )
        avg_train_errors = np.average( t_errors, axis=0 )
        avg_valid_errors = np.average( [ res['valid_errors'] for res in results ], axis=0 )
        avg_valid_loss = np.average( [ res['valid_loss'] for res in results ], axis=0 )
        print( 'average result at epoch %i: train loss %s, train error %s, valid loss %s, valid error %s' %
                 ( epoch, str_display(avg_train_loss), str_display(avg_train_errors), str_display(avg_valid_loss), str_display(avg_valid_errors) ) )

        ## test on the validation data
        validLoss, validErr, validAcc = ValidateAllData(validSeqData, fullValidate, modelSpecs)

	if config.UseRefState(modelSpecs):
        	refLoss, refErr = ValidateAllData(validSeqData, quickValidate, modelSpecs, forRefState=True)
        	print 'valid loss: ', validLoss, 'valid err: ', validErr, 'ref loss: ', refLoss, 'ref err: ', refErr
	else:
        	print 'valid loss: ', validLoss, 'valid err: ', validErr

	#print validAcc
	print "validAcc: ", [ str_display(vAcc[:, 0]) for vAcc in validAcc ], ' for top ', modelSpecs['topRatios']

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


##calculate the weight for each minibatch where weightList is the weight matrix and base is the normalization
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

def TrainModel(modelSpecs=None, trainSeqData=None, validSeqData=None, predSeqData=None):
    	if modelSpecs is None:
        	print 'Please provide a model specification for training'
        	exit(1)

    	if (not trainSeqData) or (not validSeqData):
        	print 'Please provide train and validation data for model training'
		exit(1)

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
    	print 'The model has ', modelSpecs['numParams'], ' parameters.'

    	param4mean_shapes = [ p.get_value().shape for p in params4mean ]
    	param4mean_sizes = map( np.prod, param4mean_shapes )
    	modelSpecs['numParams4mean'] = sum(param4mean_sizes)
    	print 'The model has ', modelSpecs['numParams4mean'], ' parameters for mean.'

    	param4var_shapes = [ p.get_value().shape for p in params4var ]
    	param4var_sizes = map( np.prod, param4var_shapes )
    	modelSpecs['numParams4var'] = sum(param4var_sizes)
    	print 'The model has ', modelSpecs['numParams4var'], ' parameters for variance.'

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

        modelSpecs['validation_frequency'] = min( max( len(validSeqData), len(trainSeqData)/12 ), len(trainSeqData)/6)

    	if weightList is not None and len(weightList)>0:
		errors = distancePredictor.errors(labelList, weightList)
    	else:
        	errors = distancePredictor.errors(labelList)

    	topAcc = distancePredictor.TopAccuracyByRange(labelList)

	chkpoint, restart = InitializeChkpoint(params, modelSpecs)

	## train the model for one stage, which consists of a few epochs with the same learning rate
	## if startFromBest is True, then start from the previously best model parameters
	def RunOneStage(epoch_start, epoch_end, loss4train, loss4validate, gparams, lr, pdecay, algorithm='Adam', startFromBest=(False, False) ):

                #updates, other_params = UpdateAlgorithm(modelSpecs['algorithm'], params, gparams, param_shapes, np.float32(lr) )
                updates, other_params = UpdateAlgorithm(algorithm, params, gparams, pdecay, lr=np.float32(lr), l2reg=modelSpecs['L2reg'] )

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

                        RunOneEpoch(epoch, chkpoint, params, other_params, trainSeqData, validSeqData, train, quickValidate, fullValidate, modelSpecs)
                        if epoch>=14 and (epoch - chkpoint['best_epoch'] > modelSpecs['patience']):
                                break

			## check to see if we want to run more or fewer epochs 
			change = AdjustEpochs()
			if change != 0:
				epoch = max(0, epoch - change )
				chkpoint['best_epoch'] = max(0, chkpoint['best_epoch'] - change )
				if change >0:
					print 'Running ', change, ' extra epochs...'
				elif change <0:
					print 'Running ', -change, ' fewer epochs...'

		## make a copy of check point file at the end of each stage
		if modelSpecs.has_key('checkpointFile') and (modelSpecs['checkpointFile'] is not None):
			dst = modelSpecs['checkpointFile'] + '-epoch' + str(epoch) 
			from shutil import copyfile
			copyfile(modelSpecs['checkpointFile'], dst)

                return epoch


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
                gparams = [ T.grad(cost, p, consider_constant=weightList,disconnected_inputs='warn') if p not in params4var_set else T.zeros_like(p) for p in params ]
                pdecay = [ p if p not in params4var_set else T.zeros_like(p) for p in params ]

	for stage, lr, epoch_end in zip(xrange(len(numEpochs4stages) ), modelSpecs['lrs'], numEpochs4stages):
		if epoch >= epoch_end:
			continue

                print 'training for mean using a learning rate ', lr, ' ...'
		startFromBest = (stage>0 and epoch == numEpochs4stages[stage-1] )
                epoch_start = epoch 
                epoch = RunOneStage(epoch_start, epoch_end, loss4train, loss4validate, gparams, lr, pdecay=pdecay, algorithm=modelSpecs['algorithm'], startFromBest=(startFromBest, startFromBest) )
		

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
                	gparams = [ T.grad(cost, p, consider_constant=weightList, disconnected_inputs='raise') if p in params4var_set else T.zeros_like(p) for p in params ]
                	pdecay = [ p if p in params4var_set else T.zeros_like(p) for p in params ]

		for stage, lr, epoch_end in zip(xrange(len(lrs)), lrs, numEpochs4stages):
			if epoch >= epoch_end:
				continue

			print 'training for variance using a learning rate ', lr, ' ...'
			startFromBest = ( (stage==0 and epoch == previousEpochs4Stages[-1])  or (stage>0 and epoch == numEpochs4stages[stage-1] ) )
			epoch_start = epoch
                	epoch = RunOneStage(epoch_start, epoch_end, loss4train, loss4validate, gparams, lr, pdecay=pdecay, algorithm=modelSpecs['algorithm'], startFromBest=(startFromBest, startFromBest and (stage>0) ) )

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

    	#test on prediction data if it is given. Here the prediction data shall contain ground truth.
     	if predSeqData is not None:
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

        	predLoss, predErr, predAcc = ValidateAllData(predSeqData, fullValidate, modelSpecs)
		if config.UseRefState(modelSpecs):
        		refLoss, refErr = ValidateAllData(predSeqData, quickValidate, modelSpecs, forRefState=True)
        		print 'pred loss: ', predLoss, 'pred err: ', predErr, 'ref loss: ', refLoss, 'ref err: ', refErr
		else:
        		print 'pred loss: ', predLoss, 'pred err: ', predErr
        	resultModel['predLoss'] = predLoss
        	resultModel['predErr'] = predErr

        	print "predAcc: ", [ str_display(pAcc[:, 0]) for pAcc in predAcc ], 'for top ', modelSpecs['topRatios']
        	resultModel['predAcc'] = predAcc

	## training is done, remove the checkpoint file since it has been copied at the end of each stage
	if modelSpecs.has_key('checkpointFile') and (modelSpecs['checkpointFile'] is not None):
		try:
			os.remove(modelSpecs['checkpointFile'])
		except IOError:
			print 'WARNING: error in deleting the check point file: ', modelSpecs['checkpointFile']
		

    	return resultModel



##automatically generating a model file name
def GenerateModelFileName(resultModel):

	prefix = ''
	verStr ='Model' +  str(int(resultModel['Version']*10)) + '-'
        if resultModel.has_key('UseTemplate') and resultModel['UseTemplate']:
                prefix += ('TPL' + verStr)
        else:
                prefix += ('Seq' + verStr)

        prefix += resultModel['network'] + '4'
        prefix += resultModel['responseStr'].replace(';', '-').replace(':','_') 

	if config.UseRawCCM(resultModel):
		prefix += '-CCM'
	if config.UseCCMZ(resultModel):
		prefix += '-CCMZ'
	if resultModel.has_key('NoWeight4Range') and resultModel['NoWeight4Range']:
		prefix +='-NoWR'
	if resultModel.has_key('NoWeight4Label') and resultModel['NoWeight4Label']:
		prefix +='-NoWL'

	arch = ''
	#arch +=  'L1D' + str( (sum(resultModel['conv1d_repeats']) + len(resultModel['conv1d_hiddens']) )*2 -1 )  
	arch += 'L2D' + str( (sum(resultModel['conv2d_repeats']) + len(resultModel['conv2d_hiddens']) )*2 -1 )  
	#arch += 'Log' + str( (sum(resultModel['logreg_hiddens']) + len(resultModel['logreg_hiddens']) )      )  

	"""
	if resultModel['network'].startswith('DilatedResNet'):
		arch += 'W1D' + str(resultModel['conv1d_hwsz']) + 'W2D' 
		arch += '-'.join(map(str, resultModel['conv2d_hwszs']) )
		arch += 'Dilation' + '-'.join(map(str, resultModel['conv2d_dilations']) )
	else:
		arch += 'W1D' + str(resultModel['halfWinSize_seq']) + 'W2D' + str(resultModel['halfWinSize_matrix'])
	"""
	#arch += 'I1D' + str(resultModel['n_in_seq']) + 'I2D' + str(resultModel['n_in_matrix'])

	MID = str(datetime.datetime.today()).split()[0].replace('-', '') + '-' + resultModel['ModelID'] + str(os.getpid())
	#bias = resultModel['LRbias'] + 'LRbias'
	##epoch = 'E' + str(resultModel['numEpochs'])
	##epoch = 'E' 
	suffix = '.pkl'

	datastr = os.path.basename(resultModel['dataset'][0]).split('.')[0:-2]
	datastr = ''.join(datastr)

	#components = [ prefix, arch+bias+epoch, datastr, pid, resultModel['algorithm'] ]
	components = [ prefix, arch, datastr, MID, resultModel['algorithm'] ]
	filename = os.path.join('LocalModels/', '-'.join(components) + suffix )

	if not os.path.exists('LocalModels/'):
    		os.makedirs('LocalModels/')

	return filename

## each file shall end with .txt or .list
def ParseListFile(listFiles):

	names = []
	for f in listFiles:
		if not os.path.isfile(f):
			print 'ERROR: cannot find the file: ', f
			exit(1)

		if not ( f.endswith('.txt') or f.endswith('.list') ):
			print 'ERROR: the expected list file shall only end with .txt or .list: ', f
			exit(1)

		fh = open(f, 'r')
		c = [ line.strip() for line in list(fh) ]
		fh.close()
		names.extend(c)

	return set(names)

## ratio is represented as a string
## if ratio is an integer, it shall be interpreted as the number of samples, otherwise the ratio of samples
def SampleProteinNames(whole, ratio, exclude=None):
	if ratio.isdigit():
		numSamples = np.int32(ratio)
	else:
		numSamples = np.int32( len(whole) * np.float32(ratio) )

	if exclude is not None:
		names = list( set(whole) - set(exclude) )
	else:
		names = copy.deepcopy( list(whole) )

	if numSamples > (1.1*len(names) ):
		print 'ERROR: the number of proteins to be sampled greatly exceeds the total number of proteins available: ', numSamples
		print 'The total number of proteins: ', len(names)
		exit(1)
	elif numSamples > len(names):
		print 'WARNING: the number of proteins to be sampled exceeds the total number of proteins available: ', numSamples
		print 'The total number of proteins: ', len(names)
		print 'All the available proteins will be used'
		numSamples = len(names)

	if numSamples < 1:
		print 'ERROR: the number of proteins to be sampled is 0 '
		exit(1)

	random.shuffle(names)

	return set( names[:numSamples] )

	

def main(argv):

    	modelSpecs = InitializeModelSpecs()
	modelSpecs = ParseCommandLine.ParseArguments(argv, modelSpecs)

	startTime = datetime.datetime.now()
	
	## load the datasets. Data is a list of proteins and each protein is represented as a dict()
	Data = DataProcessor.LoadDistanceFeatures(modelSpecs['dataset'], modelSpecs=modelSpecs )
        print '#Data: ', len(Data)
	allProteins = [ d['name'] for d in Data ]

	## Each protein in trainData contains three or four components: seqFeatures, matrixFeatures, embedFeatures and label matrix
        ## embedFeatures is derived from seqFeatures. Users only provide seqFeatures, matrixFeatures, and distance matrix.
        modelSpecs['n_in_seq'] = Data[0]['seqFeatures'].shape[1]
        modelSpecs['n_in_matrix'] = Data[0]['matrixFeatures'].shape[2] + Data[0]['matrixFeatures_nomean'].shape[2]
        if Data[0].has_key('embedFeatures'):
                modelSpecs['n_in_embed'] = Data[0]['embedFeatures'].shape[1]

	print 'Preparing training data...'
	## extract trainData from Data
	if len(modelSpecs['trainFile'])==1 and IsNumber(modelSpecs['trainFile'][0]):
		trainProteinSet = SampleProteinNames(allProteins, modelSpecs['trainFile'][0])
	else:
		trainProteinSet = ParseListFile(modelSpecs['trainFile'])
	trainData = [ d for d in Data if d['name'] in trainProteinSet ]
	print '#trainData: ', len(trainData)
	modelSpecs['trainProteins'] = trainProteinSet
	modelSpecs['numOfTrainProteins']= len(trainData)


	## the results are saved to modelSpecs
	LabelUtils.CalcLabelDistributionAndWeight(trainData, modelSpecs)

	##the results are also saved to modelSpecs
	FeatureUtils.CalcExpectedValueOfFeatures(trainData, modelSpecs)

        print 'Preparing batch data for training...'
        groupSize = modelSpecs['minibatchSize']
        trainSeqDataset = DataProcessor.SplitData2Batches(data=trainData, numDataPoints=groupSize, modelSpecs=modelSpecs)
        print "#trainData minibatches:", len(trainSeqDataset)

	"""
	if modelSpecs.has_key('ScaleLoss4Cost') and (modelSpecs['ScaleLoss4Cost'] is True):
		##calculate the average weight per minibatch
		maxDeviation = DataProcessor.CalcAvgWeightPerBatch(trainSeqDataset, modelSpecs)
		print 'maxWeightDeviation=', maxDeviation
	"""

	print 'Preparing validation data ...'

	if len(modelSpecs['validFile'])==1 and IsNumber(modelSpecs['validFile'][0]):
		validProteinSet = SampleProteinNames(allProteins, modelSpecs['validFile'][0], exclude=trainProteinSet)
	else:
		validProteinSet = ParseListFile(modelSpecs['validFile'])

	overlapSet = validProteinSet.intersection(trainProteinSet)
	if len(overlapSet) > min(10, 0.1*len(validProteinSet) ):
		print 'WARNING: maybe too much overlap between the validation and training sets: ', len(overlapSet)
		
	validData = [ d for d in Data if d['name'] in validProteinSet ]
	modelSpecs['validProteins'] = validProteinSet
	print '#validData: ', len(validData)
	modelSpecs['numOfValidProteins']= len(validData)
        validSeqDataset = DataProcessor.SplitData2Batches(data=validData, numDataPoints=groupSize, modelSpecs=modelSpecs)
        print "#validData minibatches:", len(validSeqDataset)


        predSeqDataset = None
	if modelSpecs['predFile'] is not None:
		if len(modelSpecs['predFile'])==1 and IsNumber(modelSpecs['predFile'][0]):
			predProteinSet = SampleProteinNames(allProteins, modelSpecs['predFile'][0], exclude=trainProteinSet.union(validProteinSet) )
			predData = [ d for d in Data if d['name'] in predProteinSet ]
		else:
			## divide all pred files into two groups, one is for list files and the other for PKL files
			listFiles = [ pfile for pfile in modelSpecs['predFile'] if pfile.endswith('.txt') or pfile.endswith('.list') ]
			pklFiles = [ pfile for pfile in modelSpecs['predFile'] if pfile.endswith('.pkl') ]
            		predDataPKL = DataProcessor.LoadDistanceFeatures(pklFiles, modelSpecs=modelSpecs, forTrainValidation=False )
			predProteinSet = ParseListFile(listFiles)
			predDataTXT = [ d for d in Data if d['name'] in predProteinSet ]
			predData = predDataPKL + predDataTXT

	    	print '#predData: ', len(predData)
	 	predSeqDataset = DataProcessor.SplitData2Batches(data=predData, numDataPoints=624, modelSpecs=modelSpecs)
                print "#predData minibatches:", len(predSeqDataset)


	beforeTrainTime = datetime.datetime.now()
	print 'time spent on loading data:', beforeTrainTime - startTime

        result = TrainModel(modelSpecs=modelSpecs, trainSeqData=trainSeqDataset, validSeqData=validSeqDataset, predSeqData=predSeqDataset)

        ##merge ModelSpecs and result
        resultModel = modelSpecs.copy()
        resultModel.update(result)

	modelFile = GenerateModelFileName(resultModel)
	print 'Writing the resultant model to ', modelFile
        cPickle.dump(resultModel, file(modelFile, 'wb'), cPickle.HIGHEST_PROTOCOL)

	afterTrainTime = datetime.datetime.now()
	print 'time spent on training:', afterTrainTime - beforeTrainTime

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

   	main(sys.argv[1:])
