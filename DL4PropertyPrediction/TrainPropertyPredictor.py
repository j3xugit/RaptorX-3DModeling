import numpy as np
import cPickle
import theano.tensor as T
import theano
import os
import os.path
import sys
import time
import datetime
import random
import gzip

from Adams import Adam, AMSGrad, AdamW, AdamWAMS
from Optimizers import SGDM2, Nesterov
from utils import Compatible

import config
from config import Response2LabelType, Response2LabelName

from Initialize import InitializeModelSpecs
import ParseCommandLine

import DataProcessor
from Model4PropertyPrediction import BuildModel


def str_display(ls):
	if not isinstance(ls, (list, tuple, np.ndarray)):
		str_ls = '{0:.4f}'.format(ls)
		return str_ls
		
	str_ls = ['{0:.4f}'.format(v) for v in ls ]
	str_ls2 = '[' + ' '.join(str_ls) + ']'
	return str_ls2

## pdecay (decay of params) is only used for AdamW and AdamWAMS
def UpdateAlgorithm(alg, params, gparams, pdecay=None, lr=None, l2reg=None):
    	## currently the code mainly work for Adam, AdamW and SGDM. Other algorithms are not fully tested.

	other_params = []
    	#alg = modelSpecs['algorithm']
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
		assert ( lr is not None)
		assert ( l2reg is not None)
	    	updates, other_params = AdamW(params, gparams, pdecay=pdecay, l2reg=l2reg, lr = lr)

	elif alg == 'AdamWAMS':
		assert ( lr is not None)
		assert ( l2reg is not None)
	    	updates, other_params = AdamWAMS(params, gparams, pdecay=pdecay, l2reg=l2reg, lr = lr)

    	else:
        	raise NotImplementedError

	return updates, other_params

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
       		print 'restarting from a check point file: ', checkPointFile
		saved_fh = open(checkPointFile, 'rb')
		chkpoint = cPickle.load( saved_fh )
		saved_fh.close()

		if not chkpoint.has_key('param_values') or not chkpoint.has_key('other_param_values'):
			print 'Error: the check point file does not have param_values or other_param_valuee! Please do not use this file to restart training. '
			exit(-1)

		if Compatible(params, chkpoint['param_values']):
			[ p.set_value(v) for p, v in zip(params, chkpoint['param_values']) ]
		else:
			print 'Error: param_values in the chkpoint file is incompatible with params!'
			exit(-1)

		restart = True

    	else:
		print 'could not find the check point file: ', checkPointFile
		print 'training the model from scratch ...'

	return chkpoint, restart

def QuickValidateAllData(SeqDataset, validate, modelSpecs):
       	losses = []
	errs = []
	if modelSpecs['UseSampleWeight']:
		w4losses = []
		w4errors = []
	else:
		w4losses = None
		w4errors = None

       	for onebatch in SeqDataset:
        	los, err = validate( *onebatch )
    		losses.append(los)
    		errs.append(err)

		##two different batches may have different number of residues and different distribution of labels 
		##so we shall normalize the loss and errors by the weight of different batches
		if modelSpecs['UseSampleWeight']:
			weights = onebatch[ len(onebatch) - len(modelSpecs['responses']) : ]
			w4loss = []
			w4error = []
			for res, w in zip(modelSpecs['responses'], weights):
				wSum = np.sum(w)
				w4loss.append(wSum)
				w4error.extend( [ wSum] * config.responseValueDims[ Response2LabelType(res) ] )
			w4losses.append(w4loss)
			w4errors.append(w4error)
		
	losses = np.array(losses)
	errs = np.array(errs)

       	return np.average(losses, axis=0, weights=w4losses), np.average(errs, axis=0, weights=w4errors)

## the subroutine for one epoch of training
def RunOneEpoch(epoch, chkpoint, params, other_params, trainSeqData, validSeqData, train, validate, modelSpecs):
       	random.shuffle(trainSeqData)
	n_train_batches = len(trainSeqData)

	if epoch > 12:
                validation_frequency = modelSpecs['validation_frequency']
        elif epoch > 10:
                validation_frequency = np.int32(1.5 * modelSpecs['validation_frequency'])
        elif epoch > 5:
                validation_frequency = np.int32(1.8 * modelSpecs['validation_frequency'])
        else:
                validation_frequency = 2 * modelSpecs['validation_frequency']


	t_loss = []
	t_errors = []
	results = []

       	for minibatch_index, onebatch in zip(xrange(len(trainSeqData) ), trainSeqData):
          	#print 'minibatch_index = ', minibatch_index

		input = onebatch
                train_loss, train_errors, param_L2 = train(*input)

            	t_loss.append(train_loss)
	    	t_errors.append(train_errors)

	    	if ( minibatch_index + 1 ) % validation_frequency  == 0:

			valid_loss, valid_errors = QuickValidateAllData(validSeqData, validate, modelSpecs)

			print( 'epoch %2d, minibatch %4d/%4d, train loss %s, train error %s, paramL2 %.2f, valid loss %s, valid error %s' %
		       		( epoch, minibatch_index + 1, n_train_batches, str_display(train_loss), str_display(train_errors), param_L2, str_display(valid_loss), str_display(valid_errors) ) )
		       		#( epoch, minibatch_index + 1, n_train_batches, train_loss, train_errors, param_L2, valid_loss, valid_errors ) )

			result = dict()
			result['train_loss'] = train_loss
			result['train_errors'] = train_errors
			result['valid_loss'] = valid_loss
			result['valid_errors'] = valid_errors

			results.append(result)

			if np.mean(valid_loss) < np.mean(chkpoint['best_validation_loss'] ):
                    		chkpoint['best_validation_loss'] = valid_loss
		    		chkpoint['best_epoch'] = epoch
		    		chkpoint['bestParamValues'] = [param.get_value () for param in params ]
		    		chkpoint['bestOtherParamValues'] = [param.get_value () for param in other_params ]

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

       	## statistics
	avg_train_loss = np.average( t_loss, axis=0  )	
	avg_train_errors = np.average( t_errors, axis=0 )	
	avg_valid_errors = np.average( [ res['valid_errors'] for res in results ], axis=0 )	
	avg_valid_loss = np.average( [ res['valid_loss'] for res in results ], axis=0 )	
	print( 'average result at epoch %i: train loss %s, train error %s, valid loss %s, valid error %s' %
           	 ( epoch, str_display(avg_train_loss), str_display(avg_train_errors), str_display(avg_valid_loss), str_display(avg_valid_errors) ) )

   	## test on the validation data
        #validLoss, validErr, validConfusion = ValidateAllData(validSeqData)
        validLoss, validErr = QuickValidateAllData(validSeqData, validate, modelSpecs)
        print 'valid loss: ', validLoss, 'valid err: ', validErr
	#print 'validConfusion: ', [ vcf[0] for vcf in validConfusion ]



def TrainModel(modelSpecs=None, trainSeqData=None, validSeqData=None, predSeqData=None):
    	if modelSpecs is None:
        	print 'Please provide a model specification for training'
        	sys.exit(-1)

    	if (not trainSeqData) or (not validSeqData):
        	print 'Please provide train or validation data for model training'
		exit(-1)

	propertyPredictor, x, xmask, labelList, weightList = BuildModel(modelSpecs)

        allvariables = [x, xmask]
	assert (labelList is not None and len(labelList)>0 )
	allvariables.extend(labelList)
	allvariables.extend(weightList)

	print '#variables (including response) used in the model: ', len(allvariables)

    	params = propertyPredictor.params
	params4var = propertyPredictor.params4var
	params4mean = list( set(params) - set(params4var) )

    	param_shapes = [ p.get_value().shape for p in params ]
    	param_sizes = map( np.prod, param_shapes )
    	modelSpecs['numParams'] = sum(param_sizes)
    	print 'The model has ', modelSpecs['numParams'], ' parameters.'

	regularizer = 0
    	paramL1 = propertyPredictor.paramL1
    	if modelSpecs.has_key('L1reg') and modelSpecs['L1reg']>0:
        	regularizer += (modelSpecs['L1reg'] * paramL1)

    	paramL2 = propertyPredictor.paramL2
    	if modelSpecs.has_key('L2reg') and modelSpecs['L2reg']>0:
        	regularizer += (modelSpecs['L2reg'] * paramL2)

    	modelSpecs['validation_frequency'] = max( min( len(validSeqData), len(trainSeqData)/8 ), len(trainSeqData)/12 ) 

    	if weightList is not None and len(weightList)>0:
		errors = propertyPredictor.errors(labelList, weightList)
    	else:
        	errors = propertyPredictor.errors(labelList)

	chkpoint, restart = InitializeChkpoint(params, modelSpecs)

	##def RunOneStage(epoch_start, epoch_end, loss, gparams, lr, pdecay=None, algorithm='Adam', startFromBest=(False, False) ):

	## startPosition is a tuple of two entries. The first entry is for params and the second for other_params
	def RunOneStage(epoch_start, epoch_end, loss, gparams, lr, pdecay=None, algorithm='Adam', startPosition=(None, None) ):

		updates, other_params = UpdateAlgorithm(modelSpecs['algorithm'], params, gparams, pdecay, lr=np.float32(lr), l2reg=modelSpecs['L2reg'] )

		"""
		if startFromBest[0] and chkpoint.has_key('bestParamValues') and  Compatible(params, chkpoint['bestParamValues']):
    			[ u.set_value(v) for u, v in zip(params, chkpoint['bestParamValues']) ]

		if startFromBest[1] and chkpoint.has_key('bestOtherParamValues') and Compatible(other_params, chkpoint['bestOtherParamValues']):
    			[ u.set_value(v) for u, v in zip(other_params, chkpoint['bestOtherParamValues']) ]
		"""

		if (startPosition[0] is not None) and Compatible(params, startPosition[0]):
    			[ u.set_value(v) for u, v in zip(params, startPosition[0] ) ]

		if (startPosition[1] is not None) and Compatible(other_params, startPosition[1] ):
    			[ u.set_value(v) for u, v in zip(other_params, startPosition[1] ) ]

    		train =	theano.function(allvariables, [loss, errors, paramL2], updates=updates, on_unused_input='warn')
    		validate = theano.function(allvariables, [loss, errors], on_unused_input='warn')

		epoch = epoch_start
    		while (epoch < epoch_end):
        		epoch += 1
			RunOneEpoch(epoch, chkpoint, params, other_params, trainSeqData, validSeqData, train, validate, modelSpecs)
			if epoch - chkpoint['best_epoch'] > modelSpecs['patience']:
				break
		return epoch

	assert ( len(modelSpecs['numEpochs'] ) > 0 )
        numEpochs4stages = np.cumsum(modelSpecs['numEpochs'] )
        epoch = chkpoint['epoch']

	## train parameters not related to variance and correlation
	##if epoch < modelSpecs['numEpochs']:
	if epoch < numEpochs4stages[-1] :
    		if weightList is not None and len(weightList)>0:
        		loss = propertyPredictor.loss(labelList, True, weightList=weightList)
    		else:
			loss = propertyPredictor.loss(labelList, True)

		cost = T.sum( T.mul(loss, modelSpecs['w4responses']) ) / np.sum(modelSpecs['w4responses']) 
		if modelSpecs['algorithm'] not in set(['AdamW', 'AdamWAMS']):
			cost = cost + regularizer

		params4var_set = set(params4var)
		gparams = [ T.grad(cost, p, consider_constant=weightList, disconnected_inputs='raise' ) if p not in params4var_set else T.zeros_like(p) for p in params ]
		pdecay = [ p if p not in params4var_set else T.zeros_like(p) for p in params ]

	for lr, epoch_end in zip(modelSpecs['lrs'], numEpochs4stages):
                if epoch >= epoch_end:
                        continue

                epoch_start = epoch
                print 'training for mean using a learning rate ', lr, ' ...'
		if restart:
			startPosition = (chkpoint['param_values'], chkpoint['other_param_values'])
		else:
			startPosition = (chkpoint['bestParamValues'], chkpoint['bestOtherParamValues'])

                epoch = RunOneStage(epoch_start, epoch_end, loss, gparams, lr, pdecay=pdecay, algorithm=modelSpecs['algorithm'], startPosition=startPosition )


	## train parameters only specific to variance and correlation
	numEpochs4var = modelSpecs['numEpochs4var']
	lrs = modelSpecs['lrs4var']

	if len(params4var) >0:
		assert ( len(numEpochs4var) > 0 )
                assert ( len(lrs) > 0 )

                numEpochs4stages = np.cumsum( numEpochs4var ) + numEpochs4stages[-1]

		if epoch < numEpochs4stages[-1]:
			print 'Training the parameters specific for correlation and variance ...'

    			if weightList is not None and len(weightList)>0:
        			loss = propertyPredictor.loss(labelList, weightList = weightList)
    			else:
				loss = propertyPredictor.loss(labelList)

			cost = T.sum( T.mul(loss, modelSpecs['w4responses']) ) / np.sum(modelSpecs['w4responses']) 

			if modelSpecs['algorithm'] not in set(['AdamW', 'AdamWAMS']):
				cost = cost + regularizer

			params4var_set = set(params4var)
			gparams = [ T.grad(cost, p, consider_constant=weightList, disconnected_inputs='raise' ) if p in params4var_set else T.zeros_like(p) for p in params ]
			pdecay = [  p if p in params4var_set else T.zeros_like(p) for p in params ]

		stage = 0
		for lr, epoch_end in zip(lrs, numEpochs4stages):
                        if epoch >= epoch_end:
                                continue
			stage += 1
			epoch_start = epoch
			print  'training for variance using a learning rate ', lr, ' ...'
			if restart:
				startPosition = (chkpoint['param_values'], chkpoint['other_param_values'])
			elif stage == 1:
				startPosition = (chkpoint['bestParamValues'], None)
			else:
				startPosition = (chkpoint['bestParamValues'], chkpoint['bestOtherParamValues'] )

			epoch = RunOneStage(epoch_start, epoch_end, loss, gparams, lr, pdecay=pdecay, algorithm=modelSpecs['algorithm4var'], startPosition=startPosition  )


    	if weightList is not None and len(weightList)>0:
        	loss = propertyPredictor.loss(labelList, weightList = weightList)
    	else:
		loss = propertyPredictor.loss(labelList)


    	validate = theano.function(allvariables, [loss, errors], on_unused_input='warn')

    	## set model parameters for valiation and possibly prediction
    	for param, value in zip(params, chkpoint['bestParamValues']):
       		param.set_value(value)

    	resultModel = {}
    	resultModel['dateTrained']=datetime.datetime.now()

    	## test on the validation data
    	validLoss, validErr = QuickValidateAllData(validSeqData, validate, modelSpecs)
    	print 'valid loss: ', validLoss, 'valid err: ', validErr
    	resultModel['validLoss'] = validLoss
    	resultModel['validErr'] = validErr

    	resultModel['paramValues']=chkpoint['bestParamValues']

    	bestParamL2norm = np.sum([(v**2).sum() for v in chkpoint['bestParamValues'] ])
    	resultModel['bestParamL2norm']=bestParamL2norm

    	bestParamL1norm = np.sum([abs(v).sum() for v in chkpoint['bestParamValues'] ])
    	resultModel['bestParamL1norm']=bestParamL1norm

    	print 'best param L1 norm: ', bestParamL1norm, 'L2 norm: ', bestParamL2norm

    	#test on prediction data if it is given
     	if predSeqData is not None:

        	predLoss, predErr = QuickValidateAllData(predSeqData, validate, modelSpecs)
        	print 'pred loss: ', predLoss, 'pred err: ', predErr
        	resultModel['predLoss'] = predLoss
        	resultModel['predErr'] = predErr

    	return resultModel


##automatically generating a model file name
def GenerateModelFileName(resultModel):

	prefix = ''
	if resultModel.has_key('UseTemplate') and resultModel['UseTemplate']:
		prefix += 'TPL'
	else:
		prefix += 'Seq'

	prefix += resultModel['network'] + '4'
	prefix += '.'.join(resultModel['responses']).replace('Discrete', '')

	arch =  'L' + str( (sum(resultModel['conv1d_repeats']) + len(resultModel['conv1d_hiddens']) )*2 -1 )  
	arch += 'Log' + str( (sum(resultModel['logreg_hiddens']) + len(resultModel['logreg_hiddens']) )      )  
	arch += 'W' + str(resultModel['halfWinSize_seq']) 
	arch += 'I' + str(resultModel['n_in_seq'])

	epoch = resultModel['algStr']
	suffix = '.pkl'

	datastr = os.path.basename(resultModel['trainFile'][0]).split('.')[0:-2]
        datastr = ''.join(datastr)


	pid = str(os.getpid())
	components = [ prefix, arch+epoch, datastr, pid]
	
	filename = os.path.join('LocalModels/', ('-'.join(components)).replace(';', '+') + suffix )
	if not os.path.isdir('LocalModels'):
		os.mkdir('LocalModels')

	return filename

def main(argv):

    	#modelSpecs = config.InitializeModelSpecs()
    	modelSpecs = InitializeModelSpecs()
	modelSpecs = ParseCommandLine.ParseArguments(argv, modelSpecs)

	startTime = datetime.datetime.now()

        ##trainData and validData are a list. Each element corresponds to one protein, which is a dict() 
	trainData = DataProcessor.LoadPropertyFeatures(modelSpecs['trainFile'], modelSpecs=modelSpecs )
        validData = DataProcessor.LoadPropertyFeatures(modelSpecs['validFile'], modelSpecs=modelSpecs )
        print '#trainData: ', len(trainData), '#validData: ', len(validData)

	## where to add code to assign weight to each residue? We need to deal with the residues without 3D coordinates for angle and SS prediction
	##a, b = DataProcessor.CalcLabelDistributionAndWeight(trainData, modelSpecs)

	modelSpecs['numOfTrainProteins']= len(trainData)

	beforeBatchTime = datetime.datetime.now()
	print 'time spent on data loading: ', beforeBatchTime - startTime

        print 'Preparing batch data for training...'
        groupSize = modelSpecs['minibatchSize']
        trainSeqDataset, _ = DataProcessor.SplitData2Batches(data=trainData, numDataPoints=groupSize, modelSpecs=modelSpecs)
        validSeqDataset, _ = DataProcessor.SplitData2Batches(data=validData, numDataPoints=groupSize, modelSpecs=modelSpecs)
        #validSeqDataset = DataProcessor.SplitData2Batches(data=validData, numDataPoints=20000, modelSpecs=modelSpecs)
        print "#trainData minibatches:", len(trainSeqDataset), "#validData minibatches:", len(validSeqDataset)

        predSeqDataset = None
	if modelSpecs['predFile'] is not None:
            	predData = DataProcessor.LoadPropertyFeatures(modelSpecs['predFile'], modelSpecs=modelSpecs, forTrainValidation=False )
	    	print '#predData: ', len(predData)
	 	predSeqDataset, _ = DataProcessor.SplitData2Batches(data=predData, numDataPoints=40, modelSpecs=modelSpecs)
                print "#predData minibatches:", len(predSeqDataset)

	## Each protein in trainData contains three or four components: seqFeatures and label
        modelSpecs['n_in_seq'] = trainData[0]['seqFeatures'].shape[1]

	beforeTrainTime = datetime.datetime.now()

	print 'time spent on generating batch data:', beforeTrainTime - beforeBatchTime

        result = TrainModel(modelSpecs=modelSpecs, trainSeqData=trainSeqDataset, validSeqData=validSeqDataset, predSeqData=predSeqDataset)

        ##merge ModelSpecs and result
        resultModel = modelSpecs.copy()
        resultModel.update(result)

	modelFile = GenerateModelFileName(resultModel)
	print 'Writing the resultant model to ', modelFile
        cPickle.dump(resultModel, file(modelFile, 'wb'), cPickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
   	recursionlimit = sys.getrecursionlimit()
   	print 'recursionlimit = ', recursionlimit
   	sys.setrecursionlimit(3*recursionlimit)

   	main(sys.argv[1:])
