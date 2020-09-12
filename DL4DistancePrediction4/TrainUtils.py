import os
import sys
import random
import datetime
import copy
import numpy as np
import math

import config
import DataProcessor
import FeatureUtils
import LabelUtils

import theano

from utilsNoT import SampleBoundingBox

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

## this function is used to skip a portion of one epoch. It reads a file SkipOneEpoch.pid to obtain which epoch to skip where pid is the process id of the training program. 
## For example, if the file has a float 16.4, then the 17th epoch will be skipped after training its 40% 
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



##automatically generating a model file name
def GenerateModelFileName(resultModel):

        prefix = ''
        verStr ='Model' +  str(int(resultModel['Version']*10)) + '-'
        if config.UseTemplate(resultModel):
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

        datastr = os.path.basename(resultModel['trainFile']).split('.')[0:-2]
        datastr = ''.join(datastr)

        #components = [ prefix, arch+bias+epoch, datastr, pid, resultModel['algorithm'] ]
        components = [ prefix, arch, datastr, MID, resultModel['algorithm'] ]
        filename = os.path.join('LocalModels/', '-'.join(components) + suffix )

        if not os.path.exists('LocalModels/'):
                os.makedirs('LocalModels/')

        return filename


def PrepareInput4Train(data, modelSpecs, floatType=np.float32, forRefState=False, UseSharedMemory=False):
	if not bool(data):
		print 'ERROR: the input data for PrepareInput4Train2 is empty'
		exit(1)

	allowedLen = int(math.floor(math.sqrt(modelSpecs['maxbatchSize']) ) )
	bounds =[]
	for d in data:
		if d['seqLen'] < allowedLen:
			bounds.append( None )
			continue
		box = SampleBoundingBox( (d['seqLen'], d['seqLen']), modelSpecs['maxbatchSize'] )
		bounds.append(box)
	
	#print allowedLen
	#print bounds
	onebatch, _= DataProcessor.AssembleOneBatch(data, modelSpecs, forRefState=forRefState, bounds=bounds, floatType=floatType, bUseSharedMemory=UseSharedMemory)

	## determine the bounding box. 
	maxSeqLen = max([ d['seqLen'] for d in data ])
	#print maxSeqLen

	if maxSeqLen > allowedLen and len(data)>1:
		print 'ERROR: one minibatch has more than one large proteins: ', [ d['name'] for d in data ]
		exit(1)
		
	if maxSeqLen <= allowedLen:
		box = np.array([0, 0, maxSeqLen, maxSeqLen]).astype(np.int32)
	else:
		## in this case, len(data) == 1 and len(bounds) == 1
		assert bounds[0] is not None
		box = np.array(bounds[0]).astype(np.int32)

	onebatch.append(box)

	if config.TrainByRefLoss(modelSpecs):
                if forRefState:
                        onebatch.append(np.int32(-1) )
                else:
                        onebatch.append(np.int32(1) )
	return onebatch

def PrepareInput4Validate(data, modelSpecs, floatType=np.float32, forRefState=False, UseSharedMemory=False):
	if not bool(data):
		print 'ERROR: the input data for PrepareInput4Validate is empty'
		exit(1)

	if UseSharedMemory:
		## when shared memory is used, there is no explicit limit on the size of an ndarray
		maxAllowedLen = np.iinfo(np.int32).max
	else:
		## when the real content of a large matrix is passed through Queue, its size shall be <2GB
		maxAllowedLen = 800

	maxSeqLen = max([ d['seqLen'] for d in data ])

	if maxSeqLen > maxAllowedLen and len(data)>1:
		print 'ERROR: when one validation protein has length > ', maxAllowedLen, ', it shall form a minibatch by itself'
		exit(1)

	##determine the bounding box. 
	if maxSeqLen <= maxAllowedLen:
		bounds = None
	else:
		bounds = []
		for d in data:
			seqLen = d['seqLen']
			if seqLen > maxAllowedLen:
				## cut off a submatrix along the diagonal line so that the top accurcy function in our deep model works correctly
				top = 0
				bottom = maxAllowedLen
				#left = seqLen - maxAllowedLen
				left = 0
				#right = seqLen
				right = maxAllowedLen
				box = [top, left, bottom, right]
				bounds.append(box)
			else:
				bounds.append(None)

	onebatch, _= DataProcessor.AssembleOneBatch(data, modelSpecs, forRefState=forRefState, bounds=bounds, floatType=floatType, bUseSharedMemory=UseSharedMemory)
	if maxSeqLen <= maxAllowedLen:
		box = np.array([0, 0, maxSeqLen, maxSeqLen]).astype(np.int32)
	else:
		## in this case, len(bounds)==1 and len(data) == 1
		assert len(bounds)==1
		assert bounds[0] is not None
		box = np.array(bounds[0]).astype(np.int32)

	onebatch.append( box )

	return onebatch

def PrepareInput4Prediction(data, modelSpecs, floatType=np.float32, UseSharedMemory=False, forRefState=False):
	if not bool(data):
		print 'ERROR: the input data for PrepareInput4Prediction is empty'
		exit(1)

	onebatch, _= DataProcessor.AssembleOneBatch(data, modelSpecs, forRefState=forRefState, floatType=floatType, bUseSharedMemory=UseSharedMemory)
	maxSeqLen = max( [ d['seqLen'] for d in data ] )	
	box = np.array([0, 0, maxSeqLen, maxSeqLen]).astype(np.int32)
	onebatch.append( box )

	return onebatch

## this function loads data from disk and/or remote storage server to a shared queue
## sharedQ is the shared Queue with a limited size
##validSeqData is a list of batches. Each batch has location info for a list of validation proteins
def ValidDataLoader(sharedQ, validSeqData, modelSpecs, assembleData=True, UseSharedMemory=False):
	for batch in validSeqData:
		 ## Load real data for one batch
                data = DataProcessor.LoadRealData(batch, modelSpecs, returnMode='list')

                ## add code here to make sure that the data has the same input dimension as the model specification
                FeatureUtils.CheckModelNDataConsistency(modelSpecs, data)

		if assembleData:
			data = PrepareInput4Validate(data, modelSpecs, floatType=np.float16, UseSharedMemory=UseSharedMemory)
		#print 'putting data to validDataLoader queue...'
		sharedQ.put(data)

## this function loads data from disk and/or remote storage server to a shared queue
## sharedQ is the shared Queue with a limited size
##validSeqData is a list of batches. Each batch has location info for a list of validation proteins
def ValidDataLoader2(sharedQ, stopValidDataLoader, validSeqData, modelSpecs, assembleData=True, UseSharedMemory=False):

	bUseCCMFnorm, bUseCCMsum, bUseCCMraw, bUseFullMI, bUseFullCov = config.ParseExtraCCMmode(modelSpecs)
	if any([bUseCCMraw, bUseFullMI, bUseFullCov]):
		## when full coevolution matrices are used, we shall use float16 to save memory
		floatType = np.float16
	else:
		floatType = theano.config.floatX

	#print 'validDataLoader has event: ', stopValidDataLoader
	for batch in validSeqData:
		if stopValidDataLoader.is_set() or os.getppid()==1:
			#print 'validDataLoader receives the stop signal'
			break

		 ## Load real data for one batch
                data = DataProcessor.LoadRealData(batch, modelSpecs, returnMode='list')

                ## add code here to make sure that the data has the same input dimension as the model specification
                FeatureUtils.CheckModelNDataConsistency(modelSpecs, data)

		if assembleData:
			data = PrepareInput4Validate(data, modelSpecs, floatType=floatType, UseSharedMemory=UseSharedMemory)
		#print 'putting data to validDataLoader queue...'
		sharedQ.put(data)

	print 'validDataLoader has finished loading data'
	sharedQ.close()
	
## similar to ValidDataLoader, but for train data	
def TrainDataLoader(sharedQ, trainMetaData, modelSpecs, assembleData=True, UseSharedMemory=False):
	## here we use labelPool to cache the labels of all the training proteins
	## one protein may have multiple sets of input features due to MSA sampling or sequnence-template alignment
	## but it can only have one set of label matrices, so it is worth to save all label matrices in RAM.
	labelPool = dict()
	labelMatrixPool = dict()

	while True:
		trainDataLocation = DataProcessor.SampleProteinInfo(trainMetaData)
		numOriginals = len(trainDataLocation)
		trainSeqData = DataProcessor.SplitData2Batches(trainDataLocation, numDataPoints=modelSpecs['minibatchSize'], modelSpecs=modelSpecs)
		random.shuffle(trainSeqData)
		for batch in trainSeqData:
			data = []
			for protein in batch:
				name = protein['name']
				if labelPool.has_key(name):
					## label is already in the pool
					d = DataProcessor.LoadRealData(protein, modelSpecs, loadLabel=False, returnMode='list')
					d['atomLabelMatrix'] = labelPool[name]
				else:
					d = DataProcessor.LoadRealData(protein, modelSpecs, returnMode='list')
					assert d.has_key('atomLabelMatrix')
					labelPool[name] = d['atomLabelMatrix']

				if config.UseSampleWeight(modelSpecs):
					if not labelMatrixPool.has_key(name): 
						labelWeightMatrix = LabelUtils.CalcLabelWeightMatrix(LabelMatrix=d['atomLabelMatrix'], modelSpecs=modelSpecs, floatType=np.float16)
						labelMatrixPool[name] = labelWeightMatrix
						d['labelWeightMatrix'] = labelWeightMatrix
					else:
						d['labelWeightMatrix'] = labelMatrixPool[name]

				data.append(d)

			FeatureUtils.CheckModelNDataConsistency(modelSpecs, data)
			if assembleData:
				data = PrepareInput4Train(data, modelSpecs, floatType=np.float16, UseSharedMemory=UseSharedMemory)
			#print 'putting data to trainDataLoader queue...'
			sharedQ.put(data)

		#print 'TrainDataLoader with #PID ', os.getpid(), ' currently has ', len(labelPool), ' label matrices  and ', len(labelMatrixPool), ' label weight matrices'

## similar to ValidDataLoader, but for train data	
def TrainDataLoader2(sharedQ, stopTrainDataLoader, trainMetaData, modelSpecs, assembleData=True, UseSharedMemory=False):
	#print 'trainDataLoader has event: ', stopTrainDataLoader

	bUseCCMFnorm, bUseCCMsum, bUseCCMraw, bUseFullMI, bUseFullCov = config.ParseExtraCCMmode(modelSpecs)
	if any([bUseCCMraw, bUseFullMI, bUseFullCov]):
		## when full coevolution matrices are used, we shall use float16 to save memory
		floatType = np.float16
	else:
		floatType = theano.config.floatX

	## here we use labelPool to cache the labels of all the training proteins
	## one protein may have multiple sets of input features due to MSA sampling or sequnence-template alignment
	## but it can only have one set of label matrices, so it is worth to save all label matrices in RAM.
	labelPool = dict()
	labelWeightPool = dict()

	while True:
		if stopTrainDataLoader.is_set() or os.getppid()==1:
			#print 'trainDataLoader receives the stop signal'
			break

		trainDataLocation = DataProcessor.SampleProteinInfo(trainMetaData)
		numOriginals = len(trainDataLocation)
		trainSeqData = DataProcessor.SplitData2Batches(trainDataLocation, numDataPoints=modelSpecs['minibatchSize'], modelSpecs=modelSpecs)
		random.shuffle(trainSeqData)

		#i = 0
		for batch in trainSeqData:
			if stopTrainDataLoader.is_set() or os.getppid()==1:
				#print 'trainDataLoader receives the stop signal'
				break

			data = []
			for protein in batch:
				name = protein['name']
				if labelPool.has_key(name):
					## label is already in the pool
					d = DataProcessor.LoadRealData(protein, modelSpecs, loadLabel=False, returnMode='list')
					d['atomLabelMatrix'] = labelPool[name]
				else:
					d = DataProcessor.LoadRealData(protein, modelSpecs, returnMode='list')
					assert d.has_key('atomLabelMatrix')
					labelPool[name] = d['atomLabelMatrix']

				if config.UseSampleWeight(modelSpecs):
					if not labelWeightPool.has_key(name): 
						labelWeightMatrix = LabelUtils.CalcLabelWeightMatrix(LabelMatrix=d['atomLabelMatrix'], modelSpecs=modelSpecs, floatType=np.float16)
						labelWeightPool[name] = labelWeightMatrix
						d['labelWeightMatrix'] = labelWeightMatrix
					else:
						d['labelWeightMatrix'] = labelWeightPool[name]

				data.append(d)

			FeatureUtils.CheckModelNDataConsistency(modelSpecs, data)
			if assembleData:
				data = PrepareInput4Train(data, modelSpecs, floatType=floatType, UseSharedMemory=UseSharedMemory)
			#print 'putting data to trainDataLoader queue...'
			sharedQ.put(data)

			"""
			i += 1
			if i%100 == 0:
				print '#batches of train data loaded: ', i
			"""

		#print 'TrainDataLoader with #PID ', os.getpid(), ' currently has ', len(labelPool), ' label matrices  and ', len(labelMatrixPool), ' label weight matrices'
	print 'TrainDataLoader has finished loading data'
	sharedQ.close()

## similar to ValidDataLoader, but for train data	
def TrainDataLoader3(sharedQ, sharedLabelPool, sharedLabelWeightPool, stopTrainDataLoader, trainMetaData, modelSpecs, assembleData=True, UseSharedMemory=False):
	#print 'trainDataLoader has event: ', stopTrainDataLoader

	## here we use labelPool to cache the labels of all the training proteins
	## one protein may have multiple sets of input features due to MSA sampling or sequnence-template alignment
	## but it can only have one set of label matrices, so it is worth to save all label matrices in RAM.
	labelPool = dict()
	labelWeightPool = dict()

	## load the labels of all training proteins
	trainDataLocation = DataProcessor.SampleProteinInfo(trainMetaData)
	for loc in trainDataLocation:
		d = DataProcessor.LoadRealData(loc, modelSpecs, loadFeature=False, returnMode='list')
		name = d['name']
		labelPool[name] = d['atomLabelMatrix']
		labelWeightMatrix = LabelUtils.CalcLabelWeightMatrix(LabelMatrix=d['atomLabelMatrix'], modelSpecs=modelSpecs, floatType=np.float16)
		labelWeightPool[name] = labelWeightMatrix

	print 'TrainDataLoader with #PID ', os.getpid(), ' has loaded ', len(labelPool), ' label matrices  and ', len(labelWeightPool), ' label weight matrices'
	## update labelPool and labelWeightPool to the shared dict()
	sharedLabelPool.update(labelPool)
	sharedLabelWeightPool.update(labelWeightPool)
	print 'TrainDataLoader with #PID ', os.getpid(), ' has update the shared labelPool and labelWeightPool'

	while True:
		if stopTrainDataLoader.is_set() or os.getppid()==1:
			print 'trainDataLoader receives the stop signal'
			break

		trainDataLocation = DataProcessor.SampleProteinInfo(trainMetaData)
		numOriginals = len(trainDataLocation)
		"""
		maxLen = 900
		trainDataLocation, numExcluded = DataProcessor.FilterByLength(trainDataLocation, maxLen)
		print 'Exclude ', numExcluded, ' train proteins longer than ', maxLen, ' AAs'
		"""
		trainSeqData = DataProcessor.SplitData2Batches(trainDataLocation, numDataPoints=modelSpecs['minibatchSize'], modelSpecs=modelSpecs)
		random.shuffle(trainSeqData)
		for batch in trainSeqData:
			if stopTrainDataLoader.is_set() or os.getppid()==1:
				print 'trainDataLoader receives the stop signal'
				break

			names = [ p['name'] for p in batch ]
			data = []
			for protein in batch:
				d = DataProcessor.LoadRealData(protein, modelSpecs, loadLabel=False, returnMode='list')
				data.append(d)

			FeatureUtils.CheckModelNDataConsistency(modelSpecs, data)
			if assembleData:
				data = PrepareInput4Train(data, modelSpecs, floatType=np.float16, UseSharedMemory=UseSharedMemory)
			#print 'putting data to trainDataLoader queue...'
			sharedQ.put( (data, names) )

	print 'TrainDataLoader has finished loading data'
	sharedQ.close()


def FetchOneBatch(sharedQ, modelSpecs, forTrain=True, assembleData=False, UseSharedMemory=False):
	data = sharedQ.get()
	if not assembleData:
		## data is already assembled, so just return it
		return data

	if forTrain:
		## prepare input data for train
		data = PrepareInput4Train(data, modelSpecs, floatType=np.float16)
	else:
		## prepare input data for validation
		data = PrepareInput4Validate(data, modelSpecs, floatType=np.float16)

	return data

