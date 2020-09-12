import os
import random
import cPickle
import numpy as np
import json
import copy

import theano

## shared memory implementation of numpy ndarray
from shared_ndarray import SharedNDArray

from Alignment import AlignmentUtils
from Common import NativeUtils
from Common.SequenceUtils import SeqOneHotEncoding

import DistanceUtils
import OrientationUtils
import FeatureUtils
import LabelUtils
from utilsNoT import RowWiseOuterProduct
from utilsNoT import SampleBoundingBox, SplitList 
from utilsNoT import GenRandomString as randomString

import config
from config import Response2LabelName, Response2LabelType, ParseResponse, GetResponseProbDims, GetResponseValueDims

## dataFile can be loaded by json
def LoadMetaData(dataFile):
        if not os.path.isfile(dataFile):
                print 'ERROR: meta data file does not exist: ', dataFile
                exit(1)
        with open(dataFile) as fh:
                metaData = json.load(fh)
        return metaData

## this function samples one set of input features from numSamples groups in metaData
## one group has one or multiple related proteins and/or one protein with multiple templates. 
## At each training epoch we sample one sequence or one sequence-template pair from a group.
def SampleProteinInfo(metaData, numSamples=None):
	## MetaData is a dict() with the following keys: groups, SeqDir, DistDir, OriDir and a set of infeatureDir
	## MetaData['groups'] is a list of groups, each group is a dict() with keys: group name, proteins, sequences, weights and templates
	## when not None, templates are used to indicate if template-based information shall be used for this protein or not
	
	## DistDir and OriDir are the folders for true distance and orientation information
	## infeatureDir is a list of folders for input features

	## input features, ground truth and template information shall be saved in PKL format

	## when numSamples is not None, obtain this number of samples. Otherwise obtain one protein from each group.
	featureFolders = metaData['featureFolders']

	if numSamples is not None and numSamples < len(metaData['groups']):
		groups = random.sample(metaData['groups'], min(numSamples, len(metaData['groups'])) )
	else:
        	groups = metaData['groups']

        proteins = []
        for group in groups:

                protein = dict()

		if metaData.has_key('nativeDir'):
			protein['nativeLocation'] = metaData['nativeDir']

		if metaData.has_key('distDir'):
                	protein['distLocation'] = metaData['distDir']

		if metaData.has_key('oriDir'):
                	protein['oriLocation'] = metaData['oriDir']
	
		if metaData.has_key('aliDir'):
			protein['aliLocation'] = metaData['aliDir']

		if metaData.has_key('tplDir'):
			protein['tplLocation'] = metaData['tplDir']

		if metaData.has_key('pdbDir'):
			protein['pdbLocation'] = metaData['pdbDir']

		## sample one protein from the group by weights
		prob = group['weights']/np.float32( np.sum(group['weights']) )
		idx = np.random.choice(np.arange(len(prob)), p=prob)
                protein['name'] = group['proteins'][idx]
                protein['sequence'] = group['sequences'][idx]
		protein['templates'] = group['templates'][idx]
                protein['seqLen'] = len(protein['sequence'])

		## sample one folder for input feature
                idx = np.random.choice( np.arange(len(featureFolders) ) )
                protein['featureLocation'] = featureFolders[idx]

                proteins.append(protein)

        return proteins

## split one metaData into numParts subsets with equal number of groups
def SplitMetaData(rawMetaData, numParts):
	groups = copy.deepcopy( rawMetaData['groups'] )
	random.shuffle(groups)

	subgroups = SplitList(groups, numParts)

	metaDatas = [ dict() ] * numParts
	for sg, md in zip(subgroups, metaDatas):
		md['groups'] = sg
		## copy other information
		for k, v in rawMetaData.iteritems():
			if k == 'groups':
				continue
			md[k] = v
	return metaDatas

def CollectTemplateMatrixFeatures(d, modelSpecs):
	#print 'Using template distance and orientation matrix...'
	if not d.has_key('tplDistMatrix'):
		print 'ERROR: the data for ', d['name'], ' has no tplDistMatrix, which is needed since you want to use template information'
		exit(1)

	if not d['tplDistMatrix'].has_key('CbCb'):
		print 'ERROR: tplDistMatrix shall have distance matrices for atom pairs CbCb'
		exit(1)

	CbCbMatrix = d['tplDistMatrix']['CbCb']

	templateMatrixFeatures = [ ]
	for response in modelSpecs['responses']:
		labelName, labelType, subType = config.ParseResponse(response)

		if labelName in config.allAtomPairNames:
			## process template dist matrix
			tplDistMatrix = d['tplDistMatrix'][labelName]
                	strengthMatrix = np.copy(tplDistMatrix)
                	np.putmask(strengthMatrix, tplDistMatrix<2.0, 2.0)
                	strengthMatrix = 2.0/strengthMatrix
                	## for invalid entry (no valid coordinates or insertion in the sequence), assign the strength to 0
                	np.putmask(strengthMatrix, tplDistMatrix<0, 0)

                	templateMatrixFeatures.append(strengthMatrix)
			continue

		if labelName in config.allOrientationNames:
			## process template orientation matrix
			if not d.has_key('tplOriMatrix'):
				print 'ERROR: the data for ', d['name'], ' has no tplOriMatrix, which is needed since you want to use template orientation information'
				exit(1)
			oriMatrix = d['tplOriMatrix'][labelName]
			## discretize and convert to one-hot encoded matrix
			oriBins=config.GetCutoffs(response)
			invalidEntrySeparated=subType.endswith('Plus') or subType.endswith('Minus')  
			#labelMatrix, _, _ = OrientationUtils.DiscretizeOrientationMatrix(oriMatrix, bins=oriBins, distThreshold4Orientation=modelSpecs['distThreshold4Orientation'], distMatrix=CbCbMatrix, invalidEntrySeparated=invalidEntrySeparated )
			labelMatrix, _, _ = OrientationUtils.DiscretizeOrientationMatrix(oriMatrix, bins=oriBins, distThreshold4Orientation=modelSpecs['distThreshold4Orientation'], distMatrix=CbCbMatrix, invalidEntrySeparated=invalidEntrySeparated, asResponse=False)
			oneHotLabelMatrix = LabelUtils.GetOneHotLabelMatrix(labelMatrix, numLabels=config.GetResponseProbDims(response))

			templateMatrixFeatures.append(oneHotLabelMatrix)
			continue

		print 'WARNING: template information for response not implemented: ', response
			
	return templateMatrixFeatures

def CollectMatrixFeatures(d, modelSpecs, returnMode='array'):
	##collecting pairwise features...
	pairfeatures_nomean = [] # a set of pairwise features for which we do not calculate their expected value
	pairfeatures = []

	if not config.NoOldLocationFeatures(modelSpecs):
		posFeature = FeatureUtils.LocationFeature(d)
        	pairfeatures_nomean.append(posFeature)
		
		cbrtFeature = FeatureUtils.CubeRootFeature(d)
		pairfeatures_nomean.append(cbrtFeature)

	if config.UseNewLocationFeatures(modelSpecs):
		posFeatures = FeatureUtils.NewLocationFeature(d)
        	pairfeatures_nomean.extend(posFeatures)

	if config.UseCCMZ(modelSpecs):
		if not d.has_key('ccmpredZ'):
	    		print 'ERROR: CCMpredZ is requested, but the data for protein ', d['name'], ' does not have it!'
	    		exit(1)
		else:
	    		pairfeatures.append( d['ccmpredZ'] )

	if config.UseRawCCM(modelSpecs):
		if not d.has_key('ccmpred'):
	    		print 'ERROR: Raw CCMpred is requested, but the data for protein ', d['name'], ' does not have it!'
	    		exit(1)
	    	pairfeatures.append( d['ccmpred'] )

        if config.UsePSICOV(modelSpecs):
		if not d.has_key('psicovZ'):
			print 'ERROR: psicovZ is requested, but the data for protein ', d['name'], ' does not have it'
			exit(1)
	    	pairfeatures.append(d['psicovZ'])

	if config.UseContactPotential(modelSpecs):
		if not d.has_key('OtherPairs'):
			print 'ERROR: pairwise contact potential is requested, but the data for protein ', d['name'], ' does not have it'
			exit(1)
	    	pairfeatures.append( d['OtherPairs'][:,:,0] )

	if config.UseMI(modelSpecs):
		if not d.has_key('OtherPairs'):
			print 'ERROR: mutual information is requested, but the data for protein ', d['name'], ' does not have it'
			exit(1)
	    	pairfeatures.append( d['OtherPairs'][:,:,1:3] )

	bUseCCMFnorm, bUseCCMsum, bUseCCMraw, bUseFullMI, bUseFullCov = config.ParseExtraCCMmode(modelSpecs)
	if bUseCCMFnorm:
		if not d.has_key('CCMFnorm') or not d.has_key('CCMFnormZ'):
			print 'ERROR: CCM Fnorm and/or FnormZ are requested, but the data for protein ', d['name'], ' does not have it'
			exit(1)
		pairfeatures.append(d['CCMFnorm'])
		pairfeatures.append(d['CCMFnormZ'])

	if bUseCCMsum:
		if not d.has_key('sumCCM'):
			print 'ERROR: CCM summary are requested, but the data for protein ', d['name'], ' does not have it'
			exit(1)
		pairfeatures.append(d['sumCCM'])

	if bUseCCMraw:
		if not d.has_key('rawCCM'):
			print 'ERROR: CCM raw matrix is requested, but the data for protein ', d['name'], ' does not have it'
			exit(1)
		pairfeatures.append(d['rawCCM'])

	if bUseFullMI:
		if not d.has_key('fullMI'):
			print 'ERROR: full MI matrix is requested, but the data for protein ', d['name'], ' does not have it'
			exit(1)
		pairfeatures.append(d['fullMI'])

	if bUseFullCov:
		if not d.has_key('fullCov'):
			print 'ERROR: full covariance matrix is requested, but the data for protein ', d['name'], ' does not have it'
			exit(1)
		pairfeatures.append(d['fullCov'])

	##add template-based distance and orientation matrices
	if config.UseTemplate(modelSpecs):
		pairfeatures_nomean.extend( CollectTemplateMatrixFeatures(d, modelSpecs) )

	if returnMode.lower() == 'array':
        	matrixFeature = np.dstack( tuple(pairfeatures) )
		if len(pairfeatures_nomean) > 0:
        		matrixFeature_nomean = np.dstack( tuple(pairfeatures_nomean) )
		else:
			seqLen = matrixFeature.shape[0]
			matrixFeature_nomean = np.zeros( (seqLen, seqLen, 0), dtype=config.MyFloat)

        	#print 'matrixFeature.shape: ', matrixFeature.shape
		return matrixFeature, matrixFeature_nomean
	else:
		return pairfeatures, pairfeatures_nomean

def CollectSequentialFeatures(d, modelSpecs, oneHotEncoding, returnMode='array'):
	##collecting sequential features...
	seqMatrices = []

	if config.UseOneHotEncoding(modelSpecs):
		seqMatrices.append(oneHotEncoding)

	if modelSpecs.has_key('UseSS') and (modelSpecs['UseSS'] is True ):
	   	seqMatrices.append( d['SS3'] )

	if modelSpecs.has_key('UseACC') and (modelSpecs['UseACC'] is True ) :
	    	seqMatrices.append( d['ACC'] )

	if modelSpecs.has_key('UsePSSM') and (modelSpecs['UsePSSM'] is True ) :
        	seqMatrices.append( d['PSSM'] )

	if modelSpecs.has_key('UseDisorder') and modelSpecs['UseDisorder'] is True:
	    	seqMatrices.append( d['DISO'] )

	if config.UseRgAsSequentialFeature(modelSpecs):
		seqMatrices.append(Rg(d, AsSequentialFeature=True))

	##membrane protein specific features
        useMPSpecificFeatures = modelSpecs.has_key('UseMPSpecificFeatures') and (modelSpecs['UseMPSpecificFeatures'] is True)
	if useMPSpecificFeatures:
		if d.has_key('MemAcc'):
			seqMatrices.append(d['MemAcc'])
		else:
			print 'ERROR: The data does not have a feature called MemAcc'
			exit(1)

		if d.has_key('MemTopo'):
	    		seqMatrices.append(d['MemTopo'])
		else:
			print 'ERROR: The data does not have a feature called MemTopo'
			exit(1)

	## Add sequence-template similarity score here. This is used to predict distance matrix from a sequence-template alignment. 
	## this is mainly used for homology modeling
	if config.UseTemplate(modelSpecs):
	  	#print 'Using template similarity score...'
		if not d.has_key('tplSimScore'):
			print 'ERROR: the data has no key tplSimScore, which is needed since you specify to use template information'
			exit(1)

		numTemplateSimFeatures = d['tplSimScore'].shape[1]
		if numTemplateSimFeatures<10 or numTemplateSimFeatures>15:
			print 'WARNING: The number of features for query-template similarity may be incorrect. Please double check it!'

	    	seqMatrices.append( d['tplSimScore'] )

	if returnMode.lower() == 'array':
		seqFeature = np.concatenate(seqMatrices, axis=1)
		return seqFeature
	else:
		return seqMatrices

## data: the real, raw data loaded from files
def ExtractFeaturesNLabels(data, modelSpecs, forTrainValidation=True, returnMode='array'):
        ## each protein has sequential and  pairwise features as input and distance matrix as label
        proteinFeatures=[]
        counter = 0
        
        for d in data:
                oneprotein = dict()
                oneprotein['name'] = d['name'] 
                
                ## convert the primary sequence to a one-hot encoding
                oneHotEncoding = SeqOneHotEncoding(d['sequence'])
        
                ## prepare features for embedding. Currently we may embed a pair of residues or a pair of residue+secondary structure
                if config.EmbeddingUsed(modelSpecs):
                        if modelSpecs['seq2matrixMode'].has_key('Seq+SS'):
                                embedFeature = RowWiseOuterProduct(oneHotEncoding, d['SS3'])
                        else:
                                embedFeature = oneHotEncoding
                        oneprotein['embedFeatures'] = embedFeature
        
                seqFeature = CollectSequentialFeatures(d, modelSpecs, oneHotEncoding, returnMode=returnMode)
                matrixFeature, matrixFeature_nomean = CollectMatrixFeatures(d, modelSpecs, returnMode=returnMode)
        
                oneprotein['sequence'] = d['sequence']
                oneprotein['seqLen'] = len(d['sequence'])
                oneprotein['seqFeatures'] = seqFeature
                oneprotein['matrixFeatures'] = matrixFeature
                oneprotein['matrixFeatures_nomean'] = matrixFeature_nomean

                if forTrainValidation:
                        oneprotein['atomLabelMatrix'] = LabelUtils.CollectLabels(d, modelSpecs)

                ##at this point, finish collecting features and labels for one protein
                proteinFeatures.append(oneprotein)

                counter += 1
                if (counter %500 ==100):
                        print 'assembled features and labels for ', counter, ' proteins.'

        return proteinFeatures


def LoadDistanceFeatures(files=None, modelSpecs=None, forTrainValidation=True):
    	if files is None or len(files)==0:
       		print 'the feature file is empty'
		exit(1)

	data = []
	for f in files:
		with open(f, 'rb') as fh:
			rawD = cPickle.load(fh)
			d = ExtractFeaturesNLabels(rawD, modelSpecs, forTrainValidation)
			data.extend(d)

	return data

def CalcMinMaxMatrixSize(bounds, seqLens):
	## calculate max allowed length for pairwise features
	maxSeqLen = max(seqLens)
	minSeqLen = min(seqLens)

	if bounds is not None:
		boundLens = []
		for box, seqLen in zip(bounds, seqLens):
			if box is None:
				boundLens.append(seqLen)
			else:
				top, left, bottom, right = box
				assert (bottom-top) == (right-left)
				assert (bottom-top) <= seqLen
				boundLens.append(bottom-top)

		## for one protein, when its box is not None, boundLen is at least the box size
		## when its box is None, boundLen is at least its seqLen
		maxMatrixSize = max(boundLens)
		minMatrixSize = min(boundLens)
	else:
		maxMatrixSize = maxSeqLen
		minMatrixSize = minSeqLen

	return minMatrixSize, maxMatrixSize

## use shared memory for some very large arrays
def AssembleOneBatch(data, modelSpecs, forRefState=False, bounds=None, floatType=theano.config.floatX, bUseSharedMemory=False):
	if not data:
		print 'WARNING: the list of data is empty'
		return None

	numSeqs = len(data)
	seqLens = [ d['seqLen'] for d in data ]
	names = [ d['name'] for d in data ]

	## use maxSeqLen and minSeqLen for sequential features
	## we do not crop sequential features at this step since the theano deep model will do so after 1D convolution operation
	maxSeqLen = max(seqLens)
	minSeqLen = min(seqLens)
	#print 'maxSeqLen= ', maxSeqLen, 'minSeqLen= ', minSeqLen

	numSeqFeatures = FeatureUtils.DetermineNumSeqFeatures(data[0]['seqFeatures'])
        X1d = np.zeros(shape=(numSeqs, maxSeqLen, numSeqFeatures), dtype=floatType)

	numMatrixFeatures = FeatureUtils.DetermineNumMatrixFeatures(data[0]['matrixFeatures']) + FeatureUtils.DetermineNumMatrixFeatures(data[0]['matrixFeatures_nomean'])
	## we use maxMatrixSize and minMatrixSize for pairwise features
	## we crop pairwise features at this step to save memory and computational time
	minMatrixSize, maxMatrixSize = CalcMinMaxMatrixSize(bounds, seqLens)

	if bUseSharedMemory:
		shmX2d = SharedNDArray((numSeqs, maxMatrixSize, maxMatrixSize, numMatrixFeatures), dtype=floatType, name='/RaptorX-' + str(os.getppid()) + '-X2d-' + randomString(6) )
		X2d = shmX2d.array
		X2d[:] = 0
	else:
       		X2d = np.zeros(shape=(numSeqs, maxMatrixSize, maxMatrixSize, numMatrixFeatures), dtype=floatType)

	X1dem = None
	if data[0].has_key('embedFeatures'):
		numEmbedFeatures = data[0]['embedFeatures'].shape[1]
            	X1dem = np.zeros(shape=(numSeqs, maxSeqLen, numEmbedFeatures), dtype=floatType)

	## Y shall be a list of 2D or 3D matrices, each for one response
        Y = []
        if data[0].has_key('atomLabelMatrix'):
		for response in modelSpecs['responses']:
			labelName, labelType, _ = ParseResponse(response)
			dataType = np.int16
			if not config.IsDiscreteLabel(labelType):
				dataType = floatType
			rValDims = GetResponseValueDims(response)
			if rValDims == 1:
            			y=np.zeros(shape=(numSeqs, maxMatrixSize, maxMatrixSize), dtype=dataType )
				Y.append(y)
				
			else:
            			y=np.zeros(shape=(numSeqs, maxMatrixSize, maxMatrixSize, rValDims), dtype=dataType )
				Y.append(y)

	## when Y is empty, weight is useless. So When Y is None, weight shall also be None
	weightMatrix = []
	if bool(Y) and config.UseSampleWeight(modelSpecs):
        	weightMatrix = [ np.zeros(shape=(numSeqs, maxMatrixSize, maxMatrixSize), dtype=floatType) ] * len( modelSpecs['responses'] )

	## for mask. we do not used shared ndarray for them since they are small
	M1d = np.zeros(shape=(numSeqs, maxSeqLen - minSeqLen), dtype=np.int8 )
	M2d = np.zeros(shape=(numSeqs, maxMatrixSize - minMatrixSize, maxMatrixSize), dtype=np.int8 )

	if bounds is not None:
		boxes = bounds
	else:
		boxes = [None]* len(data)

        for j, d, box in zip(range(len(data)), data, boxes):
            	seqLen = d['seqLen']

		## posInSeq, posInX and posInY are the starting position of one protein in the final output tensor
		posInSeq = -seqLen

		## here X and Y refer to x-axis and y-axis
		if box is not None:
			top, left, bottom, right = box
			posInX =  -(bottom-top)
			posInY =  -(right-left)
		else:
			posInX = -seqLen
			posInY = -seqLen

		if forRefState:
			## this code needs reexamination, it may not be correct when d['seqFeatures']/d['matrixFeatures'] is represented as a list of arrays instead of a single array
            		X1d[j, posInSeq:, : ] = np.array([ modelSpecs['seqFeatures_expected'] ] * seqLen).reshape((seqLen, -1))

			tmp = [ modelSpecs['matrixFeatures_expected'] ] * (seqLen * seqLen)
			tmp2 = np.array(tmp).reshape((seqLen, seqLen, -1) )
            		tmp3 = np.concatenate((tmp2, d['matrixFeatures_nomean']), axis=2)
			if box is not None:
            			X2d[j, posInX:, posInY:, : ] = tmp3[top:bottom, left:right, ]
			else:
            			X2d[j, posInX:, posInY:, : ] = tmp3
		else:
			if isinstance(d['seqFeatures'], np.ndarray):
            			X1d[j, posInSeq:, : ] = d['seqFeatures']
			else:
				startPos = 0
				for f in d['seqFeatures']:
					if len(f.shape) == 1:
            					X1d[j, posInSeq:, startPos:startPos+1] = f[:, np.newaxis]
						startPos += 1
					elif len(f.shape) == 2:
            					X1d[j, posInSeq:, startPos:startPos+f.shape[1] ] = f
						startPos = startPos + f.shape[1]
					else:
						print 'wrong shape in sequential feature: ', f.shape
						exit(1)

			# add 2D features in matrixFeatures to holder staring from the start position
			# holder is a 3D array and start is the starting position in the 3rd dimension
			def Add2DFeatures(matrixFeatures, holder, start):
				if isinstance(matrixFeatures, np.ndarray):
					features = [ matrixFeatures ]
				else:
					features = matrixFeatures

				startPos = start
				#for f in matrixFeatures:
				for f in features:
					if len(f.shape) == 2:
						endPos = startPos + 1
						if box is None:
                                                	holder[:, :, startPos: endPos ] = f[:, :, np.newaxis]
                                        	else:
                                                	holder[:, :, startPos: endPos ] = f[top:bottom, left:right, np.newaxis]
					elif len(f.shape) == 3:
						endPos = startPos + f.shape[2]
						if box is None:
                                                        holder[:, :, startPos: endPos ] = f
                                                else:
                                                        holder[:, :, startPos: endPos ] = f[top:bottom, left:right, :]
					else:
						print 'wrong shape in matrixFeatures: ', f.shape
						exit(1)
					startPos = endPos

				return endPos

			end = Add2DFeatures(d['matrixFeatures'], X2d[j, posInX:, posInY:, : ], 0)
			Add2DFeatures(d['matrixFeatures_nomean'], X2d[j, posInX:, posInY:, : ], end)

            	M1d[j, posInSeq: ].fill(1)
	    	M2d[j, posInX:, posInY: ].fill(1)

	    	if X1dem is not None:
			## embed feature is always represented as a single array, so the code shall be correct
			if forRefState:
				X1dem[j, posInSeq:, : ] = np.array([ modelSpecs['embedFeatures_expected'] ] * seqLen ).reshape((seqLen, -1))
			else:
				X1dem[j, posInSeq:, : ] = d['embedFeatures']

		
		for y, response in zip(Y, modelSpecs['responses']):
			if box is not None:
                		tmp = d['atomLabelMatrix'][response][top:bottom, left:right]
			else:
                		tmp = d['atomLabelMatrix'][response]
			if len(y.shape) == 3:
                		y[j, posInX:, posInY: ] = tmp
			else:
                		y[j, posInX:, posInY:, ] = tmp

		if bool(weightMatrix):
			if d.has_key('labelWeightMatrix'):
				labelWeightMatrix = d['labelWeightMatrix']
			else:
				labelWeightMatrix = LabelUtils.CalcLabelWeightMatrix(d['atomLabelMatrix'], modelSpecs, floatType=floatType)

		for w, response in zip(weightMatrix, modelSpecs['responses']):
			if box is not None:
                		w[j, posInX:, posInY: ] = labelWeightMatrix[response][top:bottom, left:right]
			else:
                		w[j, posInX:, posInY: ] = labelWeightMatrix[response]
	 
	if bUseSharedMemory:
		onebatch = [X1d, shmX2d, M1d, M2d] 
	else: 
		onebatch = [X1d, X2d, M1d, M2d] 

	if X1dem is not None:
	    	onebatch.append(X1dem)

	onebatch.extend(Y)
	onebatch.extend(weightMatrix)

	return onebatch, names

## this function add labels and their weights to an already-assembled batch
## batch is generated by AssembleOneBatch
def AddLabel2OneBatch(names, batch, modelSpecs, sharedLabelPool, sharedLabelWeightPool, floatType=theano.config.floatX):

	numSeqs = len(names)
	for name in names:
		if (not sharedLabelPool.has_key(name)) or (not sharedLabelWeightPool.has_key(name)):
			print 'the label or label weight matrix does not exist for protein ', name
			exit(1)

	seqLens = [ sharedLabelWeightPool[name].shape[0] for name in names ]

	## get the boundingbox for this batch
	if not config.TrainByRefLoss(modelSpecs):
                box = batch[-1]
        else:
                box = batch[-2]

	top, left, bottom, right = box
	assert bottom-top == right-left
	boxsize = bottom - top

	if boxsize < max(seqLens) and numSeqs>1:
		## make sure that there is only one protein in this batch
		print 'ERROR: when one batch has a large protein, it can only have one protein'
		exit(1)

	## we crop pairwise labels at this step to save memory and computational time
	maxMatrixSize = min(boxsize, max(seqLens) )

	## Y shall be a list of 2D or 3D matrices, each for one response
        Y = []
	for response in modelSpecs['responses']:
		labelName, labelType, _ = ParseResponse(response)
		dataType = np.int16
		if not config.IsDiscreteLabel(labelType):
			dataType = floatType
		rValDims = GetResponseValueDims(response)
		if rValDims == 1:
            		y=np.zeros(shape=(numSeqs, maxMatrixSize, maxMatrixSize), dtype=dataType )
			Y.append(y)
				
		else:
            		y=np.zeros(shape=(numSeqs, maxMatrixSize, maxMatrixSize, rValDims), dtype=dataType )
			Y.append(y)

	## when Y is empty, weight is useless. So When Y is empty, weight shall also be empty
	weightMatrix = []
	if bool(Y) and config.UseSampleWeight(modelSpecs):
        	weightMatrix = [ np.zeros(shape=(numSeqs, maxMatrixSize, maxMatrixSize), dtype=floatType) ] * len( modelSpecs['responses'] )

        for j, name, seqLen in zip(range(len(names)), names, seqLens):

		## we align all matrices in the bottom/right corner
		## posInX and posInY are the starting position of one protein in the final output tensor
		## here X and Y refer to x-axis and y-axis
		posInX =  -min(boxsize, seqLen)
		posInY =  -min(boxsize, seqLen)

		for y, response in zip(Y, modelSpecs['responses']):

			if boxsize < seqLen:
                		tmp = sharedLabelPool[name][response][top:bottom, left:right]
			else:
                		tmp = sharedLabelPool[name][response]
			if len(y.shape) == 3:
                		y[j, posInX:, posInY: ] = tmp
			else:
                		y[j, posInX:, posInY:, ] = tmp

		labelWeightMatrix = sharedLabelWeightPool[name]
		for w, response in zip(weightMatrix, modelSpecs['responses']):
			if boxsize <seqLen:
                		w[j, posInX:, posInY: ] = labelWeightMatrix[response][top:bottom, left:right]
			else:
                		w[j, posInX:, posInY: ] = labelWeightMatrix[response]

	## the input batch contains bounding box
	tail = 1

	## check to see if the input batch contains one flag for RefState
	if config.TrainByRefLoss(modelSpecs):
		tail += 1

	newbatch = batch[: -tail ]
	newbatch.extend(Y)
	newbatch.extend(weightMatrix)
	newbatch.extend(batch[-tail: ])

	return newbatch

## this function prepares one batch of data for training, validation and prediction
## data is a list of proteins with seq and pairwise features and possibly labels
## when forRefState is set to True, then generate a batch of data for the calculation of reference probability

## bounds is a list of bounding boxes, each for one protein in data. Each box is a list 4 entries: top, left, bottom, right
## if bounds is None, then do not crop any protein data.
## if one entry in bounds is None, the corresponding protein is small and thus, does not need crop
## for validation/prediction data, bounds shall be set to None (or very large) since we do not crop validation data
## for train data, if one protein has size larger than allowed, its boundingbox shall not be None

##Note that right now each bounding box shall be a square instead of rectangle, although it does not have to be on the diagnonal line of the protein contact/distance matrix

## this function returns one minibatch and its corresponding proteins names. The minibatch is a list of tensors, each tensor represents one type of features/labels
## for sequential features, the output tensor is not cropped, regardless of bounds
## for pairwise features/labels, the output tenosr is cropped if bounds is not None

"""
def AssembleOneBatch(data, modelSpecs, forRefState=False, bounds=None, floatType=theano.config.floatX):
"""

## exclude those proteins with length > maxLen. Use this function for validation proteins only
def FilterByLength(data, maxLen=800):
	numOriginals = len(data)
	newData = [ d for d in data if d['seqLen'] < maxLen ]
	numNews = len(newData)

	return newData, numOriginals - numNews
	
##this function only splits a list of raw data into batch, but does not actually generate a minibatch
##You shall use AssembleOneBatch to generate a minibatch from the raw data
##one minibatch shall be generated right before calling the train, validate and predict functions
def SplitData2Batches(data=None, numDataPoints=1000000, modelSpecs=None):

    	if data is None:
        	print 'The raw data for batch assignment is empty!'
		exit(1)

    	if numDataPoints < 10:
		print 'Please specify the number of data points in a minibatch'
		exit(1)

    	## sort proteins by length from large to small
    	data.sort(key=lambda x: x['seqLen'], reverse=True)

	##seqDataset stores the resultant data
    	batches = []

    	i = 0
    	while i < len(data):

        	currentSeqLen = data[i]['seqLen']
		numSeqs = min( len(data) - i, max(1, numDataPoints/np.square(currentSeqLen) ) )
		#print 'This batch contains ', numSeqs, ' sequences'

		oneBatch = data[i : i+numSeqs]
		batches.append(oneBatch)

		i += numSeqs

    	return batches

## mode=0, the traditional method, i.e., SplitData2Batches, in which each protein has the same weight

## mode=1, each protein has weight 1
## mode=2, a large protein has weight L*L/matrixSize where L is its length and a small protein has weight 1
## mode=3, a protein has weight L*L/matrixSize where L is its length
## here one protein is large if L >= sqrt(matrixSize), otherwise small
## when mode in [1, 2, 3], modelSpecs['minibatchSize'] shall be interpreted as matrixSize measured by the number of residue pairs
## and modelSpecs['maxbatchSize'] shall be interpreted as the batch size in terms of the number of residue pairs
def NewSplitData2Batches(data=None, modelSpecs=None, numBatchBound=15000):
	if not modelSpecs.has('batchMode'):
		mode = 0
	else:
		mode = modelSpecs['batchMode']

	if mode == 0:
		return SplitData2Batches(data, numDataPoints=modelSpecs['minibatchSize'], modelSpecs=modelSpecs)

    	if data is None:
        	print 'The raw data for batch assignment is empty!'
		exit(1)

	matrixSize = modelSpecs['minibatchSize']

	if mode == 1:
		w = [1.] * len(data)
	elif mode == 2:
		w = [ min(1., np.square( d['seqLen'] )*1./matrixSize ) for d in data ] 
	elif mode == 3:
		w = [ np.square( d['seqLen'] )*1./matrixSize for d in data ]
	else:
		print 'unsupported batch assignment mode'
		exit(1) 

	numProteinsPerBatch = modelSpecs['maxbatchSize'] / matrixSize
	numBatches = min(numBatchBound, np.sum(w) / numProteinsPerBatch )

	## convert w to probability
	w = w / np.sum(w)
	
	indices = np.arange(len(data) )
	sampling = np.random.choice(indices, size=(numBatches, numProteinsPerBatch), replace=True, p=w)

	batches = []
	for row in sampling:
		oneBatch = [ data[i] for i in row ]
		batches.append(oneBatch)

	return batches


## this function may need correction. here batches is a list of real data
def CalcAvgWeightPerBatch(batches, modelSpecs):
	if not modelSpecs['UseSampleWeight']:
		return None

	numResponses = len(modelSpecs['responses'])
	allWeights = []

	for b in batches:
		oneBatchWeight = []
		for wMatrix in b[-numResponses: ]:
			bounds = SampleBoundingBox( (wMatrix.shape[1], wMatrix.shape[2]), modelSpecs['maxbatchSize'] )
			new_wMatrix = wMatrix[:,  bounds[0]:bounds[2], bounds[1]:bounds[3] ]
			wSum = np.sum(new_wMatrix) 
			oneBatchWeight.append(wSum)

		allWeights.append(oneBatchWeight)

	avgWeights = np.average(allWeights, axis=0)

	modelSpecs['batchWeightBase'] = np.array(avgWeights).astype(theano.config.floatX)

	maxWeights = np.amax(allWeights, axis=0)
	minWeights = np.amin(allWeights, axis=0)

	## reutrn the maximum deviation
	return maxWeights/avgWeights, minWeights/avgWeights
	

## proteinInfo is a single location or a list of locations for real data
## Each entry is a Python dict() with keys: name, seqLen, one savefolder for features and two savefolders for ground truth
def LoadRealData(proteinInfo, modelSpecs, loadFeature=True, loadLabel=True, returnMode='array'):

	if not isinstance(proteinInfo, (tuple, list, np.ndarray)):
		locations = [ proteinInfo]
	else:
		locations = proteinInfo

	labelNames = []
	for response in modelSpecs['responses']:
		labelName, _, _ = ParseResponse(response)
		labelNames.append(labelName)

	oriNeeded = bool(  set(labelNames).intersection( set(config.allDistLabelNames) ) ) and loadLabel

	## Cb-Cb distMatrix is needed to discretize orientation matrix
	distNeeded = bool(  set(labelNames).intersection( set(config.allDistLabelNames) ) ) and loadLabel or oriNeeded

        proteins = []
        for protein in locations:
                rawData = dict()
                rawData['name'] = protein['name']
                rawData['sequence'] = protein['sequence']
		rawData['length'] = len(rawData['sequence'])

		## load data
                if loadFeature:
                        feature = FeatureUtils.LoadFeaturePKL(protein['name'], location=protein['featureLocation'], modelSpecs=modelSpecs)
                        rawData.update(feature)

		if loadFeature and config.UseTemplate(modelSpecs):
			templates = protein['templates']
			if len(templates) < 1:
				#print 'ERROR: template is requested, but none available for protein: ', protein['name']
				#exit(1)
				feature = AlignmentUtils.GenNullAlignmentFeatures( len(protein['sequence']) )
			else:
				## randomly select one template
				template = random.choice(templates)
				feature = AlignmentUtils.GenerateAlignmentFeatures2( (protein['name'], template), queryData=rawData, aliDir=protein['aliLocation'], tplDir=protein['tplLocation'], modelSpecs=modelSpecs)

			## feature is a Dict() with keys: tplSimScore, tplDistMatrix and tplOriMatrix. Their values have types array, dict, and dict
			rawData.update(feature)

		if distNeeded or oriNeeded:
                        truth = NativeUtils.LoadGroundTruth(protein['name'], location=protein['nativeLocation'])
			NativeUtils.AddGroundTruth(rawData, truth)

		"""
                if distNeeded:
                        atomDistMatrix = DistanceUtils.LoadNativeDistMatrix(protein['name'], location=protein['distLocation'])
                        DistanceUtils.AddDistMatrix(rawData, atomDistMatrix)

                if oriNeeded:
                        atomOriMatrix = OrientationUtils.LoadNativeOriMatrix(protein['name'], location=protein['oriLocation'])
                        OrientationUtils.AddOrientationMatrix(rawData, atomOriMatrix)
		"""

		## extract information needed for this specific deep model
        	if loadFeature:
                	result = ExtractFeaturesNLabels([rawData], modelSpecs, forTrainValidation=loadLabel, returnMode=returnMode)
			result = result[0]
		elif loadLabel:
			result = dict()
			result['name'] = protein['name']
			result['atomLabelMatrix'] = LabelUtils.CollectLabels(rawData, modelSpecs)
		else:
			print 'ERROR: in LoadRealData(), you shall load input features and/or labels'
			exit(1)

                proteins.append(result)

	if not isinstance(proteinInfo, (tuple, list, np.ndarray)):
		return proteins[0]

        return proteins

