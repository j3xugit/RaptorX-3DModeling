import os
import sys
import cPickle
import numpy as np

import config
from config import Response2LabelName, Response2LabelType, ParseResponse, GetResponseProbDims, GetResponseValueDims

import DataProcessor
import DistanceUtils
import OrientationUtils
import RangeNWeight

from utilsNoT import SampleBoundingBox

import theano

## this function converts the real ground truth (distance and orientation) matrices into label matrices
## d is real data for one protein. it is a dict() with keys: name, atomDistMatrix and atomOrientationMatrix
def CollectLabels(d, modelSpecs):
	atomLabelMatrix = dict()

	for response in modelSpecs['responses']:
		labelName, labelType, subType = ParseResponse(response)
		if labelName in config.allDistLabelNames:
			if not d.has_key('atomDistMatrix'):
				print 'ERROR: atomic distance matrix is needed for ', d['name'], ' to be used for training and/or validation'
				exit(1)

			atomDistMatrix = d['atomDistMatrix']

			if not atomDistMatrix.has_key(labelName):
                               	print 'In the raw data, ', d['name'], ' does not have distance matrix for ', labelName
                               	exit(1)
			distm = atomDistMatrix[labelName]

			if labelName in [ 'HB', 'Beta']:
				## no need to discretize since they are already binary matrices
	    			atomLabelMatrix[response] = distm

			elif labelType.startswith('Discrete'):
				## process the other atom pairs such as Cb-Cb, Ca-Ca
	    			labelMatrix, _, _  = DistanceUtils.DiscretizeDistMatrix(distm, config.distCutoffs[subType], invalidDistanceSeparated = (subType.endswith('Plus') or subType.endswith('Minus')) )
	    			atomLabelMatrix[response] = labelMatrix

			elif labelType.startswith('LogNormal'):
	    			atomLabelMatrix[response] = DistanceUtils.LogDistMatrix(distm)
			elif labelType.startswith('Normal'):
	    			atomLabelMatrix[response] = distm
			else:
				print 'ERROR: unsupported labelName or labelType in response: ', response
				exit(1)

		elif labelName in config.allOrientationNames:
			if (not d.has_key('atomOrientationMatrix') or (not d.has_key('atomDistMatrix')) ):
				print 'ERROR: both atomic orientation and distance matrices are needed for ', d['name'], ' as training and/or validation data'
				exit(1)

			atomDistMatrix = d['atomDistMatrix']
			if not atomDistMatrix.has_key('CbCb'):
				print 'ERROR: Cb-Cb distance matrix is needed to discretize orientation matrix for protein ', d['name']
				exit(1)

			atomOrientationMatrix = d['atomOrientationMatrix']
			if not atomOrientationMatrix.has_key(labelName):
                               	print 'ERROR: In the raw data, ', d['name'], ' does not have orientation matrix for ', labelName
                               	exit(1)

			orim = atomOrientationMatrix[labelName]
			if labelType.startswith('Discrete'):
				if labelName in config.allDihedralNames:
					bins = config.dihedralCutoffs[subType]
				else:
					bins = config.angleCutoffs[subType]
				labelMatrix, _, _ = OrientationUtils.DiscretizeOrientationMatrix(orim, bins=bins, distThreshold4Orientation=modelSpecs['distThreshold4Orientation'], distMatrix=atomDistMatrix['CbCb'], invalidEntrySeparated= (subType.endswith('Plus') or subType.endswith('Minus')) )
				atomLabelMatrix[response] = labelMatrix
			else:
				print 'ERROR: unsupported labelType in response: ', response
				exit(1)
		else:
			print 'ERROR: unsupported labelName in response: ', response
			exit(1)

	return atomLabelMatrix

## this function reads protein names and extracts label matrices from a list of files that may save both input features and native distance/orientation matrices
def LoadLabelMatrices(files=None, modelSpecs=None):
    	if files is None or len(files)==0:
       		print 'ERROR: the feature file is empty'
		exit(1)

	data = []
	for file in files:
		with open(file, 'rb') as fh:
			data.append(cPickle.load(fh))

    	proteinFeatures=[]
    	counter = 0

   	for d in data:
        	oneprotein = dict()
        	oneprotein['name'] = d['name']
		oneprotein['sequence'] = d['sequence']
		oneprotein['seqLen'] = len(d['sequence'])
		oneprotein['atomLabelMatrix'] = CollectLabels(d, modelSpecs)

        	proteinFeatures.append(oneprotein)

		counter += 1
		if (counter %500 ==1):
            		print 'assembled features and labels for ', counter, ' proteins.'

    	return proteinFeatures


## data is a list of proteins. Each protein is a dict with at least 'atomLabelMatrix' as the key
def CalcLabelDistribution(data, modelSpecs):
	## collect all discrete label matrices
	allLabelMatrices = dict()
	for response in modelSpecs['responses']:
		labelType = Response2LabelType(response)
		if labelType.startswith('LogNormal') or labelType.startswith('Normal'):
			continue

		allLabelMatrices[response] = [ d['atomLabelMatrix'][response] for d in data ]


	## calculate the discrete label distribution
	allRefProbs = dict()
	for response in modelSpecs['responses']:
		labelName, labelType, subType = config.ParseResponse(response)
		if labelType.startswith('LogNormal') or labelType.startswith('Normal'):
			allRefProbs[response] = np.array([1.] * numRanges).reshape((-1, 1)).astype(np.float32)
			continue

		if modelSpecs.has_key('UseBoundingBox4RefProbs') and (modelSpecs['UseBoundingBox4RefProbs'] is True):
			## here we sample a sub label matrix using BoundingBox to account for the real training scenario
			newLabelMatrices = []
			for lMatrix in allLabelMatrices[response]:
				bounds = SampleBoundingBox( (lMatrix.shape[0], lMatrix.shape[1]),  modelSpecs['maxbatchSize'] )
				new_lMatrix = lMatrix[ bounds[0]:bounds[2], bounds[1]:bounds[3] ].astype(np.int32)
				newLabelMatrices.append(new_lMatrix)
			if labelName in config.allOrientationNames:
				allRefProbs[response] = OrientationUtils.CalcLabelProb(data = newLabelMatrices, numLabels = GetResponseProbDims(response), numRanges=RangeNWeight.GetNumRanges(modelSpecs))
			else:
				allRefProbs[response] = DistanceUtils.CalcLabelProb(data = newLabelMatrices, numLabels = GetResponseProbDims(response), numRanges=RangeNWeight.GetNumRanges(modelSpecs))
		else:
			if labelName in config.allOrientationNames:
				allRefProbs[response] = OrientationUtils.CalcLabelProb(data=[ m.astype(np.int32) for m in allLabelMatrices[response] ], numLabels = GetResponseProbDims(response), numRanges=RangeNWeight.GetNumRanges(modelSpecs))
			else:
				allRefProbs[response] = DistanceUtils.CalcLabelProb(data=[ m.astype(np.int32) for m in allLabelMatrices[response] ], numLabels = GetResponseProbDims(response), numRanges=RangeNWeight.GetNumRanges(modelSpecs) )

        modelSpecs['labelDistributions'] = allRefProbs
	return allRefProbs


def CalcLabelWeight(modelSpecs):
	print 'Calculating label weight ...'

	numRanges = RangeNWeight.GetNumRanges(modelSpecs)

	RangeNWeight.SetWeight4Range(modelSpecs)
        #print 'weight for range: ', modelSpecs['weight4range']

	RangeNWeight.SetWeight43C2C(modelSpecs)
       	#print 'LRbias= ', modelSpecs['LRbias']
	#print 'weight43C= ', modelSpecs['weight4Discrete3C']


	allRefProbs = modelSpecs['labelDistributions']
	##for discrete labels, we calculate their weights by inferring from the weight intialized to 3 bins: 0-8, 8-15 and >15 or -1, which makes inference easier
        modelSpecs['weight4labels'] = dict()

	for response in modelSpecs['responses']:
		labelName, labelType, subType = config.ParseResponse(response)
		numLabels = GetResponseProbDims(response)

		if config.IsContinuousLabel(labelType):
			## just need to assign range weight for continuous response
			modelSpecs['weight4labels'][response] = modelSpecs['weight4continuous']
			continue

		if not config.IsDiscreteLabel(labelType):
			print 'ERROR: unsupported response in CalcLabelWeight: ', response
			exit(1)

		if  labelName in config.allOrientationNames or config.NoWeight4Label(modelSpecs):
			modelSpecs['weight4labels'][response] = np.multiply( np.ones((numRanges, numLabels), dtype=np.float32), modelSpecs['weight4range'])

		elif labelName in ['HB', 'Beta']:
			## if the response is for HB and Beta-Pairing
			if subType.startswith('2C'):
				modelSpecs['weight4labels'][response] = modelSpecs['weight4' + response]
			else:
				print 'ERROR: unsupported label subtype in CalcLabelWeight: ', response
				exit(1)

		elif labelName in config.allAtomPairNames:
			## calculate label weight for atom pairs Cb-Cb, Ca-Ca, Cg-Cg, CaCg, and NO
			if subType.startswith('2C'):
				print 'ERROR: 2C is not supported for contact/distance prediction any more'
				exit(1)
			elif subType.startswith('3C'):
				## if 3C is used for the response
				modelSpecs['weight4labels'][response] = modelSpecs['weight4Discrete3C']
			else:
				modelSpecs['weight4labels'][response] = DistanceUtils.CalcLabelWeight(modelSpecs['weight4Discrete3C'], allRefProbs[response], config.distCutoffs[subType] )

		else:
			print 'ERROR: unsupported label name in CalcLabelWeight: ', response
			exit(1)

		## set the weight of the label for the invalid entry (distance or orientation) to 0
		if subType.endswith('Minus'):
			modelSpecs['weight4labels'][response][:,-1] = 0

	"""
	## for log
	for response in modelSpecs['responses']:
		print 'weight4labels for response: ', response
		print modelSpecs['weight4labels'][response]
	"""

	return modelSpecs['weight4labels']

##this function calculates the label distribution of the training proteins and then label weight for long-, medium-, short- and near-range labels
## to assign weight to a specific label matrix, please use another function CalcLabelWeightMatrix()
## data is a list of dict(), which has at least one key 'atomLabelMatrix' 
## the weight factor for a continuous distance label (i.e., regression) is a 4*1 matrix
## the weight factor for a discrete distance label (i.e., classificaton) is a 4*numLabels matrix
def CalcLabelDistributionAndWeight(data=None, modelSpecs=None):
	CalcLabelDistribution(data, modelSpecs)
	CalcLabelWeight(modelSpecs)
	
## do not mix this function with the previous function
## trainData is a list of proteins, each represented as dict()
def CalcLabelDistributionNWeightBySampling(trainMetaData, modelSpecs):
        trainDataLocation = DataProcessor.SampleProteinInfo(trainMetaData, numSamples=10000)

        ## only load ground truth but not input features to save memory and speed up
        labelData = []
	for loc in trainDataLocation:
        	p = DataProcessor.LoadRealData(loc, modelSpecs, loadFeature=False)
                labelData.append(p)

        CalcLabelDistributionAndWeight(labelData, modelSpecs)


## this function assigns weight to each element of a specific label matrix
## same label may have different weights depending on if the residue pair is in near-range, short-range, medium-range or long-range.
## labelMatrices is a dictionary and has an entry for each response. 
## This function returns a dictionary object for labelWeightMatrix
def CalcLabelWeightMatrix(LabelMatrix=None, modelSpecs=None, floatType=theano.config.floatX):
    	if LabelMatrix is None:
		return None

	RangeBoundaries = RangeNWeight.GetRangeBoundaries(numRanges=RangeNWeight.GetNumRanges(modelSpecs))

	shape = LabelMatrix.values()[0].shape
	a = np.mgrid[0:shape[0], 0:shape[1] ]
	seqSeparation = abs(a[0] - a[1])
	masks = []
	for i, bound in zip(range(len(RangeBoundaries)), RangeBoundaries):
		if i==0:
			mask = (seqSeparation>=bound).astype(np.int16)
		else:
			mask = ( (seqSeparation>=bound) * (seqSeparation<RangeBoundaries[i-1]) ).astype(np.int16)
		masks.append(mask)
	
	for response in modelSpecs['responses']:
		if not modelSpecs['weight4labels'].has_key(response):
			print 'ERROR: Cannot find the weight factor for response ', response
			exit(1)

	##the below procedure is not very effective. We shall improve it later.
	labelWeightMatrices = dict()
	for response in modelSpecs['responses']:
		labelName, labelType, subType = config.ParseResponse(response)

		## wMatrix is a matrix with dimension config.numRanges * numLabels, where numRanges corresponds to ER, LR, MR, SR, and NR, respectively
		wMatrix = (modelSpecs['weight4labels'][response]).astype(np.float32)
		assert wMatrix.shape[0] == len(RangeBoundaries)

		if config.IsContinuousLabel(labelType):
			low, high = config.GetLabelMinMaxValues(labelName)
			## if the label is real value, then for each range, there is only a single weight for all values
			tmpWeightMatrices = []
			for w in wMatrix:
				## the below two sentences may be incorrect. Need further examination.
				M0s = np.zeros_like(LabelMatrix[response], dtype=np.int16)
				tmp = w[ M0s ]
				if labelName in config.allOrientationNames:
					## set the weight for the invalid orientation entry to 0. 
					np.putmask(tmp, LabelMatrix[response] > high, 0 )
					np.putmask(tmp, LabelMatrix[response] < low, 0 )
				else:
					## set the weight for the invalid distance entry to 0. An invalid entry in the label matrix is indicated by a negative value,e.g., -1
					np.putmask(tmp, LabelMatrix[response] < low, 0 )
				tmpWeightMatrices.append(tmp)
		else:
			tmpWeightMatrices = [ w[LabelMatrix[response] ] for w in wMatrix ]

		labelWeightMatrices[response] = sum( [ m * w for m, w in zip(masks, tmpWeightMatrices) ] ).astype(floatType)

    	return labelWeightMatrices

"""
## need revision
##hopefully this is a better implementation of CalcLabWeightMatrix, but needs test for efficiency and correctness
def CalcLabelWeightMatrix2(LabelMatrix=None, modelSpecs=None):
        if LabelMatrix is None:
                return None

        labelType = modelSpecs['distLabelType']
        if not modelSpecs.has_key('weight4' + labelType):
                print 'Cannot find the label weight for the label type ', labelType, '. Please make sure that it has been generated already.'
                exit(-1)

        ##wMatrices is a dictionary. Each item is a matrix of size 4*numLabels for each atomPairType.
        wMatrices = modelSpecs['weight4' + labelType]

        ##the below procedure is not very effective. We shall improve it later.
        labelWeightMatrices = dict()
        for apt in modelSpecs['atomPairTypes']:
		size = LabelMatrix[apt].shape

		## for long-range weight
		labelWeightMatrices[apt] = np.choose(LabelMatrix[apt], wMatrices[apt][0])

		## set the diagonal line to 0
		np.fill_diagonal(labelWeightMatrices[apt], 0)
		
		##for medium-range weight
		for offset in np.arange(23, 11, -1):
			i = np.arange(0, size[0]-offset)
			j = i + offset
			labelWeightMatrices[apt][ i, j ] = np.choose(LabelMatrix[apt].diagonal(offset), wMatrices[apt][1] )
		for offset in np.arange(-23, -11, 1):
			i = np.arange(-offset, size[0])
			j = i + offset
			labelWeightMatrices[apt][ i, j ] = np.choose(LabelMatrix[apt].diagonal(offset), wMatrices[apt][1] )

		##for short-range weight
		for offset in np.arange(11, 5, -1):
			i = np.arange(0, size[0]-offset)
			j = i + offset
			labelWeightMatrices[apt][ i, j ] = np.choose(LabelMatrix[apt].diagonal(offset), wMatrices[apt][2] )
		for offset in np.arange(-11,-5, 1):
			i = np.arange(-offset, size[0])
			j = i + offset
			labelWeightMatrices[apt][ i, j ] = np.choose(LabelMatrix[apt].diagonal(offset), wMatrices[apt][2] )

		if modelSpecs['rangeMode'] != 'All':
			for offset in np.arange(5, 1, -1):
				np.fill_diagonal(labelWeightMatrics[apt][0:size[0]-offset, offset:size[1]], 0)
			for offset in np.arange(1, 5, 1):
				np.fill_diagonal(labelWeightMatrics[apt][offset:size[0], 0:size[1]-offset ], 0)
				
			continue

		##for near-range weight
		for offset in np.arange(6, 1, -1):
			i = np.arange(0, size[0]-offset)
			j = i + offset
			labelWeightMatrices[apt][ i, j ] = np.choose(LabelMatrix[apt].diagonal(offset), wMatrices[apt][3] )
		for offset in np.arange(-6,-1, 1):
			i = np.arange(-offset, size[0])
			j = i + offset
			labelWeightMatrices[apt][ i, j ] = np.choose(LabelMatrix[apt].diagonal(offset), wMatrices[apt][3] )

        return labelWeightMatrices

"""

## Get one-hot encoding of labelMatrix. LabelMatrix has dimension L*L and the resultant matrix is binary and has dimension L*L*numLabels
def GetOneHotLabelMatrix(labelMatrix, numLabels):
	nRows, nCols = labelMatrix.shape
	result = np.zeros((nRows, nCols, numLabels), dtype=np.int16)
	mg = np.mgrid[0:nRows, 0:nCols]
	result[ mg[0], mg[1], labelMatrix] =1

	return result
