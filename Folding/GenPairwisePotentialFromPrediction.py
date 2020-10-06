import numpy as np
import scipy
import cPickle
import sys
import os
import getopt
import copy
import random

from DL4DistancePrediction4 import config
from DL4DistancePrediction4 import DistanceUtils
from DL4DistancePrediction4 import RangeNWeight

def Usage():
    	print 'python GenPairwisePotentialFromPrediction.py [-a labelNames] [-w weightScheme] [-r refType | -f refFile] [-o] [-l minPotential] [-u maxPotential] [ -s savefile ] predidctedPairwiseDistribution_PKL'
	print '  This scripts convert predicted distance/orientation distribution to potential and save it in PKL format'
    	print '      The input is a PKL file for predicted distance/orientation probability, e.g. target.predictedDistMatrix.pkl or target.mergedDistMatrix.pkl where target is a protein name'
    	print '      This file is a tuple of 6 or 7 items: name, primary sequence, predDistProb, predContactMatrix, labelWeight, labelDistribution and refDistProb (optional)'
	print '		labelWeight: the weight of labels used in training and labelDistribution is the background probability of labels in training data'
	print '		refDistProb: the distance distribution predicted by the DL model from expected input features'
	print ''
	print '  -a: labelNames, e.g., AllOri, CbCb+CaCa+NO+AllOri, CbCb, CaCa, CaCa+CbCb, AllAP+AllOri (default)'
	print ''
	print '  -w: weight scheme for distance/orientation potential (default 3): 0 (no weight), 1 (use weight for distance), 2(use weight for orientation) and 3(for both)'
	print '		for distance, weight = 1-Prob(disorder); for orientation, weight = Prob(dist<20) where P() is predicted probability'

	print '  -r: reference state (and some parameters) for distance potential, including DFIRE (default), DOPE and SimuRW (case insensitive)'
	print '		DFIRE+18+1.61: 18 (default) is the dist cutoff and 1.61 (default) is the exponent'
	print '			if the exponent is a value>10, then a random value between 1.57 and 1.63 will be used'
	print '		DOPE+20+1.0:  20 is the dist cutoff and 1.0 is the scale factor for radius of gyration'
	print '		SimuRW+20: 20 is the dist cutoff'
	print '  -f: file needed for reference state SimuRW'
	print '  -o: do not use reference for orientation (default Use)'
	print '  -s: save file (shall end with .pkl); when empty, a PKL file will be automatically generated under current work directory with name derived from methods and parameter setting'
	print '		the resultant file is a tuple containing name, sequence, pairPotential, cutoffs and validProbability'
	print '		pairPotential is a dict(), cutoffs defines the discretization boundary and validProbability may be used as weight for corresponding potential'
	print ''
	print '  -l, -u: min and max potential for one pair of atoms, default values: -30 and 30, respectively. They are only used for potential check and do not impact the real potential values'

#check value of one potential matrix
def CheckPotentialValues(m, minPotential=-40., maxPotential=40):
	assert m is not None

	mWithoutDiag = m[~np.eye(m.shape[0],dtype=bool)].reshape(m.shape[0],-1)
	if not np.isfinite(mWithoutDiag).all():
		print 'ERROR: the potential matrix may contain an NaN or infinite value'
		return False

	numLabels = m.shape[2]
	for i in np.arange(numLabels):
		## here we only check the potential values of two residues with sequence separation >=k
		potential = np.triu(m[:, :, i], k=4)
		if np.amax(potential ) > maxPotential:
                        print 'WARNING: possibly too big potential value: ', np.amax(potential), ' for label ', i
                if np.amin(potential ) < minPotential:
                        print 'WARNING: possibly too small potential value: ', np.amin(potential), ' for label ', i
	return True

## this function is obsolete
## predDistMatrix is the predicted distance distribution by a deep learning model.
## labelWeight is the weight of the label used in training a deep learning model
## labelDistribution is the label distribution in the training data
## when labelWeight != 1, the prob values in predDistMatrix shall be corrected before being converted into distance potential
def FixPredictedDistProb(predDistMatrix, labelWeight, labelDistribution):
	newPredDistMatrix = dict()
	for response in predDistMatrix.keys():
    		fixedProb = DistanceUtils.FixDistProb( predDistMatrix[response], labelWeight[response], labelDistribution[response])
		newPredDistMatrix[response] = fixedProb

    	return newPredDistMatrix

def CalcOrientationPotential(predOriMatrix, useRef=True, useWeight=True, labelWeight=None, labelDistribution=None, minPotential=-30, maxPotential=30):
	potentials = dict()
	## validProbs save the Prob(d<20) for all atom/residue pairs
	validProbs = dict()
	for response, predProb in predOriMatrix.iteritems():
		labelName, labelType, subType = config.ParseResponse(response)
                if labelName not in config.allOrientationNames:
                        #print 'WARNING: unsupported response for orientation potential: ', response
                        continue
                if not config.IsDiscreteLabel(labelType):
                        print 'WARNING: the orientation label is not discrete: ', response
                        continue
		numLabels = config.GetResponseProbDims(response)
		if subType.endswith('Plus') or subType.endswith('Minus'):
			largestValidLabel = numLabels -2
		else:
			largestValidLabel = numLabels -1

		probOfValid = predProb[:, :, :largestValidLabel]
		#potential = np.zeros_like(probOfValid)
		potential = -np.log(probOfValid)

		if useRef:
			refOfValid = labelDistribution[response][:,:largestValidLabel ]
			refPot = -np.log(refOfValid)
			for i in range(predProb.shape[0]):
				for j in range(predProb.shape[1]):
					if i==j:
						continue
					offset = abs(i-j)
					rangeIndex = RangeNWeight.GetRangeIndex(offset, numRanges=refOfValid.shape[0])
                        		if rangeIndex < 0:
                        			continue
					potential[i, j] -= refPot[rangeIndex]

		## shift potential by its mean
		potential -= np.mean(potential, axis=2, keepdims=True)

		validProb = np.sum(probOfValid, axis=2)

		## multiply the potential by the probability of distance<20A
		if useWeight:
			potential *= validProb[:, :, np.newaxis]
			
		"""
		x = np.sum(probOfValid, axis=2)
		np.fill_diagonal(x, 0)
		np.fill_diagonal(x[1:], 0)
		np.fill_diagonal(x[:,1:], 0)
		ind = np.unravel_index(np.argsort(-x, axis=None), x.shape)

		ratio = min(topRatio, predProb.shape[0]-3)
		topK = np.int32(ratio * predProb.shape[0])
		"""
		"""
		for i, j in zip(ind[0][:topK], ind[1][:topK]):
			offset = abs(i-j)
			rangeIndex = RangeNWeight.GetRangeIndex(offset, numRanges=refOfValid.shape[0])
                        if rangeIndex < 0:
                        	continue
			potential[i, j] = -np.log(probOfValid[i, j])
			if useRef:
				potential[i, j] -= refPot[rangeIndex]
			potential[i, j] -= np.mean(potential[i, j])
			if weighted:
				potential[i, j] *= x[i, j]
		"""

		CheckPotentialValues(m=potential)
		potentials[response] = potential.astype(np.float32)
		validProbs[response] = validProb.astype(np.float32)

	return potentials, validProbs

def CalcPotentialByDFIRE(predDistMatrix, alpha=1.61, largestDistance=18, useWeight=False, minPotential=-30, maxPotential=30):
	potentials = dict()

	## validProbs saves the prob of one atom/residue pair likely have valid coordinates
	validProbs = dict()
	for response in predDistMatrix.keys():
		labelName, labelType, subType = config.ParseResponse(response)
		if labelName not in config.allAtomPairNames:
			#print 'WARNING: unsupported response for DFIRE potential: ', response
			continue
		if not config.IsDiscreteLabel(labelType):
                        print 'WARNING: the distance label is not discrete: ', response
			continue

		cutoff = config.GetCutoffs(response)

		## determine the last distance bin
		rc = min(cutoff[-1], largestDistance) - 0.001
		if (rc<10.0):
			print 'ERROR: the largest distance cutoff for DFIRE is too small: ', rc
			exit(1)
		rc_index = DistanceUtils.LabelsOfOneDistance(rc, cutoff)

		binwidths = [ d2 - d1 for d1, d2 in zip(cutoff[:-1], cutoff[1:]) ]
		bincenters = [ (d2 + d1)/2. for d1, d2 in zip(cutoff[:-1], cutoff[1:]) ]

		## calculate reference potential defined as alpha*log (r/rc) + log(\delta r/ \delta rc)
		## \delta(r) is binwidths and r is the bincenters
		refPot = alpha * np.log( bincenters / bincenters[rc_index]) + np.log( binwidths / binwidths[rc_index] )

		## idx is the index for a bin
		def CalcApproxRefPot(idx=0):
                        points = np.arange(cutoff[idx] + 0.5/2, cutoff[idx+1], 0.5)
                        values = np.power(points / bincenters[rc_index], alpha)
			avg = np.average(values)
                        tmpRefPot = np.log(avg) + np.log( binwidths[idx] / binwidths[rc_index] )
			return tmpRefPot

		## get a more accurate estimation of reference for the bin with a large width
		for i in range(len(binwidths)):
			if binwidths[i] >= 1:
				refPot[i] = CalcApproxRefPot(i)
		
		## calculate the observed potential defined as log( p(r) /p(rc) ) where p(r) is the predicted distance probability
		predProb = predDistMatrix[response]
		predProbRC = predProb[:, :, rc_index : rc_index+1]
		#obsPot = np.log(predProb / (sys.float_info.min + predProbRC))
		obsPot = np.log(predProb / predProbRC)

		## calculate the final potential, which is the difference between reference potential and observed potential
		potential = np.zeros_like(predDistMatrix[response])
		potential[:, :, :rc_index ] = refPot[: rc_index] - obsPot[:, :, :rc_index]

		if subType.endswith('Plus'):
			validProb = 1 - predProb[:, :, -1]
		else:
			validProb = np.ones((predProb.shape[0], predProb.shape[1]), dtype=np.float32)

		##if useWeight=True and the prob of being disorder exists, adjust potential by the prob of not being in disorder status
		if useWeight and subType.endswith('Plus'):
			potential *= validProb[:, :, np.newaxis]

		## remove the potential for the last distance bin, which corresponds to disorder status
		if subType.endswith('Plus'):
			potential = potential[:, :, :-1]

		CheckPotentialValues(m=potential)

		potentials[response] = potential.astype(np.float32)
		validProbs[response] = validProb.astype(np.float32)

	return potentials, validProbs

def CalcPotentialByDOPE(predDistMatrix, largestDistance=20, rgScale=1., useWeight=False, minPotential=-30., maxPotential=30.):
	potentials = dict()
	validProbs = dict()
	for response in predDistMatrix.keys():
		labelName, labelType, subType = config.ParseResponse(response)
                if labelName not in config.allAtomPairNames:
                        #print 'WARNING: unsupported response for DOPE potential: ', response
                        continue
		if not conifg.IsDiscreteLabel(labelType):
			continue

                cutoff = config.GetCutoffs(response)

		## determine the last distance bin
                rc = min(cutoff[-1], largestDistance) - 0.001
                if (rc<10.0):
                        print 'ERROR: the largest distance cutoff for DOPE is too small: ', rc
                        exit(1)
                rc_index = DistanceUtils.LabelsOfOneDistance(rc, cutoff)

		binwidths = [ d2 - d1 for d1, d2 in zip(cutoff[:-1], cutoff[1:]) ]
		bincenters = [ (d2 + d1)/2. for d1, d2 in zip(cutoff[:-1], cutoff[1:]) ]

		## a is the radius of reference sphere and rg is the estimated radius of gyration
		length = predDistMatrix[response].shape[0]
		rg = 0.395*length**(3./5)+7.257	
		a = np.sqrt(5./3) * rg * rgScale

		""" calculate n(r,a) defined in the DOPE paper. Below is the original formulation.
		## rc is the upper bound of distance between two atoms
		rc = bincenters[-1]
		if rc <= 2*a:
			#nra = 6. * np.square(bincenters * (bincenters - 2*a)) * (bincenters + 4*a) / np.power(rc,3) /(np.power(rc, 3) - 18 * np.square(a)*rc + 32 * np.power(a, 3))
		else:
			#nra = 3* np.square(bincenters * (bincenters - 2*a)) * (bincenters + 4*a) / 16. / np.power(a, 6)
		"""
		## calculate n(r,a) described in the DOPE paper. Ignore the constant factor and the denominator since they are same for all distance bins
		nra = np.square(bincenters * (bincenters - 2*a)) * (bincenters + 4*a) 

		def CalcApproxRefPot(idx=0):
			points = np.arange(cutoff[idx] + 0.5/2, cutoff[idx+1], 0.5)
			values = np.square(points * (points - 2*a)) * (points + 4*a) 
			tmpNra = np.average(values)	
			return tmpNra

		## get a more accurate estimation of nra for the first several bins if their binwidth is > 0.5
		for i in range(len(binwidths)):
			if binwidths[i] >= 1:
				nra[i] = CalcApproxRefPot(i) 

		## calculate reference potential defined as log (nra(r)/nra(rc)) + log(\delta r/ \delta rc)
		## \delta(r) is equal to binwidths
		refPot = np.log( nra / nra[rc_index] * binwidths / binwidths[rc_index] )
		
	  	## calculate the observed potential defined as log( p(r) /p(rc) ) where p(r) is the predicted distance probability
                predProb = predDistMatrix[response]
                predProbRC = predProb[:, :, rc_index : rc_index+1]
                obsPot = np.log(predProb / predProbRC)

                ## calculate the final potential, which is the difference between reference and observed potential
                potential = np.zeros_like(predDistMatrix[response])
                potential[:, :, :rc_index ] = refPot[: rc_index] - obsPot[:, :, :rc_index]

		if subType.endswith('Plus'):
			validProb = 1 - predProb[:, :, -1]
		else:
			validProb = np.ones((predProb.shape[0], predProb.shape[1]), dtype=np.float32)

		##if useWeight and the prob of disroder exists, adjust potential by prob of not beining in disorder status
		if useWeight and subType.endswith('Plus'):
			potential *= validProb[:, :, np.newaxis]

		## remove the potential for the last distance bin, which corresponds to disorder status
		if subType.endswith('Plus'):
			potential = potential[:, :, :-1]

		CheckPotentialValues(m=potential)

		potentials[response] = potential.astype(np.float32)
		validProbs[response] = validProb.astype(np.float32)

	return potentials, validProbs
		
## this function is not well developed. Some revision is needed.
def CalcPotentialBySimuRW(predDistMatrix, refFile, largestDistance=20, sequence=None, useWeight=False, minPotential=-30., maxPotential=30.):
	f=open(refFile, 'rb')
	refData = cPickle.load(f)
	f.close()

	potentials = dict()
        for response in predDistMatrix.keys():
		labelName, labelType, _ = config.ParseResponse(response)
		if labelName not in config.allAtomPairNames:
                        #print 'WARNING: unsupported response for SimuRW potential: ', response
			continue
		if not conifg.IsDiscreteLabel(labelType):
			continue

                predProb = predDistMatrix[response]

		## the first row of refProb corresponds to offset=1
                refProb = refData[response]
		if labelName != 'CbCb':
			print 'distance label name not supported yet: ', labelName
			exit(1)

		if not subType.endswith('34C'):
			print 'distance label type not supported yet: ', subType
			exit(1)

		cutoff = config.GetCutoffs(response)

		length = predProb.shape[0]
		numLabels = predProb.shape[2]
		assert numLabels == refProb.shape[1]

		## maxAllowedDist[offset] is the maximum physically feasible distance between two Cb atoms when their sequence separation is equal to offset
		maxAllowedDist = [ (offset * 3.8 + 3.06) for offset in range(length) ]
		maxAllowedDist[0] = 0
		eps = 0.00001
		maxAllowedDist[2] = 10.5 - eps
		maxAllowedDist[3] = 13.0 - eps
		maxAllowedDist[4] = 15.5 - eps
		maxAllowedDist[5] = 17.5 - eps
		maxAllowedDist[6] = 19.5 - eps

		potential = np.zeros_like(predProb)

		for i in range(0, length):
			for j in range(i+2, length):
				offset = j-i
				## find the distance bin into which the maxAllowedDist falls
				lastDistBin = DistanceUtils.LabelsOfOneDistance(maxAllowedDist[offset], cutoff)
				if lastDistBin < (numLabels - 1):
					## merge the pred prob and ref prob in the bins from lastDistBin to the end
					pred = predProb[i, j,  : lastDistBin+1]
					ref = refProb[offset-1][:lastDistBin+1]

					potential[i, j, :lastDistBin+1] = -np.log( pred / ref )
					potential[i, j, lastDistBin+1: ] = maxPotential
				else:
					## determine the last distance bin
                			rc = min(cutoff[-1], largestDistance) - 0.001
                			if (rc<10.0):
                        			print 'ERROR: the largest distance cutoff for SimuRW is too small: ', rc
                        			exit(1)
                			rc_index = DistanceUtils.LabelsOfOneDistance(rc, cutoff)

					refProbLen = refProb.shape[0]
					#idx4rc = numLabels - 2
					potential[i, j] = -np.log( predProb[i, j] / refProb[min(offset, refProbLen) -1 ] )
					potential[i, j] -= potential[i, j, rc_index]
					potential[i, j, rc_index + 1: ] = 0

				## only valid for symmetric atom pairs
				potential[j, i] = potential[i, j]

		if useWeigt and subType.endswith('Plus'):
			potential *= (1-predProb[:, :, -1])

                CheckPotentialValues(potential)

		potentials[response] = potential

        return potentials

## This function calculates DFIRE/DOPE potential together with Orientation potential
def CalcDistOriPotential(predData, labelNames=['CaCa', 'CbCb', 'NO'] + ['Ca1Cb1Cb2Ca2','N1Ca1Cb1Cb2','Ca1Cb1Cb2'], distPotType='DFIRE', param4Potential=1.61, largestDistance=18, useWeight4Dist=True, useRef4Ori=True, useWeight4Ori=True, minPotential=-30, maxPotential=30):
	assert distPotType.upper() in ['DFIRE', 'DOPE']

	predProbMatrix, labelWeight, labelDistribution = predData

	validDistribution = dict()
   	validLabelWeight = dict()
        validLabelDistribution = dict()

        existingLabelNames = []
        for response, pred in predProbMatrix.iteritems():
                labelName,_, _ = config.ParseResponse(response)
                if labelName not in labelNames:
                        continue
                existingLabelNames.append(labelName)
                validDistribution[response] = pred
                validLabelWeight[response] = labelWeight[response]
                validLabelDistribution[response] = labelDistribution[response]

        missingLabelNames = list(set(labelNames) - set(existingLabelNames))
        if len(missingLabelNames)>0:
                print 'WARNING: the predicted probability file does not have information for the following label names: ', missingLabelNames

        pairPotential = dict()
        validProb = dict()

	if distPotType == 'DOPE':
		distPotential, distValidProb = CalcPotentialByDOPE(validDistribution, largestDistance=rc, rgScale=param4Potential, useWeight=useWeight4Dist, minPotential=minPotential, maxPotential=maxPotential)
	else:
        	distPotential, distValidProb = CalcPotentialByDFIRE(validDistribution, alpha=param4Potential, largestDistance=largestDistance, useWeight=useWeight4Dist, minPotential=minPotential, maxPotential=maxPotential)
        pairPotential.update(distPotential)
	validProb.update(distValidProb)

        oriPotential, oriValidProb = CalcOrientationPotential(validDistribution, useRef=useRef4Ori, useWeight=useWeight4Ori, labelWeight=validLabelWeight, labelDistribution=validLabelDistribution, minPotential=minPotential, maxPotential=maxPotential)
        pairPotential.update(oriPotential)
        validProb.update(oriValidProb)

	cutoffs = dict()
	for response in pairPotential.keys():
		cutoffs[response] = config.GetCutoffs(response)

	return pairPotential, cutoffs, validProb, distPotential, oriPotential


allRefTypesWithFiles = [ ref.upper() for ref in ['SimuRW'] ]
allRefTypes = [ 'DFIRE', 'DOPE' ] 

def main(argv):

    	inputFile = None
    	targetName = None
	labelNames = config.allAtomPairNames + config.allOrientationNames
	potentialFileSuffix = 'pkl'
	minPotential = -30.0
	maxPotential = 30.0

	UseWeight4Orientation = True
	UseWeight4Distance = True

	## the largest dist cutoff
	rc = 18

	alpha4DFIRE = 1.61
	alpha4DFIREstr = '1.61'

	rgScale4DOPE = 1.

	## reference 
	reference = 'DFIRE'

	##
	UseRef4Orientation = True

	## refFile for SimuRW
	refFile = None

	#savefolder = os.getcwd()
	savefile=""

	if len(argv) < 1:
		Usage()
		exit(1)

    	try:
        	opts, args = getopt.getopt(argv,"a:w:r:l:u:f:s:o",["labelNames=", "useWeight=", "refState=", "minPotential=", "maxPotential=", "refFile=", "savefile=", "noRef4Orientation="])
        	#print opts, args
    	except getopt.GetoptError:
        	Usage()
        	exit(1)

    	if len(args) != 1:
        	Usage()
        	exit(1)

	inputFile = args[0]

    	for opt, arg in opts:
		if opt in ("-a", "--labelNames"):
			labelNames = config.ParseLabelNames(arg)

		elif opt in ("-w", "--useWeight"):
			scheme = np.int32(arg)
			UseWeight4Orientation = (2 & scheme)>0
			UseWeight4Distance = (1 & scheme)>0

		elif opt in ("-r", "--refState"):
			fields = arg.split('+')
			reference = fields[0].upper()
			if reference not in allRefTypes:
				print 'ERROR: allowed reference types: ', allRefTypes
				exit(1)

			if len(fields) > 1:
				if fields[1].isdigit():
					rc = np.int32(fields[1])
				else:
					rc = np.float32(fields[1])

				if reference  == 'DFIRE':
					if len(fields) > 2:
						alpha4DFIREstr = fields[2]
						alpha4DFIRE = np.float32(fields[2])

				elif reference == 'DOPE':
					if len(fields) > 2:
						rgScale4DOPE = np.float32(fields[2])
				elif reference == 'SimuRW'.upper():
					#rc = np.float32(fields[1])
					print 'Using SimuRW potential'
				else:
					print 'ERROR: unsupported reference format: ', arg
					exit(1)
				

		elif opt in ("-f", "--refFile"):
			refFile = arg
			if not os.path.isfile(refFile):
				print 'the provided file for reference state is not valid: ', refFile
				exit(1)

		elif opt in ("-o", "--noRef4Orientation"):
			UseRef4Orientation = False

		elif opt in ("-s", "--savefile"):
			savefile = arg

		elif opt in ("-l", "--minPotential"):
			minPotential = np.float32(arg)
		elif opt in ("-u", "--maxPotential"):
			maxPotential = np.float32(arg)

		else:
	    		Usage()
	    		exit(1)

    	if inputFile is None:
		print 'ERROR: Please provide an input file'
		exit(1)
    	if not os.path.isfile(inputFile):
		print 'ERROR: The input file does not exist: ', inputFile
		exit(1)

	if reference in allRefTypesWithFiles and refFile is None:
		print 'ERROR: The file for user-sepcified reference state is empty'
		exit(1)

	if reference == 'DFIRE':
		if alpha4DFIRE > 10:
			## take a random value between 1.57 and 1.63
			alpha4DFIRE=random.uniform(1.57, 1.63)

		print 'alpha for DFIRE potential is ', alpha4DFIRE
		if alpha4DFIRE<1.55 or alpha4DFIRE>1.75:
			print 'ERROR: alpha4DFIRE shall be between 1.55 and 1.75'
			exit(1)

	if reference == 'DOPE':
		print 'rgScale for DOPE potential is', rgScale4DOPE
		if rgScale4DOPE > 1.2 or rgScale4DOPE <0.8:
			print 'ERROR: rgScale4DOPE shall be between 0.8 and 1.2'
			exit(1)

	if UseWeight4Distance:
		print 'Use weight for distance potential'
	if UseWeight4Orientation:
		print 'Use weight for orientation potential'
	if not UseRef4Orientation:
		print 'Do not use reference for orientation'


    	content = DistanceUtils.LoadRawDistProbFile(inputFile)
	assert len(content) >=6

    	name, sequence, predictedProb, predictedContactProb, labelWeight, labelDistribution = content[:6]
	assert labelWeight is not None, "labelWeight shall not be empty"
	predData = (predictedProb, labelWeight, labelDistribution)

        targetName = os.path.basename(inputFile).split('.')[0]
	print 'Generating potential for ', targetName, 'with the following labels: ', labelNames

	filenames = [ targetName, 'pairPotential']

	if reference == 'DFIRE':
		pairPotential, cutoffs, validProb, distPotential, oriPotential = CalcDistOriPotential(predData, labelNames, distPotType='DFIRE', param4Potential=alpha4DFIRE, largestDistance=rc, useWeight4Dist=UseWeight4Distance, useRef4Ori=UseRef4Orientation, useWeight4Ori=UseWeight4Orientation, minPotential=minPotential, maxPotential=maxPotential)
		filenames.extend([reference, str(rc), alpha4DFIREstr])
	elif reference == 'DOPE':
		pairPotential, cutoffs, validProb, distPotential, oriPotential = CalcDistOriPotential(predData, labelNames, distPotType='DOPE', param4Potential=rgScale4DOPE, largestDistance=rc, useWeight4Dist=UseWeight4Distance, useRef4Ori=UseRef4Orientation, useWeight4Ori=UseWeight4Orientation, minPotential=minPotential, maxPotential=maxPotential)
		filenames.extend([reference, str(rc), str(rgScale4DOPE)])
	else:
		print 'ERROR: unimplemented potential type: ', reference
		exit(1)

	if bool(oriPotential) and UseRef4Orientation:
		filenames.append('Ref4O')

	wStr=None
	if (bool(distPotential) and UseWeight4Distance) and (bool(oriPotential) and UseWeight4Orientation):
		wStr = 'Wt4OD'
	elif bool(oriPotential) and UseWeight4Orientation:
		wStr = 'Wt4O'
	elif bool(distPotential) and UseWeight4Distance:
		wStr = 'Wt4D'

	if wStr is not None:
		filenames.append(wStr)

	filenames.append('pkl')
	if savefile == "":
		savefile = '.'.join(filenames)

	## save the result
        with open(savefile, 'wb') as fh:
		cPickle.dump((name, sequence, pairPotential, cutoffs, validProb), fh, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    	main(sys.argv[1:])
