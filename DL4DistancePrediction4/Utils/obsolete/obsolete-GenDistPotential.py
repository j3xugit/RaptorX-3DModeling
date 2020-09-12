import numpy as np
import cPickle
import sys
import os
import getopt
import copy

from DL4DistancePrediction4 import config
from DL4DistancePrediction4 import DistanceUtils
from DL4DistancePrediction4 import ContactUtils
from DL4DistancePrediction4 import RangeNWeight

def Usage():
    	print 'python GenDistPotential.py -i predidctedDistDistribution_PKL [-a atomPairType] [-r refType | -f refFile] [-t] [-l minPotential] [-u maxPotential] [ -s minSeqSep]'
    	print '  -i: the PKL file containing the predicted distance probability, e.g. target.predictedDistMatrix.pkl or target.mergedDistMatrix.pkl where target is a protein name'
    	print ' 	This file contains a tuple of 6 or 7 items: name, primary sequence, predDistProb, predContactMatrix, labelWeight, labelDistribution and refDistProb (optional)'
	print '		labelWeight: the weight of labels used in training and labelDistribution is the background probability of labels in training data'
	print '		refDistProb: the distance distribution predicted by the DL model from expected input features'
	print ''
	print '  -a: atom pair type, e.g., All, CbCb (default), CaCa, CaCa+CbCb where All indicates 5 atom pair types'
	print ''
	print '  -r: reference state for distance potential, including DFIRE (default), DOPE, SimuRW, case insensitive'
	print '		Some parameters may be used for a reference state, e.g., '
	print '		DFIRE+18+1.61: 18 is the dist cutoff and 1.61 is the exponent'
	print '		DOPE+20+1.1 :  20 is the dist cutoff and 1.1 is the scale factor for radius of gyration'
	print '		SimuRW+20: 20 is the dist cutoff'
	print '  -f: file needed for reference state SimuRW'
	print ''
	print '  -l, -u: min and max potential for one pair of atoms, default values: -30 and 30, respectively. They are only used for potential check and do not impact the real potential values'
    	print '  -t: if specified, output the distance potential matrix in text format; otherwise in PKL file (default)'
	print '  -s: When generating text output, output the potential for a pair of atoms whose sequence separation >= this value(default 3). '

eps = np.finfo(np.float).eps

#check value of one potential matrix
def CheckPotentialValues(m=None, minPotential=-30., maxPotential=30):
	if m is None:
		return
	numLabels = m.shape[2]
	for i in np.arange(numLabels):
		potential = np.triu(m[:, :, i], k=5)
		if np.amax(potential ) > maxPotential:
                        print 'WARNING: possibly too big potential value: ', np.amax(potential), ' for label ', i
                if np.amin(potential ) < minPotential:
                        print 'WARNING: possibly too small potential value: ', np.amin(potential), ' for label ', i


## predDistMatrix is the predicted distance distribution by a deep learning model.
## labelWeight is the weight of the label used in training a deep learning model
## labelDistribution is the label distribution in the training data
## when labelWeight is not equal to 1, the prob values in predDistMatrix shall be corrected before being converted into distance potential
def FixPredictedDistProb(predDistMatrix, labelWeight, labelDistribution):
	newPredDistMatrix = dict()
	for response in predDistMatrix.keys():
    		fixedProb = DistanceUtils.FixDistProb( predDistMatrix[response], labelWeight[response], labelDistribution[response])
		newPredDistMatrix[response] = fixedProb

    	return newPredDistMatrix

def CalcPotentialByDFIRE(predDistMatrix, alpha=1.61, largestDistance=15, minPotential=-20, maxPotential=20):
	potentials = dict()
	for response in predDistMatrix.keys():
		labelName, labelType, subType = config.ParseResponse(response)
		if labelName not in config.allAtomPairNames:
			print 'WARNING: unsupported response for DFIRE potential: ', response
			continue
		if not conifg.IsDiscreteLabel(labelType):
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

		## idx is the index for binwidth
		def CalcApproxRefPot(idx=0):
                        points = np.arange(cutoff[idx] + 0.5/2, cutoff[idx+1], 0.5)
                        values = np.power(points / bincenters[rc_index], alpha)
			avg = np.average(values)
                        tmpRefPot = np.log(avg) + np.log( binwidths[idx] / binwidths[rc_index] )
			return tmpRefPot

		## get a more accurate estimation of reference for the first bin
		[ refPot[i] = CalcApproxRefPot(i) for i in range(len(binwidths)) if binwdiths[i] >= 1 ]
		
		## calculate the observed potential defined as log( p(r) /p(rc) ) where p(r) is the predicted distance probability
		predProb = predDistMatrix[response]
		predProbRC = predProb[:, :, rc_index : rc_index+1]
		obsPot = np.log(predProb / predProbRC)

		## calculate the final potential, which is the difference between reference potential and observed potential
		potential = np.zeros_like(predDistMatrix[response])
		potential[:, :, :rc_index ] = refPot[: rc_index] - obsPot[:, :, :rc_index]

		CheckPotentialValues(m=potential)

		potentials[response] = potential

	return potentials

def CalcPotentialByDOPE(predDistMatrix, largestDistance=20, rgScale=1., minPotential=-20., maxPotential=20.):
	potentials = dict()
	for response in predDistMatrix.keys():
		labelName, labelType, subType = config.ParseResponse(response)
                if labelName not in config.allAtomPairNames:
                        print 'WARNING: unsupported response for DOPE potential: ', response
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
		[ nra[i] = CalcApproxRefPot(i) for i in range(len(binwidths)) if binwidths[i] >= 1 ]

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

		CheckPotentialValues(m=potential)

		potentials[response] = potential

	return potentials
		

def CalcPotentialBySimuRW(predDistMatrix, userRef, largestDistance=20, sequence=None, minPotential=-20., maxPotential=20.):
	f=open(userRef, 'rb')
	refData = cPickle.load(f)
	f.close()

	potentials = dict()
        for response in predDistMatrix.keys():
		labelName, labelType, _ = config.ParseResponse(response)
		if labelName not in config.allAtomPairNames:
                        print 'WARNING: unsupported response for SimuRW potential: ', response
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

                CheckPotentialValues(potential)

		potentials[response] = potential

        return potentials

def CalcPotentialByEmpSD(predDistMatrix, userRef, largestDistance=20, sequence=None, minPotential=-20., maxPotential=20.):
	f=open(userRef, 'rb')
	refData = cPickle.load(f)
	f.close()

	potentials = dict()
        for response, predProb in predDistMatrix.iteritems():
		labelName, labelType, _ = config.ParseResponse(response)
		if labelName not in config.allAtomPairNames:
			continue
		if not conifg.IsDiscreteLabel(labelType):
			continue

                refProbList = refData[response][1]

		length = predProb.shape[0]
		if length < 400:
			refProbs = [ ref for sz, freq, ref in refProbList if sz<=1.3*length and sz>=length/1.3 ]
		else:
			refProbs = [ ref for sz, freq, ref in refProbList if sz>=350 ]

		print '#refProbMatrix: ', len(refProbs), ' for proteins with length= ', length

		refProb = np.average(refProbs, axis=0)
		potential = - np.log ( predProb / refProb )

		rc = largestDistance
		cutoff = config.GetCutoffs(response)
		lastDistBin = DistanceUtils.LabelsOfOneDistance(rc, cutoff)
		
		lastCol = potential[:, :, lastDistBin]
		potential = potential - lastCol
		potential{;, :, lastDistBin: ] = 0

		CheckPotentialValues(potential)
                potentials[response] = potential

        return potentials

def CalcPotentialByEmpSI(predDistMatrix, userRef, largestDistance=20, sequence=None, minPotential=-20., maxPotential=20.):
	f=open(userRef, 'rb')
	refData = cPickle.load(f)
	f.close()

	potentials = dict()
        for response, prdProb in predDistMatrix.iteritems():
		labelName, labelType, _ = config.ParseResponse(response)
		if labelName not in config.allAtomPairNames:
			continue
		if not conifg.IsDiscreteLabel(labelType):
			continue

                refProb = refData[response][0]
		potential = - np.log ( predProb / refProb )
		
		rc = largestDistance
		cutoff = config.GetCutoffs(response)
		lastDistBin = DistanceUtils.LabelsOfOneDistance(rc, cutoff)
		lastCol = potential[:, :, lastDistBin]
		potential = potential - lastCol
		potential[:, :, lastDistBin: ] =0

		CheckPotentialValues(potential)
                potentials[response] = potential

        return potentials

	

allRefTypesWithFiles = [ ref.upper() for ref in ['SimuRW', 'EmpSI', 'EmpSD'] ]
allRefTypes = [ 'DFIRE', 'DOPE' ] + allRefTypesWithFiles

def main(argv):

    	inputFile = None
    	targetName = None
	labelNames = ['CbCb']
	potentialFileSuffix = 'pkl'
	minPotential = -30.0
	maxPotential = 30.0
	minSeqSep = 3
	minSeqSepStr='3'

	## the largest dist cutoff
	rc = 18

	alpha4DFIRE = 1.61
	rgScale4DOPE = 1.

	## reference 
	reference = 'DFIRE'

	## refFile
	refFile = None

    	try:
        	opts, args = getopt.getopt(argv,"i:a:r:l:u:s:f:tn",["input=", "atomPairType=", "refState=", "minPotential=", "maxPotential=", "minSeqSep=", "refFile=", "textFormat=", "nonZero="])
        	print opts, args
    	except getopt.GetoptError:
        	Usage()
        	exit(1)


    	if len(opts) < 1:
        	Usage()
        	exit(1)

    	for opt, arg in opts:
		if opt in ("-i", "--input"):
	    		inputFile = arg

		elif opt in ("-a", "--atomPairType"):
			labelNames = config.ParseLabelNames(arg)

		elif opt in ("-r", "--refState"):
			fields = arg.split('+')
			reference = fields[0].upper()
			if reference not in allRefTypes:
				print 'allowed reference types: ', allRefTypes
				exit(1)

			if len(fields) > 1:
				if reference  == 'DFIRE':
					rc = np.float32(fields[1])
					if len(fields) > 2:
						alpha4DFIRE = np.float32(fields[2])

				elif reference == 'DOPE':
					rc = np.float32(fields[1])
					if len(fields) > 2:
						rgScale4DOPE = np.float32(fields[2])
				elif reference == 'SimuRW'.upper():
					rc = np.float32(fields[1])
				else:
					print 'WARNING: unsupported reference format: ', arg
				

		elif opt in ("-f", "--refFile"):
			refFile = arg
			if not os.path.isfile(refFile):
				print 'the provided file for reference state is not valid: ', refFile
				exit(1)

		elif opt in ("-l", "--minPotential"):
			minPotential = np.float32(arg)
		elif opt in ("-u", "--maxPotential"):
			maxPotential = np.float32(arg)

		elif opt in ("-s", "--minSeqSep"):
			minSeqSep = np.int32(arg)
			minSeqSepStr = arg
			if minSeqSep < 1:
				print 'ERROR: minSeqSep shall be at least 1'
				exit(1)

		elif opt in ("-t", "--textFormat"):
	    		potentialFileSuffix = '.txt'

		elif opt in ("-n", "--nonZero"):
			resetFlag = False	

		else:
	    		Usage()
	    		exit(1)

    	if inputFile is None:
		print 'Please provide an input file'
		exit(1)
    	if not os.path.isfile(inputFile):
		print 'The input file does not exist: ', inputFile
		exit(1)

	if reference in allRefTypesWithFiles and refFile is None:
		print 'The file for user-sepcified reference state is empty'
		exit(1)

        targetName = os.path.basename(inputFile).split('.')[0]

    	content = DistanceUtils.LoadRawDistProbFile(inputFile)
	assert len(content) >=6

    	name, sequence, predictedDistProb, predictedContactProb, labelWeight, labelDistribution = content[:6]
	assert labelWeight is not None, "labelWeight shall not be empty"

	## if needed, add code to here the predicted dist probability

	filenames = [ targetName, 'distPotential']
	if reference == 'DFIRE':
		potential = CalcPotentialByDFIRE(predictedDistProb, alpha=alpha4DFIRE, largestDistance=rc, minPotential=minPotential, maxPotential=maxPotential)
		filenames.extend([reference, str(rc), str(alpha4DFIRE), potentialFileSuffix])
	elif reference == 'DOPE':
		potential = CalcPotentialByDOPE(predictedDistProb, largestDistance=rc, rgScale=rgScale4DOPE, minPotential=minPotential, maxPotential=maxPotential)
		filenames.extend([reference, str(rc), str(rgScale4DOPE), potentialFileSuffix])
	elif reference == 'SimuRW'.upper():
		potential = CalcPotentialBySimuRW(predictedDistProb, refFile, largestDistance=rc, minPotential=minPotential, maxPotential=maxPotential)
		filenames.extend([reference, str(rc), potentialFileSuffix])
	else:
		print 'ERROR: unimplemented reference state: ', reference
		exit(1)

	potentialFileName = '.'.join(filenames)

	## save to PKL file
	if potentialFileName.endswith('.pkl'):
        	fh = open(potentialFileName, 'wb')
		potential_new = dict()
		distCutoffs = dict()
		for response, pot in potential.iteritems():
			labelName = config.Response2LabelName(response)
			if labelName not in set(labelNames):
				continue

			potential_new[response] = pot
			distCutoffs[response] = config.GetCutoffs(response)

		cPickle.dump((name, sequence, potential_new, distCutoffs), fh, protocol=cPickle.HIGHEST_PROTOCOL)
		fh.close()
		return

	## save to text file
	potentialFileName = targetName + '.distPotential.s' + minSeqSepStr + potentialFileSuffix
	fh = open(potentialFileName, 'w')
	fh.write('#TARGET\t' + targetName + '\n')
	fh.write('#SEQ\t' + sequence + '\n')
	fh.write('#DistanceBinBoundaries\t' + "Please check config.py" + '\n')

	for response, pot in potential.iteritems():
		labelName, labelType, subType = config.ParseResponse(response)
		if labelName not in set(labelNames):
			continue

		size = pot.shape
		for i in xrange(size[0]):
			rawPotStrs = []

			for j in xrange(i+ minSeqSep, size[1]):
				atom1, atom2 = config.SelectAtomPair(sequence, i, j, labelName)
				y = pot[i, j]

				rawPotStr = ' '.join(['AtomPair', atom1.upper(), str(i+1), atom2.upper(), str(j+1), subType] + [ "{:.4f}".format(e) for e in y ] )
				rawPotStrs.append(rawPotStr)

			if len(rawPotStrs) >0:
				fh.write('\n'.join(rawPotStrs) + '\n')

	fh.close()


if __name__ == "__main__":
    	main(sys.argv[1:])
