import numpy as np
import cPickle
import sys
import os

from DL4DistancePrediction4 import config, DistanceUtils
#from config import ParseResponse, GetResponseProbDims, GetResponseValueDims

def SaveBoundInListFormat(target, sequence, boundMatrix, boundFileName):
	assert (boundMatrix.shape[0] == len(sequence))
	assert (boundMatrix.shape[1] == len(sequence))

	fh = open(boundFileName, 'w')
	fh.write('TARGET ' + target + '\n')
	fh.write('REMARK this file contains estimated Cb-Cb distance information\n')
        fh.write('AUTHOR RaptorX-Contact\n')
        fh.write('METHOD deep dilated residual networks (one variant of deep CNN). Consult jinboxu@gmail.com for details. \n')
	fh.write('REMARK first 2 columns: residue numbers, 3rd column: estimated distance, 4th columns: distance standard deviation\n')
	fh.write('REMARK the other 8 columns: estimated distance (lower and upper) deviations by 4 different methods\n')

	segmentLen = 100
        for i in range(0, len(sequence), segmentLen):
                fh.write('SEQ ' + sequence[i:i+segmentLen]+'\n')


	for i in range(len(sequence)):
		for j in range(i+3, len(sequence)):
			if boundMatrix[i,j,0] < 0:
				continue

			v = [ str(v) for v in [i+1, j+1 ] ]
			for d in boundMatrix[i, j]:
				v.append( "%.3f" % d )
			line = ' '.join(v) + '\n'
			fh.write(line)
	fh.write('END\n')
	fh.close()

def DetermineProbThresholds(fixedProb, ratio_4_15A, distCutoffs):

	size = fixedProb.shape
	M1s = np.ones( (size[0], size[1]), dtype = np.int8)
	mask_LR = np.triu(M1s, 24) + np.tril(M1s, -24)
	mask_MLR = np.triu(M1s, 12) + np.tril(M1s, -12)
	mask_SMLR = np.triu(M1s, 6) + np.tril(M1s, -6)
	mask_all = np.triu(M1s, 2) + np.tril(M1s, -2)
	mask_MR = mask_MLR - mask_LR
	mask_SR = mask_SMLR - mask_MLR
	mask_NR = mask_all - mask_SMLR

	len = fixedProb.shape[0]
	maxNums = [ np.int32(x * len * 2) for x in ratio_4_15A ]

	labelOf15 = DistanceUtils.LabelsOfOneDistance(config.InteractionLimit, distCutoffs)
	fixedProb_revised = np.sum(fixedProb[:, :, labelOf15:], axis=2)

	cutoffs = []

	for mask, maxnum in zip([ mask_LR, mask_MR, mask_SR, mask_NR], maxNums ):
		res = fixedProb_revised[mask.nonzero()]
		res_sorted = res [ res.argsort() ]

		if res_sorted.shape[0] < maxnum + 1:
			cutoffs.append( res_sorted[-1] )
		else:
			cutoffs.append( res_sorted[ maxnum ] )

	return cutoffs
	        
##estimate distance bounds. In this function we assume that the distance bin has equal width
def EstimateDistanceBound(fixedProb0, response):

	labelName, distLabelType, subType = config.ParseResponse(response)
   	if 'Plus' in distLabelType:
		## merge the last interval (which represents invalid distance) to the last second one
		size = fixedProb0.shape
		fixedProb = np.zeros( (size[0], size[1], size[2]-1), dtype=fixedProb0.dtype)
		fixedProb[:, :, : -1] = fixedProb0[:, :, : -2]
		fixedProb[:, :, -1] = np.sum(fixedProb0[:, :, -2:], axis=2 )
   	else:
		fixedProb = fixedProb0
   
    	if config.GetResponseProbDims(response) < 12:
		print 'ERROR: it is not meaningful to estimate inter-resdiue distance when the number of labels is < 12'
		exit(1)

	subType = distLabelType[len('Discrete'): ]
    	distCutoffs_original = config.distCutoffs[subType]
	distCutoffs = distCutoffs_original[1:]

    	## probability thresholds for 15 Angstrom, from long-range, to medium-range, short-range and near-range
    	## if the predicted probability for 15A is larger than the threshold, we think the distance of one residue pair shall be larger than 15A
    	thresholds_4_15A = [0.75, 0.65, 0.5, 1. ]
	"""
	if distCutoffs[-1] > 17.0:
    		thresholds_4_15A = [0.85, 0.75, 0.6, 1. ]
	"""

    	## ratio threshods for 15A. we take at most ratio*L long-, medium- and short-range distance restraints for a protein of L residues
    	ratios_4_15A = [9.5, 2.2, 2.2, 4. ]

    	## determine the real probability thresholds for 15A by ratio
    	cutoffs = DetermineProbThresholds(fixedProb, ratios_4_15A, distCutoffs_original)

    	#print 'prob cutoffs determined by ratio are: ', cutoffs

    	prob_thresholds_4_15A = [ min(x, y) for x, y in zip(thresholds_4_15A, cutoffs) ]

	labelOf15 = DistanceUtils.LabelsOfOneDistance(config.InteractionLimit, distCutoffs_original)

    	#print 'the final prob cutoffs are: ', thresholds_4_15A

	halfBinWidth = np.average( [ distCutoffs[i] - distCutoffs[i-1] for i in range(1, len(distCutoffs) ) ] ) / 2.

    	mid_distance = np.array(distCutoffs - halfBinWidth).astype(np.float32)
    	upper_boundary_distance = np.array(distCutoffs).astype(np.float32)

    	numDistBins = mid_distance.shape[0]
    	mid_distance_sq = np.square(mid_distance)

    	size = fixedProb.shape
	#print size
	#print numDistBins
    	assert size[2] == numDistBins+1

    	## estimates[:, :, 0] is the expected distance if it is less than 15A
    	## estimates[:, :, 1] is the variance
    	## estimates[:, :, 2] is the lower bound
    	## estimates[:, :, 3] is the upper bound
    	## not sure why initilize to -1
    	estimates = np.full( (size[0], size[1], 10), -1, dtype=np.float32)

    	for i in range(size[0]):
		for j in range(size[1]):
	    		offset = abs(i-j)
	    		if offset < 2:
	        		continue
	    		elif offset < 6:
                		rangeIndex = 3
	    		elif offset < 12:
				rangeIndex = 2
	    		elif offset < 24:
				rangeIndex = 1
	    		else:
				rangeIndex = 0

	    		## if the prob of this residue pair suggest that the estimated distance is likely to be >15A, then do nothing
	    		if np.sum(fixedProb[i, j, labelOf15:]) > prob_thresholds_4_15A[rangeIndex]:
				continue

	    		## renormalize the distance prob distribution by setting the prob of the largest distance bin to 0
	    		newProb = fixedProb[i, j, :numDistBins]/np.sum(fixedProb[i, j, :numDistBins])

	    		dist_mean = np.average( mid_distance, weights=newProb )
            		dist_sq_mean = np.average( mid_distance_sq, weights=newProb )
            		dist_std = np.sqrt( dist_sq_mean - np.square(dist_mean) + np.square(halfBinWidth)*1./3 )

	    		## find the bin into which dist_mean falls into
            		binIndex = 0
	    		while dist_mean > upper_boundary_distance[binIndex] :
				binIndex = binIndex + 1

            		## now dist_mean <= upper_boundary_distance[binIndex] and dist_mean > upper_boundary_distance[binIndex - 1]
            		upperProb = np.zeros( numDistBins - binIndex, dtype=np.float32 )
	    		upperProb[0] = (upper_boundary_distance[binIndex] - dist_mean)/(2*halfBinWidth) * newProb[binIndex]
	    		upperProb[1: ] = newProb[binIndex+1: ]

	    		lowerProb = np.zeros ( binIndex + 1, dtype=np.float32)
	    		lowerProb[0:binIndex] = newProb[0:binIndex]
	    		lowerProb[binIndex] = newProb[binIndex] - upperProb[0]

	    		## calculate the upper distance std
            		dist_var_upper = np.dot( np.square(mid_distance[binIndex+1:] - dist_mean) + np.square(halfBinWidth)*1./3, upperProb[1: ] ) + upperProb[0] * np.square(upper_boundary_distance[binIndex] - dist_mean)*1./3

	    		## the unnormalized distance deviation
	    		dist_std_upper = np.sqrt( dist_var_upper)

	    		## the normalized distance deviation
	    		dist_std_upper2 = np.sqrt( dist_var_upper / np.sum(upperProb) )
            		dist_std_upper3 = (np.dot( mid_distance[binIndex+1:] - dist_mean, upperProb[1: ] ) + upperProb[0] * (upper_boundary_distance[binIndex] - dist_mean)/2)/np.sum(upperProb)

	    		## the expected distance deviation
	    		dist_std_upper4 = np.sum(upperProb) * dist_std_upper2

	    		## calculate the lower distance std
	    		if binIndex == 0:
				left_boundary = distCutoffs[0] - 2*halfBinWidth
	    		else:
				left_boundary = upper_boundary_distance[binIndex-1]
            		dist_var_lower = np.dot( np.square(mid_distance[:binIndex] - dist_mean) + np.square(halfBinWidth)*1./3, lowerProb[:binIndex] ) + lowerProb[binIndex] * np.square(dist_mean - left_boundary)*1./3

	    		## the unnormalized distance deviation
	    		dist_std_lower = np.sqrt( dist_var_lower)

	    		## the normalized distance deviation
	    		dist_std_lower2 = np.sqrt( dist_var_lower / np.sum(lowerProb) )
            		dist_std_lower3 = ( np.dot( dist_mean - mid_distance[:binIndex], lowerProb[:binIndex] ) + lowerProb[binIndex] * (dist_mean - left_boundary)/2 ) / np.sum(lowerProb)

	    		## the expected distance deviation
	    		dist_std_lower4 = np.sum(lowerProb) * dist_std_lower2

			"""
			##only keep those residue pairs with estimated distance < 15 Angstrom
	    		if dist_mean >= 15.0:
				estimates[i, j] = np.array( [ -1. ] * 10 ).astype(np.float32)
	    		else:
	        		estimates[i, j] = np.array([dist_mean, dist_std, dist_std_lower, dist_std_upper, dist_std_lower2, dist_std_upper2, dist_std_lower3, dist_std_upper3, dist_std_lower4, dist_std_upper4]).astype(np.float32)
			"""
	        	estimates[i, j] = np.array([dist_mean, dist_std, dist_std_lower, dist_std_upper, dist_std_lower2, dist_std_upper2, dist_std_lower3, dist_std_upper3, dist_std_lower4, dist_std_upper4]).astype(np.float32)

    	return estimates

def EstimateDistanceBounds(predictedMatrix):
	bounds = dict()
	for response, pred in predictedMatrix.iteritems():
		labelName, labelType, subType = config.ParseResponse(response)
		if labelName not in config.allAtomPairNames:
			continue
		bounds[labelName] = EstimateDistanceBound(pred, response)
	return bounds

def Usage():
    	print 'python EstimateAtomDistBounds.py pkl_file [anything] '
    	print '	pkl_file: the input file, which is usually has suffix .predictedDistMatrix.pkl or .mergedDistMatrix.pkl '
    	print ' 	This file contains a tuple of at least 3 items: name, primary sequence, predicted dist prob matrix'
	print '		the Cb-Cb distance bound will be saved to a file XXX.bound.txt where XXX is the protein name'
    	print '	if there are more than one arguments, output the distance bound info in a file XXX.bound.pkl'

import getopt

def main(argv):
    	inputFile = None
	if len(argv) < 1:
		Usage()
		exit(1)
	inputFile = argv[0]

    	printBoundPKL = False
	if len(argv) >=2:
	    	printBoundPKL = True

    	if inputFile is None:
		print 'ERROR: Please provide an input file'
		exit(1)
    	if not os.path.isfile(inputFile):
		print 'ERROR: The input file does not exist: ', inputFile
		exit(1)

    	content = DistanceUtils.LoadRawDistProbFile(inputFile)
    	targetName, sequence, predictedDistProbMatrix, predictedContactProbMatrix = content[:4]

	""" Skip this step since in version 3, we use an unbiased deep model
	if labelWeight is not None:
    		fixedProb = dict()
		for apt in predictedDistProb.keys():
			#print 'shapes: ', predictedDistProb[apt].shape, np.array(labelWeight[apt]).shape, np.array(labelDistribution[apt]).shape
    			fixedProb[apt] = DistanceUtils.FixDistProb( predictedDistProb[apt], labelWeight[apt], labelDistribution[apt])
		
	else:
		## in this case, the probability values in predictedDistProb are already corrected
		fixedProb = predictedDistProb
	fixedProb = predictedDistProb

    	if printProbMatrix:
		probFileName = targetName + probFileSuffix
		fh = open(probFileName, 'wb')
		cPickle.dump(fixedProb, fh, protocol = cPickle.HIGHEST_PROTOCOL)
		fh.close()
	"""

	bounds = EstimateDistanceBounds(predictedDistProbMatrix)

	## output Cb-Cb bound in text format
	if bounds.has_key('CbCb'):
		boundFileName = targetName + '.bound.txt'
		boundMatrix = bounds['CbCb']
		SaveBoundInListFormat(targetName, sequence, boundMatrix, boundFileName)

    	if not printBoundPKL:
		return
        boundFileName = targetName + '.bound.pkl'
	with open(boundFileName, 'wb') as fh:
		cPickle.dump( (bounds, targetName, sequence), fh, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    	main(sys.argv[1:])
