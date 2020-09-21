import os
import sys
import numpy as np
import cPickle
import copy

import config
#import DataProcessor
#from DistanceUtils import LoadNativeDistMatrixFromFile
import Metrics

from utilsNoT import str_display

## load a contact matrix in text format. This matrix has L rows and L columns
def LoadContactMatrix(file=None):
    if file is None:
        print 'please provide a distance matrix file'
        exit(1)

    if not os.path.isfile(file):
        print 'The distance matrix file does not exist: ', file
        exit(1)

    content = np.genfromtxt(file, dtype=np.float32)

    return content

## this function loads contact prediction in CASP format and save it into a matrix
## here we assume the contact matrix is symmetric and 0 indicating that there is no predicted contact
def LoadContactMatrixInCASPFormat(filename):
	fh = open(filename, 'w')
 	content = [ line.strip() for line in list(fh) ]
    	fh.close()

    	if len(content) < 5:
        	print 'the input file contains fewer than 5 lines'
        	return None

    	##the first line must be "PFRMAT RR"
    	if content[0] != "PFRMAT RR":
        	print 'The first line of the input file is not PFRMAT RR'
        	return False

    	if content[1].startswith('TARGET') is not True:
        	print 'The second line of the input file shall start with TARGET.'
        	return None

    	targetName = content[1].split()[1]
    	sequence=""
    	probs = []

	for line in content[2:]:
    		if line.startswith('AUTHOR'):
        		author=line.split()[1]
			continue

    		if line.startswith('METHOD'):
			method=line.split()[1]
			continue

    		if line.startswith('MODEL'):
			modelNum=np.int32( line.split()[1] )
			assert modelNum==1, "currently only Model 1 is supported"

        	columns = line.split()
        	if len(columns) == 1:
            		sequence += columns[0]
        	elif len(columns) == 5:
            		indices = [ int(x) for x in columns[:2] ]
            		bounds = [ int(x) for x in columns[2:4] ]
            		prob = np.float32(columns[4])

            		if bounds[0] !=0 or bounds[1] !=8:
                		print 'wrong distance bound in the following line: '
                		print line
                		return None

            		if prob > 1 or prob <0 :
                		print 'The confidence score in the following line shall be between 0 and 1: '
                		print line
                		return None

            		if indices[0]<1 or indices[0]>len(sequence) or indices[1]<1 or indices[1]>len(sequence):
                		print 'The residue index in the following line is out of range: '
                		print line
                		return None

            		if indices[0] > indices[1]:
                		print 'The first index in a residue pair shall be smaller than the 2nd one:'
                		print line
                		return None

            		probs.append( indices + [ prob ] )

        	else:
            		print 'The following line in the input file has an incorrect format: '
            		print line
            		return None

	contactMatrix = np.zeros( (len(sequence), len(sequence)), dtype=np.float32)
	for p in probs:
		contactMatrix[ p[0], p[1] ] = p[2]
		contactMatrix[ p[1], p[0] ] = p[2]

	return contactMatrix, targetName, sequence


##this function assumes that the residue index starts from 1
def SaveContactMatrixInCASPFormat(target, sequence, contactMatrix, filename, distMatrix=None, probScaleFactor=config.ProbScaleFactor, author='RaptorX', method='deep dilated convolutional residual neural networks by Jinbo Xu (jinboxu@gmail.com)'):

	if probScaleFactor == 1:
		contactMatrix2 = contactMatrix
	else:
		contactMatrix2 = np.power(contactMatrix, probScaleFactor)

	## the minimum contact prob to output
	threshold4CASP = 0.01

        ## numAllowedPairs is the maximum number of pairs allowed by CASP
	numAllowedPairs = 40000

        ##here we only consider the upper triangle, so it is not accurate for nonsymmetric matrices, e.g., CaCg, NO and HB matrices
        m = np.triu(contactMatrix2, 6)
        flatten_index = np.argsort( -1. * m, axis=None)
        real_index = np.unravel_index( flatten_index, contactMatrix2.shape )

        maxNumPairs = len(real_index )

	## the minimum number of long-range residue pairs
	seqLen = len(sequence)
        minNumLRPairs = min(3*seqLen, maxNumPairs)

	## the maximum number of allowed mediumm and short-range residue paris
        maxNumMSRPairs = maxNumPairs
        if maxNumPairs > numAllowedPairs:
                maxNumMRRPairs = max(0, maxNumPairs - minNumLRPairs)

        numPairs, numMSRPairs = 0, 0

	lines = []
        ## the header information
        lines.append('PFRMAT RR')
        lines.append('TARGET ' + target)
        lines.append('AUTHOR ' + author)
        lines.append('METHOD ' + method)
	if distMatrix is None:
		lines.append('RMODE 1')
	else:
		lines.append('RMODE 2')

	lines.append('REMARK format is explained at https://predictioncenter.org/casp14/index.cgi?page=format#RR')

	## sequence information not allowed any more, so we place it in the REMARK record
        segmentLen = 50
        for i in range(0, len(sequence), segmentLen):
        	lines.append('REMARK SEQ ' + sequence[i:i+segmentLen])

        lines.append('MODEL 1')

        for i, j in zip(real_index[0], real_index[1]):
                if numPairs > numAllowedPairs:
                        break

		if i >= j:
			continue

                offset = abs(i - j)
                if offset < 6:
                        continue

		if numPairs>30000 and contactMatrix2[i, j]<threshold4CASP :
			continue

		"""
		if contactMatrix2[i, j]<0.001:
			continue
		"""

                numPairs += 1

                if offset < 24:
                        numMSRPairs += 1

                #line = ' '.join( [ str(v) for v in [i+1, j+1, 0, 8] ] +[ "%.3f" % (contactMatrix2[i, j]) ] ) + '\n'
                items = [ str(v) for v in [i+1, j+1] ] + [ "%.3f" % contactMatrix2[i, j] ] 
		if distMatrix is not None:
			distProbs = [ "{:.3f}".format(e) for e in distMatrix[i, j] ]
			items.extend(distProbs)

		line = ' '.join(items)
		lines.append(line)

	lines.append('END')

	with open(filename, 'w') as fh:
		fh.write('\n'.join(lines) )


## convert a dist prob matrix to a contact prob matrix
## labelOf8 is the cutoff label for contact
#def Distance2Contact(distProb, distLabelType):
def Distance2Contact(distProb, labelOf8=1):
        contactProb = np.sum( distProb[:, :, :labelOf8], axis=2)
        return contactProb

##calculate MCC of a predicted contact matrix using a given score cutoff
##here we consider three cases: long-range contacts, long + medium-range contacts, long + medium- + short-range contacts
def CalcMCCF1(pred=None, truth=None, probCutoff=0.5, contactCutoff=8.0):
    if pred is None:
        print 'ERROR: please provide a predicted contact matrix'
        exit(1)

    if truth is None:
        print 'ERROR: please provide a true distance matrix'
        exit(1)

    assert pred.shape == truth.shape

    ## in case the matrix is not square, e.g., interfacial contact matrix
    seqLen = pred.shape[0]
    seqLen2 = pred.shape[1]

    pred_binary = (pred>probCutoff)
    truth_binary = ( 0<truth) & (truth<contactCutoff )
    pred_truth = pred_binary * 2 + truth_binary
    numPredicted = np.sum(pred_binary)
    numTruths = np.sum(truth_binary)
    #print "#predicted=", numPredicted, "#natives=", numTruths

    mask_LR = np.triu_indices(seqLen, 24, m=seqLen2)
    mask_MLR = np.triu_indices(seqLen, 12, m=seqLen2)
    mask_SMLR = np.triu_indices(seqLen, 6, m=seqLen2)


    metrics = []
    for mask in [ mask_LR, mask_MLR, mask_SMLR]:

        res = pred_truth[mask]
	total = res.shape[0]
	count = np.bincount(res, minlength=4)
	assert (total == np.sum(count) )

	## pred=0, truth=0	
	TN = count[0]

	## pred=0, truth=1
	FN = count[1]

	## pred=1, truth=0
	FP = count[2]

	## pred=1, truth=1
	TP = count[3]

	#print TP, FP, TN, FN

	MCC = Metrics.MCC(TP, FP, TN, FN)
	F1, precision, recall = Metrics.F1(TP, FP, TN, FN)

	metrics.extend ([MCC, TP, FP, TN, FN, F1, precision, recall])

    return np.array(metrics)


##this program outputs an array of contact prediction accuracy, arranged in the order of long-, medium-, long+medium- and short-range.
## for each range, the accuracy is calculated on the top L*ratio prediction where L is the sequence length.

## pred and truth are 2D matrices. Each entry in pred is a confidence score assigned to the corresponding residue pair indicating how likely this pair forms a contact
## truth is the ground truth distance matrix. The larger the distance, the more unlikely it is a contact. It is fine that one entry has value -1.
## in this distance matrix, only the entries with value between 0 and contactCutoff are treated as contacts.

def TopAccuracy(pred=None, truth=None, ratio=[1, 0.5, 0.2, 0.1], contactCutoff=8.0):
    if pred is None:
        print 'please provide a predicted contact matrix'
        exit(1)

    if truth is None:
        print 'please provide a true distance matrix'
        exit(1)

    if pred.shape != truth.shape:
	print 'inconsistent matrix shape: pred.shape=', pred.shape, ' truth.shape=', truth.shape
	exit(1)

    pred_truth = np.dstack( (pred, truth) )

    M1s = np.ones_like(truth, dtype = np.int8)
    mask_ER = np.triu(M1s, 48)
    mask_LR = np.triu(M1s, 24)
    mask_MLR = np.triu(M1s, 12)
    mask_SMLR = np.triu(M1s, 6)
    mask_MR = mask_MLR - mask_LR
    mask_SR = mask_SMLR - mask_MLR

    seqLen = pred.shape[0]

    #accs = []
    precs = []
    recalls = []
    F1s = []

    for mask in [ mask_ER, mask_LR, mask_MR, mask_MLR, mask_SR]:

        res = pred_truth[mask.nonzero()]
	if res.size == 0:
		#accs.extend( [0.0] * len(ratio) )
		precs.extend( [0.0] * len(ratio) )
		recalls.extend( [0.0] * len(ratio) )
		F1s.extend( [0.0] * len(ratio) )
		continue

        res_sorted = res [ (-res[:,0]).argsort() ]

	allLabels = res_sorted[:, 1]
	numNatives = ( (0<allLabels) & (allLabels<contactCutoff) ).sum()

        for r in ratio:
            numTops = int(seqLen * r)
            numTops = min(numTops, numNatives, res_sorted.shape[0] )
            topLabels = res_sorted[:numTops, 1]
            numCorrects = ( (0<topLabels) & (topLabels<contactCutoff ) ).sum()
            #accuracy = numCorrects * 1./numTops
            #accs.append(accuracy)
	    if numTops > 0:
            	prec = numCorrects * 1./numTops
	    else:
		prec = 0.
	    precs.append(prec)

	    if numNatives > 0:
	    	recall = numCorrects * 1./numNatives
	    else:
	        recall = 0.
	    recalls.append(recall)

	    if (prec + recall)>0:
	    	F1 = 2. * prec * recall/(prec + recall)
	    else:
		F1 = 0 
	    F1s.append(F1)

    res = np.array([precs, recalls, F1s])
    return res

## Evaluate contact prediction for a single protein. 
## predictedContactMatrix is a dictionary in which each key is an atom pair, e.g., CbCb, CaCa, CgCg
## native is a python dict() which has distance martrices for several atom pair types
#def EvaluateSingleContactPrediction(predictedContactMatrix, nativefile):
def EvaluateSingleContactPrediction(predictedContactMatrix, native):

        #native = LoadNativeDistMatrixFromFile(nativefile)
        accuracy = dict()
        for labelName, pred in predictedContactMatrix.iteritems():
		if labelName not in config.allDistLabelNames:
			continue

		if labelName.startswith('HB'):
                       	#accuracy[response] = TopAccuracy(pred=predictedContactMatrix[response], truth=native[response], contactCutoff=config.MaxHBDistance)
                       	accuracy[labelName] = TopAccuracy(pred=pred, truth=native[labelName], contactCutoff=config.MaxHBDistance)
		else:
                       	#accuracy[response] = TopAccuracy(pred=predictedContactMatrix[response], truth=native[response])
                       	accuracy[labelName] = TopAccuracy(pred=pred, truth=native[labelName])

	return accuracy

## pred is a 2D contact matrix, each entry has a prob value
## native is a python dict() which has distance matrices for several atom pairs
#def EvaluateSingleCbCbContactPrediction(pred, nativefile):
def EvaluateCbCbContactPrediction(pred, native):

        #native = DistanceUtils.LoadNativeDistMatrixFromFile(nativefile)
        accuracy = TopAccuracy(pred, truth=native['CbCb'])
	return accuracy


## predictedContactMatrices is a dictionary of contact matrices. Each key is a protein name
## natives is a dict of native distance matrices. Each key is a protein name.
def EvaluateContactPredictionsByMatrix(predictedContactMatrices, allnatives):

        allaccuracy = dict()
        #for name, results in predictedContactMatrices.iteritems():
        for name, predContactMatrix in predictedContactMatrices.iteritems():
		if not allnatives.has_key(name) or allnatives[name] is None:
			continue

		print 'Evaluating contact prediction accuracy for protein ', name, ' ...'
		allaccuracy[name] = EvaluateSingleContactPrediction(predContactMatrix, allnatives[name])

        ## calculate average contact prediction accuracy
        allaccuracyByApt = dict()
        for name, results in allaccuracy.iteritems():
                for apt in results.keys():
                        if not allaccuracyByApt.has_key(apt):
                                allaccuracyByApt[apt] = [ results[apt] ]
                        else:
                                allaccuracyByApt[apt].append(results[apt])

	avgacc = dict()
	for apt, allacc in allaccuracyByApt.iteritems():
		avgacc[apt] = np.average( np.array(allacc), axis=0)

	return avgacc, allaccuracy

## evaluate contact pred accuracy by files
def EvaluateContactPredictions(predictedContactMatrices, nativefolder, fileSuffix='.native.pkl'):
	allnatives = dict()
	for name in predictedContactMatrices.keys():
		filename = os.path.join(nativefolder, name + fileSuffix)
		with open(filename, 'rb') as fh:
			nativeMatrix = cPickle.load(fh)
		if fileSuffix.endswith('native.pkl'):
			allnatives[name] = nativeMatrix['atomDistMatrix']
		else:
			allnatives[name] = nativeMatrix
	return EvaluateContactPredictionsByMatrix(predictedContactMatrices, allnatives)

## print the result returned by TopAccuracy and/or EvaluateSingleCbCbContactPrediction
def PrintSingleAPTContactAccuracy(name, acc, apt='CbCb'):
	print name, 'prec', apt, str_display(acc[0])
	print name, 'recall', apt, str_display(acc[1])
	print name, 'F1', apt, str_display(acc[2])

## print acc returned by EvaluateSignleContactPrediction
def PrintSingleContactAccuracy(name, acc):
        for k, v in acc.iteritems():
		PrintSingleAPTContactAccuracy(name, v, apt=k)

## this function prints out the average or individual contact accuracy calculated by EvaluateContactPredictions
def PrintAllContactAccuracy(avgAcc, allAccs):
	PrintSingleContactAccuracy('average', avgAcc)

        for name, acc in allAccs.iteritems():
		PrintSingleContactAccuracy(name, acc)

## calculate the average or sum of the highest top k predicted probabities where k = L*ratio
## pred is the predicted contact prob matrix with dimension L*L
## if calcAverage = True, then return average, otherwise sum
def TopContactProbSum(pred=None, ratio=[1], calcAverage=True):
	if pred is None:
        	print 'ERROR: please provide a predicted contact matrix'
        	exit(1)

    	M1s = np.ones_like(pred, dtype = np.int8)
    	mask_LR = np.triu(M1s, 24)
    	mask_MLR = np.triu(M1s, 12)

    	seqLen = pred.shape[0]

    	probSum = []
    	for mask in [ mask_LR, mask_MLR ]:
        	res = pred[mask.nonzero()]
		if res.size == 0:
			probSum.extend( [0.0] * len(ratio) )
			continue

        	res_sorted = res[ (-res).argsort() ]

        	for r in ratio:
            		numTops = np.rint(seqLen * r).astype(np.int32)
            		numTops = min(numTops, res_sorted.shape[0] )
			if calcAverage:
	    			probSum.append( sum(res_sorted[:numTops]) / numTops )
			else:
	    			probSum.append( sum(res_sorted[:numTops]) )

    	return np.array(probSum)

## derive a contact matrix from a distance matrix. Note that here in the resultant contact matrix, 1 represents contact and 0 non-contact.
## this representation is different from the discreitzation of distance matrix, in which 0 is used to indicate the first distance bin
def FromDistMatrix(distMatrix):
	cMatrix = np.ones_like(distMatrix, dtype=np.int16)
	np.putmask(cMatrix, distMatrix<=0, 0)
	np.putmask(cMatrix, distMatrix>=config.ContactDefinition, 0)
	return cMatrix

## calculate the number of contacts for each residue, excluding those contacts connecting two residues with sequence separation < minSeqSep
## note that the input cMatrix is a contact matrix, in which an entry of 1 indicates one contact and 0 non-contact
## here we use 4 as the default value of minSeqSep to be consistent with the old template file
def CalcContactNumber(contactMatrix, minSeqSep=4):
	m, n = contactMatrix.shape
	assert m==n

	cMatrix = copy.deepcopy(contactMatrix)
	np.fill_diagonal(cMatrix, 0)	
	for offset in range(1, minSeqSep):
		rng1 = np.arange(m-offset)
		rng2 = rng1 + offset
		cMatrix[rng1, rng2] = 0
		cMatrix[rng2, rng1] = 0

	return np.sum(cMatrix, axis=1)
