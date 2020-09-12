import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np

import config
import DistanceUtils
import ContactUtils
from config import Response2LabelName, Response2LabelType

from utilsNoT import str_display

import getopt

def Usage():

    	print 'python TPLMergePredictedDistMatrix.py [-m method | -g ground_truth_file | -c ] pkl_files '
	print '   This script merges predicted dist matrices derived from multiple templates'
	print '  -m: specify merge method: amean for arithmetic mean (default) while gmean for geometric mean'
	print '  -g: specify the native dist matrix file (with suffix .atomDistMatrix.pkl). If provided, contact prediction accuracy will be calculated'
    	print '  -c: if specified, output the Cb-Cb contact prob matrix in txt file (with suffix .gcnn for matrix output and .CASP.rr for CASP format)'
    	print '  pkl_files: a list of input .pkl files with name like XXX.predictedDistMatrix.pkl or XXX.fixedDistMatrix,pkl'
	print '		e.g., folder1/1f20A-template1.predictedDistMatrix.pkl folder2/1f20A-template2.predictedDistMatrix.pkl '

    	print '     Each input file contains a tuple of 6 items: name, primary sequence,  predicted distance prob, predicted contact prob, labelWeights, reference probabilities'
	print '     The result file is named after targetName-template1-template2-template3.predictedDistMatrix.pkl'

"""
def str_display(ls):
        if not isinstance(ls, (list, tuple, np.ndarray)):
                str_ls = '{0:.4f}'.format(ls)
                return str_ls

        str_ls = ['{0:.4f}'.format(v) for v in ls ]
        str_ls2 = ' '.join(str_ls)
        return str_ls2
"""

def MergeOneProtein(inputFiles, method):

        if inputFiles is None or len(inputFiles) < 2:
                print 'Please provide at least two predicted matrices for merge'
                exit(-1)

        seqName = None
        sequence = None

        distProbs = dict()
        contactProbs = dict()
        labelDistributions = dict()
        labelWeights = dict()
        labelWeightFlags = []

	tempNames = []
        for inputFile in inputFiles:
                content = DistanceUtils.LoadRawDistProbFile(inputFile)

                name0, sequence0, predictedDistProb, predictedContactProb, labelWeight, labelDistribution = content

                ##add code here to check all the input files have the same protein name
		seqName0 = '-'.join(name0.split('-')[0:-1])
		tempName = name0.split('-')[-1]
		tempNames.append(tempName)

		labelWeightFlags.append( labelWeight is not None )

		if seqName is None:
			seqName = seqName0
		else:
			assert seqName == seqName0

                if sequence is None:
                        sequence = sequence0
                else:
                        assert sequence == sequence0


                for apt in predictedDistProb.keys():
                        if not distProbs.has_key(apt):
                                distProbs[apt] =[]
                        distProbs[apt].append( predictedDistProb[apt] )

                for apt in predictedContactProb.keys():
                        if not contactProbs.has_key(apt):
                                contactProbs[apt] = []
                        contactProbs[apt].append( predictedContactProb[apt] )

                if labelWeight is not None:
                        for apt in labelWeight.keys():
                                if not labelWeights.has_key(apt):
                                        labelWeights[apt] = []
                                labelWeights[apt].append( labelWeight[apt] )

                for apt in labelDistribution.keys():
                        if not labelDistributions.has_key(apt):
                                labelDistributions[apt] = []
                        labelDistributions[apt].append( labelDistribution[apt] )

        ## check consistency among labelWeightFlags
        consistent  = all( flag==labelWeightFlags[0] for flag in labelWeightFlags)
        if not consistent:
                print 'ERROR: the input matrix files have inconsistent format. Some have a labelWeight while others do not.'
                exit(-1)

        ### Ms is a dictionary, each value in Ms is a list of matrices
        ### this function calculates the geometric mean of all the matrices in the same list and the renormalize the last dim of the resultant mean

        def CalcGeometricMean( Ms ):
                result = dict()
                for apt, v in Ms.iteritems():
                        result[apt] = scipy.stats.mstats.gmean(v, axis=0)
                        tmp_sum = np.sum(result[apt], axis=-1, keepdims=True)
                        result[apt] = result[apt]/tmp_sum

                return result

        ## calculate arithmetic mean
        def CalcArithmeticMean( Ms ):
                result = dict()
                for apt, v in Ms.iteritems():
                        result[apt] = np.mean(v, axis=0)

                return result

        if method == 'amean':
                distMatrixProb = CalcArithmeticMean(distProbs)
                labelDistribution = CalcArithmeticMean(labelDistributions)
        else:
                distMatrixProb = CalcGeometricMean(distProbs)
                labelDistribution = CalcGeometricMean(labelDistributions)

	contactMatrixProb = dict()
	for k in distMatrixProb.keys():
		apt = Response2LabelName(k)
		labelType = Response2LabelType(k)

		if not labelType.startswith('Discrete'):
			print 'ERROR: this labelType currently not supported in TPLMergePredicteDistMatrix.py : ', labelType
			exit(-1)

		subType = labelType[ len('Discrete'): ]
		labelOf8 = DistanceUtils.LabelsOfOneDistance(config.ContactDefinition, config.distCutoffs[subType])
		contactMatrixProb[apt] = ContactUtils.Distance2Contact(distMatrixProb[k], labelOf8)

        if labelWeightFlags[0] is True:
                labelWeight = CalcArithmeticMean(labelWeights)

	targetName = '-'.join( [ seqName ] + tempNames )
        if labelWeightFlags[0] is True:
                content4save = (targetName, sequence, distMatrixProb, contactMatrixProb, labelWeight, labelDistribution)
        else:
                content4save = (targetName, sequence, distMatrixProb, contactMatrixProb, None, labelDistribution)

        return contactMatrixProb, content4save

def main(argv):


	method = 'amean'
	methodpool = set(['amean', 'gmean'])

    	contactFileSuffix = '.gcnn'
    	printContactMatrix = False

	nativefile = None
	printContactMatrix=False

    	try:
        	opts, args = getopt.getopt(argv,"cm:g:",["contact=", "method=", "nativefile="])
        	print opts, args
    	except getopt.GetoptError:
        	Usage()
        	exit(-1)

	inputFiles = args
	if len(inputFiles) < 2:
		Usage()
		exit(-1)


    	for opt, arg in opts:
		if opt in ("-c", "--contact"):
			printContactMatrix = True
		elif opt in ("-g", "--nativefile"):
			nativefile = arg
			if not os.path.isfile(nativefile):
                        	print 'The specified file does not exist or is not accessible:', nativefile
                                exit(-1)
		elif opt in ("-m", "--method"):
			method = arg.strip().lower()

			if method not in methodpool:
				print 'ERROR: please specify a correct method for merge. It can only be in ', methodpool
				exit(-1)

		else:
	    		print Usage()
	    		exit(-1)

	contactMatrixProb, content4save = MergeOneProtein(inputFiles, method)

	targetName = content4save[0]
	sequence = content4save[1]

	savefile = targetName + '.predictedDistMatrix.pkl'

        fh = open(savefile, 'wb')
        cPickle.dump( content4save, fh, protocol = cPickle.HIGHEST_PROTOCOL)
        fh.close()


	if nativefile is not None:
                print 'nativeFile=', nativefile
                acc = ContactUtils.EvaluateSingleContactPrediction(contactMatrixProb, nativefile)
		print '******************contact prediction accuracy*********************'
                for k, v in acc.iteritems():
                        print targetName, k, str_display(v)


    	if printContactMatrix:
		for apt, m in contactMatrixProb.iteritems():
			if apt == 'CbCb':
    				contactFileName = targetName + contactFileSuffix
				contactCASPFileName = targetName + '.CASP.rr'

    				np.savetxt(contactFileName, m, fmt='%.6f', delimiter=' ')
				ContactUtils.SaveContactMatrixInCASPFormat(targetName, sequence, m, contactCASPFileName)



if __name__ == "__main__":
    	main(sys.argv[1:])
