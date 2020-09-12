import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np

import config
import DistanceUtils
import ContactUtils

import getopt

def Usage():

    	print 'python MergePredictedContactMatrix.py [-m method | -g ground_truth_folder | -c ] pkl_files '
	print '  This script builds a predicted contact matrix from a set of predicted distance matrix files '
	print '  -m: specify merge method: amean for arithmetic mean (default) while gmean for geometric mean'
	print '  -g: specify the folder containing the ground truth. If provided, contact prediction accuracy will be calculated'
    	print '  -c: output the contact prob matrix in txt file. The Cb-Cb contact file has suffix .gcnn. The Ca-Ca contact file has suffix .CaCa.gcnn .'
    	print '  pkl_files: a list of input .pkl files with name like XXX.predictedDistMatrix.pkl or XXX.fixedDistMatrix.pkl'
	print '		e.g., folder1/1f20A.predictedDistMatrix.pkl folder2/1f20A.predictedDistMatrix.pkl '
    	print '     An input file contains a tuple of 6 items: name, primary sequence, predicted distance prob matrix, predicted contact prob matrix, labelWeights, reference probabilities'
    	print '     but labelWeight is None for a fixedDistMatrix.pkl file'

	print '  This script will output one or a few files. The .predictedContactMatrix.pkl file is a python dictionary containing protein name, sequence and predicted contact matrix '
	print '  If the -c option is used, the resultant .gcnn file is a contact matrix in text format.'

def MergeContactMatrix4OneProtein(inputFiles, method):

        if inputFiles is None or len(inputFiles) < 2:
                print 'Please provide at least two predicted matrices for merge'
                exit(-1)

        targetNames = []
        for inputFile in inputFiles:
                if not os.path.isfile(inputFile):
                        print 'The input file does not exist: ', inputFile
                        exit(-1)
                targetNames.append(os.path.basename(inputFile).split('.')[0])

        assert all( name==targetNames[0] for name in targetNames)

        targetName = targetNames[0]
        sequence = None

        contactProbs = dict()

        for inputFile in inputFiles:
                content = DistanceUtils.LoadRawDistProbFile(inputFile)

                name, sequence0, predictedDistProb, predictedContactProb, labelWeight, labelDistribution = content[:6]

                ##add code here to check all the input files have the same protein name
                assert ( name == targetName)

                if sequence is None:
                        sequence = sequence0
                else:
                        assert sequence == sequence0


                for apt in predictedContactProb.keys():
                        if not contactProbs.has_key(apt):
                                contactProbs[apt] = []
                        contactProbs[apt].append( predictedContactProb[apt] )

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
                contactMatrixProb = CalcArithmeticMean(contactProbs)
        else:
                contactMatrixProb = CalcGeometricMean(contactProbs)

	content4save = dict()
	content4save['name'] = targetName
	content4save['sequence'] = sequence
	content4save['predContactMatrix'] = contactMatrixProb

        return content4save


def MergeAndSaveOneProtein(inputFiles, method, printContactMatrix):

	content4save = MergeContactMatrix4OneProtein(inputFiles, method)

	targetName = content4save['name']
	sequence = content4save['sequence']
	contactMatrixProb = content4save['predContactMatrix']

	if method == 'amean':
		savefile = targetName + '.predictedContactMatrix.pkl'
	else:
		savefile = targetName + '.' + method + '.predictedContactMatrix.pkl'
        fh = open(savefile, 'wb')
        cPickle.dump( content4save, fh, protocol = cPickle.HIGHEST_PROTOCOL)
        fh.close()

    	if printContactMatrix:
    		contactFileSuffix = '.gcnn'
		for apt, m in contactMatrixProb.iteritems():
			if apt == 'CbCb':
    				contactFileName = targetName + contactFileSuffix
				contactCASPFileName = targetName + '.CASP.rr'
			
    				np.savetxt(contactFileName, m, fmt='%1.6f', delimiter=' ')
				ContactUtils.SaveContactMatrixInCASPFormat(targetName, sequence, m, contactCASPFileName)

	return content4save


def main(argv):

	methodpool = set(['amean', 'gmean'])

    	inputFiles = None
    	targetName = None
	method = 'amean'

	nativefolder = None
    	printContactMatrix = False

    	try:
        	opts, args = getopt.getopt(argv,"cm:g:",["contact=", "method=", "nativefolder="])
        	print opts, args
    	except getopt.GetoptError:
        	Usage()
        	exit(1)

	inputFiles = args
	if len(inputFiles) < 2:
		Usage()
		exit(1)


    	for opt, arg in opts:
		if opt in ("-c", "--contact"):
			printContactMatrix = True
		elif opt in ("-g", "--nativefolder"):
			nativefolder = arg
			if not os.path.isdir(nativefolder):
                        	print 'The specified folder does not exist or is not accessible:', nativefolder
                                exit(1)
		elif opt in ("-m", "--method"):
			method = arg.strip().lower()

			if method not in methodpool:
				print 'ERROR: please specify a correct method for merge. It can only be amean or gmean.'
				exit(1)

		else:
	    		print Usage()
	    		exit(1)

	result = MergeAndSaveOneProtein(inputFiles, method, printContactMatrix)

	if nativefolder is not None:
                print 'nativeFolder=', nativefolder
		contactPredictions = dict()
		contactPredictions[ result['name'] ] = result['predContactMatrix']
                acc = ContactUtils.EvaluateContactPredictions(contactPredictions, nativefolder)
		##print acc

if __name__ == "__main__":
    	main(sys.argv[1:])
