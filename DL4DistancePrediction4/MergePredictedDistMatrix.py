import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np

import config
import DistanceUtils
import ContactUtils
from config import ParseResponse

import getopt

def Usage():
    	##print 'python MergePredictedDistMatrix.py [-m method | -g ground_truth_folder | -s savefile | -c ] pkl_files '
    	print 'python MergePredictedDistMatrix.py [-m method | -g ground_truth_folder | -s savefile ] pkl_files'
	print '	This script merges a list of predicted dist/orientation matrices for the same protein'
	print '	-m: algorithm for merge, amean for arithmetic mean (default) and gmean for geometric mean'
	print '	-g: the folder for ground truth. If provided, contact prediction accuracy will be calculated'
    	##print '	-c: output predicted Cb-Cb contact matrix in both matrix and CASP format, which shall end with .CASP.rr and .CM.txt'
    	print '	pkl_files: a list of input .pkl files with name like XXX.predictedDistMatrix.pkl separated by space'
	print '	    e.g., folder1/1f20A.predictedDistMatrix.pkl folder2/1f20A.predictedDistMatrix.pkl '
    	print '	    Each input file is a tuple of the same number of items: name, primary sequence,  predicted distance prob, predicted contact prob, labelWeights and labelDistribution'
	print '	    It is fine the input files may have different base names, but sequences in all the files shall be exactly same. The name in the first file will be saved to the resultant file'

def MergeOneProtein(inputFiles, method):
        if inputFiles is None or len(inputFiles) < 2:
                print 'ERROR: at least two predicted matrices are needed for merge'
                exit(1)

	targetNames = []
        for inputFile in inputFiles:
                if not os.path.isfile(inputFile):
                        print 'ERROR: invalid input file:', inputFile
                        exit(1)
                #targetNames.append(os.path.basename(inputFile).split('.')[0])
        #assert all( name==targetNames[0] for name in targetNames)

        #targetName = targetNames[0]
	proteinName = None
        sequence = None

        distProbs = dict()
        contactProbs = dict()
        labelDistributions = dict()
        labelWeights = dict()
	refDistProbs = dict()
        labelWeightFlags = []

	refDistMatrixProbExist = True

        for inputFile in inputFiles:
                content = DistanceUtils.LoadRawDistProbFile(inputFile)
                name, sequence0, predictedDistProb, predictedContactProb, labelWeight, labelDistribution = content[:6]
		if proteinName is None:
			proteinName = name

		## check if labelWeight exists or not
		labelWeightFlags.append( labelWeight is not None )

		## make sure all files have the same sequences
                if sequence is None:
                        sequence = sequence0
                else:
                        assert sequence == sequence0

                for response in predictedDistProb.keys():
                        if not distProbs.has_key(response):
                                distProbs[response] =[]
                        distProbs[response].append( predictedDistProb[response] )

		if len(content) >= 7:
                	refDistProb = content[6]
                	for response in refDistProb.keys():
                        	if not refDistProbs.has_key(response):
                                	refDistProbs[response] =[]
                        	refDistProbs[response].append( refDistProb[response] )
		else:
			refDistMatrixProbExist = False

                for labelName in predictedContactProb.keys():
                        if not contactProbs.has_key(labelName):
                                contactProbs[labelName] = []
                        contactProbs[labelName].append( predictedContactProb[labelName] )

                if labelWeight is not None:
                        for response in labelWeight.keys():
                                if not labelWeights.has_key(response):
                                        labelWeights[response] = [labelWeight[response]]
				elif labelWeight[response].shape ==  labelWeights[response][-1].shape:
                                	labelWeights[response].append( labelWeight[response] )

                for response in labelDistribution.keys():
                        if not labelDistributions.has_key(response):
                                labelDistributions[response] = [ labelDistribution[response] ]
			elif labelDistribution[response].shape ==  labelDistributions[response][-1].shape:
                        	labelDistributions[response].append( labelDistribution[response] )

        ## check consistency among labelWeightFlags
        consistent  = all( flag==labelWeightFlags[0] for flag in labelWeightFlags)
        if not consistent:
                print 'ERROR: input files have inconsistent format; some have labelWeight while others do not.'
                exit(1)

        ### Ms is a dictionary, each value in Ms is a list of matrices
        ### this function calculates the geometric mean of all the matrices in the same list and the renormalize the last dim of the resultant mean
        def CalcGeometricMean( Ms ):
                result = dict()
                for k, v in Ms.iteritems():
                        result[k] = scipy.stats.mstats.gmean(v, axis=0)
                        tmp_sum = np.sum(result[k], axis=-1, keepdims=True)
                        result[k] = result[k]/tmp_sum
                return result

        ## calculate arithmetic mean
        def CalcArithmeticMean( Ms ):
                result = dict()
                for k, v in Ms.iteritems():
                        result[k] = np.mean(v, axis=0)
                return result

        if method == 'amean':
                distMatrixProb = CalcArithmeticMean(distProbs)
                refDistMatrixProb = CalcArithmeticMean(refDistProbs)
                labelDistribution = CalcArithmeticMean(labelDistributions)
        else:
                distMatrixProb = CalcGeometricMean(distProbs)
                refDistMatrixProb = CalcGeometricMean(refDistProbs)
                labelDistribution = CalcGeometricMean(labelDistributions)

	## here handle the situation when different label types appear in the input files
	contactMatrixProb = dict()
	contactMatrixCount = dict()
	for response in distMatrixProb.keys():
		labelName, labelType, subType = ParseResponse(response)
		if labelName not in config.allAtomPairNames:
			continue

		labelOf8 = DistanceUtils.LabelsOfOneDistance(config.ContactDefinition, config.distCutoffs[subType])
		if not contactMatrixProb.has_key(labelName):
			contactMatrixProb[labelName] = ContactUtils.Distance2Contact(distMatrixProb[response], labelOf8)
			contactMatrixCount[labelName] = 1
		else:
			contactMatrixProb[labelName] += ContactUtils.Distance2Contact(distMatrixProb[response], labelOf8)
			contactMatrixCount[labelName] += 1

	for labelName in contactMatrixProb.keys():
		contactMatrixProb[labelName] = contactMatrixProb[labelName]/contactMatrixCount[labelName]

        if labelWeightFlags[0] is True:
                labelWeight = CalcArithmeticMean(labelWeights)
	else:
		labelWeight = None

	if refDistMatrixProbExist:
        	content4save = (proteinName, sequence, distMatrixProb, contactMatrixProb, labelWeight, labelDistribution, refDistMatrixProb)
	else:
        	content4save = (proteinName, sequence, distMatrixProb, contactMatrixProb, labelWeight, labelDistribution)

        return contactMatrixProb, content4save

def main(argv):

	savefile = None

	methodpool = set(['amean', 'gmean'])
	method = 'amean'

    	printContactMatrix = False
    	contactFileSuffix = '.CM.txt'

	nativefile = None

    	try:
        	opts, args = getopt.getopt(argv,"cm:g:s:",["contact=", "method=", "nativefile=", "savefile"])
        	#print opts, args
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
		elif opt in ("-g", "--nativefile"):
			nativefile = arg
			if not os.path.isfile(nativefile):
                        	print 'ERROR: the specified ground truth file does not exist or is not accessible:', nativefile
                                exit(1)
		elif opt in ("-m", "--method"):
			method = arg.strip().lower()
			if method not in methodpool:
				print 'ERROR: please specify a correct method for merge. It can only be in ', methodpool
				exit(1)
		elif opt in ("-s", "--savefile"):
			savefile = arg
		else:
	    		Usage()
	    		exit(1)

	contactMatrixProb, content4save = MergeOneProtein(inputFiles, method)

	targetName = content4save[0]
	sequence = content4save[1]

	if savefile is None:
		savefile ='.'.join([targetName, method, 'predictedDistMatrix.pkl'])
        with open(savefile, 'wb') as fh:
        	cPickle.dump( content4save, fh, protocol = cPickle.HIGHEST_PROTOCOL)

	if nativefile is not None:
                print 'nativeFile=', nativefile
                acc = ContactUtils.EvaluateSingleContactPrediction(contactMatrixProb, nativefile)
		print '******************contact prediction accuracy*********************'
		ContactUtils.PrintSingleAPTContactAccuracy(targetName, acc)

	###
    	if printContactMatrix:
		for labelName, m in contactMatrixProb.iteritems():
			if labelName == 'CbCb':
				savefolder = os.path.dirname(savefile)
    				contactFileName = os.path.join(savefolder, targetName + '.CM.txt')
				contactCASPFileName = os.path.join(savefolder, targetName + '.CASP.rr')

    				np.savetxt(contactFileName, m, fmt='%.6f', delimiter=' ')
				ContactUtils.SaveContactMatrixInCASPFormat(targetName, sequence, m, contactCASPFileName)
	###


if __name__ == "__main__":
    	main(sys.argv[1:])
