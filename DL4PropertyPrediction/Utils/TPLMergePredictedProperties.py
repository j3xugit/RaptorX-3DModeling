import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np
import getopt

from DL4PropertyPrediction import PropertyUtils
from DL4PropertyPrediction.config import Response2LabelName, Response2LabelType

def Usage():
    	print 'python TPLMergePredictedProperties.py [-m method | -g ground_truth_file ] pkl_files '
	print '  This script merges predicted properties derived from multiple templates'
	print '  -m: specify merge method: amean for arithmetic mean (default) or gmean for geometric mean'
	print '  -g: specify the file for native labels (with suffix .nativeProperties.pkl). If provided, prediction accuracy will be evaluated'
    	print '  pkl_files: a list of input .pkl files with name like XXX.predictedProperty.pkl'
	print '		e.g., folder1/1f20A-template1.predictedProperty.pkl folder2/1f20A-template2.predictedProperty.pkl '
    	print '     Each input file contains a tuple of 4 items: name, primary sequence,  predicted property prob, predicted property value'
	print '     The result will be saved to a file named after seqName-template1-template2.predictedProperties.pkl'

def str_display(ls):
        if not isinstance(ls, (list, tuple, np.ndarray)):
                str_ls = '{0:.4f}'.format(ls)
                return str_ls

        str_ls = ['{0:.4f}'.format(v) for v in ls ]
        str_ls2 = ' '.join(str_ls)
        return str_ls2

def MergeOneProtein(inputFiles, method):
        if inputFiles is None or len(inputFiles) < 2:
                print 'ERROR: at least two sets of predicted properties are needed for merge'
                exit(1)

        seqName = None
        sequence = None

        propertyProbs = dict()
        propertyValues = dict()

	tempNames = []
        for inputFile in inputFiles:
                content = PropertyUtils.LoadPredictedProperties(inputFile)
                name0, sequence0, predictedPropertyProb, predictedPropertyValue = content

                ##add code here to check all the input files have the same protein name
		seqName0 = '-'.join(name0.split('-')[0:-1])
		tempName0 = name0.split('-')[-1]
		tempNames.append(tempName0)

		if seqName is None:
			seqName = seqName0
		else:
			assert seqName == seqName0

                if sequence is None:
                        sequence = sequence0
                else:
                        assert sequence == sequence0

                for response in predictedPropertyProb.keys():
                        if not propertyProbs.has_key(response):
                                propertyProbs[response] =[]
                        propertyProbs[response].append( predictedPropertyProb[response] )

                for response in predictedPropertyValue.keys():
                        if not propertyValues.has_key(response):
                                propertyValues[response] =[]
                        propertyValues[response].append( predictedPropertyValue[response] )
	
	finalPropertyValue={}
	finalPropertyProb={}

	for response in propertyProbs.keys():
		finalPropertyProb[response] = np.average( propertyProbs[response], axis=0)
		labelType = Response2LabelType(response)
		if labelType.startswith('Discrete'):
                        tmpresult = np.argmax(finalPropertyProb[response], axis=1)
                        finalPropertyValue[response] = PropertyUtils.Coding2String(tmpresult, response)
		else:
			finalPropertyValue[response] = np.average( propertyValues[response], axis=0)

	targetName = '-'.join( [ seqName0 ] + tempNames )
        content4save = (targetName, sequence, finalPropertyProb, finalPropertyValue)

        return finalPropertyValue, content4save

def main(argv):

	methodpool = set(['amean', 'gmean'])
	method = 'amean'

	nativefile = None

    	try:
        	opts, args = getopt.getopt(argv,"m:g:",["method=", "nativefile="])
        	print opts, args
    	except getopt.GetoptError:
        	Usage()
        	exit(1)

	inputFiles = args
	if len(inputFiles) < 2:
		Usage()
		exit(1)

    	for opt, arg in opts:
		if opt in ("-g", "--nativefile"):
			nativefile = arg
			if not os.path.isfile(nativefile):
                        	print 'ERROR: invalid native file', nativefile
                                exit(1)

		elif opt in ("-m", "--method"):
			method = arg.strip().lower()
			if method not in methodpool:
				print 'ERROR: merge method shall be one of in', methodpool
				exit(1)
		else:
	    		Usage()
	    		exit(1)

	finalPropertyValue, content4save = MergeOneProtein(inputFiles, method)

	targetName = content4save[0]
	sequence = content4save[1]

	savefile = targetName + '.predictedProperties.pkl'

        with open(savefile, 'wb') as fh:
        	cPickle.dump( content4save, fh, protocol = cPickle.HIGHEST_PROTOCOL)

	if nativefile is not None:
                error = PropertyUtils.EvaluateSinglePropertyPrediction(finalPropertyValue, nativefile)
		print 'per-residue prediction error for ', targetName, ':', error

if __name__ == "__main__":
    	main(sys.argv[1:])
