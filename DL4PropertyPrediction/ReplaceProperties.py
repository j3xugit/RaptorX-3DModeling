import os
import sys
import numpy as np
import cPickle
import getopt

import config
import PropertyUtils
from copy import deepcopy

## this script is mainly used for multi-domain protein local property prediction assembly.
## this script replaces a propety prediction of a multi-domain protein by another more accurate template-based property prediction for a domain
## this script can also be used to mix template-based property prediction with template-free prediction

def Usage():
	print 'python ReplaceProperties.py [ -g NativeLabelFile | -n newName ] basePropertyFile subPropertyFiles'
	print '	basePropertyFile has the original complete property for the whole chain. Each subPropertyFile has property predicted for one domain'
	print ' all files shall have the same format as a .predictedProperties.pkl file'
	print ' That is, each file is a tuple of 4 items: name, primary sequence, predicted property prob, predicted property value'
	print ' The primary sequence in a subPropertyFile shall be a subsequence of the basePropertyFile'
	print '  -g: specify the file for native labels. If provided, prediction accuracy will be calculated'
	print '  -n: if specified, provide a new name to the resultant file'
	print ' The result will be saved into a file with name XXX.predictedProperties.pkl where XXX is the new name or targetName.mixed if newName is not provided'

## baseProperty and subProperty can be a string, a list or a 2D matrix, all with shape (L, ) or (L, n) where L is the sequence length
def ReplaceSubProperty(baseProperty, subProperty, startPos=-1):
	if startPos < 0:
		print 'ERROR: please specify a correct start position in ReplaceSubProperty()'
		exit(1)

	if isinstance(subProperty, basestring) or isinstance(subProperty, list):
		subLen = len(subProperty)
	else:
		subLen = subProperty.shape[0]

	if isinstance(subProperty, basestring):
		result = baseProperty[0:startPos] + subProperty + baseProperty[startPos+subLen:]
		return result
	
	result = deepcopy(baseProperty)
	result[startPos: startPos + subLen ] = subProperty
	return result	

def main(argv):
        nativefile = None
	newName=None

        try:
                opts, args = getopt.getopt(argv,"g:n:",["nativefile=", "name="])
                print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

        if len(args) < 2:
                Usage()
                exit(1)

	basePropertyFile = args[0]
	subPropertyFiles = args[1:]

        for opt, arg in opts:
                if opt in ("-g", "--nativefile"):
                        nativefile = arg
                        if not os.path.isfile(nativeffile):
                                print 'ERROR: invalid native file', nativefile
                                exit(1)

		elif opt in ("-n", "--name"):
			newName = arg
                else:
                        Usage()
                        exit(1)

	baseProperty = PropertyUtils.LoadPredictedProperties(basePropertyFile)
	targetName = baseProperty[0]
	sequence = baseProperty[1]

	subProperties = []
	for subPropertyFile in subPropertyFiles:
		subProperty = PropertyUtils.LoadPredictedProperties(subPropertyFile)
		subProperties.append(subProperty)

	## replace the property matrix or string
	newPropertyProbs = {}
	for response, m in baseProperty[2].iteritems():
		tmpResult = m
		for subProperty in subProperties:
			subSequence = subProperty[1]
			index = sequence.find(subSequence)
			if index<0:
				print 'ERROR: the sequence in the subPropertyFile is not a substring of the sequence in basePropertyFile'
				exit(1)

			if not subProperty[2].has_key(response):
				print 'ERROR: the subPropertyFile does not have property information for ', response
				exit(1)
			tmpResult = ReplaceSubProperty(tmpResult, subProperty[2][response], index)

		newPropertyProbs[response] = tmpResult

	newPropertyValues = {}
	for response, m in baseProperty[3].iteritems():
		tmpResult = m
		for subProperty in subProperties:
			subSequence = subProperty[1]
			index = sequence.find(subSequence)
			if index<0:
				print 'ERROR: the sequence in the subPropertyFile is not a substring of the sequence in basePropertyFile'
				exit(1)

			if not subProperty[3].has_key(response):
				print 'ERROR: the subPropertyFile does not have property information for ', response
				exit(1)
			tmpResult = ReplaceSubProperty(tmpResult, subProperty[3][response], index)

		newPropertyValues[response] = tmpResult

	if newName is not None:
		fileName = newName
	else:
		fileName = targetName + '-mixed'

	## save the new result
	content4save = (targetName, sequence, newPropertyProbs, newPropertyValues)

	savefile = fileName + '.predictedProperties.pkl'
        with open(savefile, 'wb') as fh:
        	cPickle.dump(content4save, fh, protocol = cPickle.HIGHEST_PROTOCOL)

	if nativefile is not None:
                print 'nativeFile=', nativefile
                acc = PropertyUtils.EvaluateSinglePropertyPrediction(newPropertyValues, nativefile)
                print '******************property prediction accuracy*********************'
                for k, v in acc.iteritems():
                	print targetName, k, str_display(v)

if __name__ == "__main__":
        main(sys.argv[1:])
