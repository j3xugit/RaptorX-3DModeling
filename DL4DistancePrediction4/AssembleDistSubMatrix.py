import os
import sys
import numpy as np
import cPickle
import getopt
from copy import deepcopy

import utilsNoT
import config
import DistanceUtils
import ContactUtils
from Common import SequenceUtils

## this script is mainly used for multi-domain protein contact/distance prediction matrix assembly.
## this script assembles a set of dist/contact submatrice of a multi-domain protein into a whole-chain dist matrix
## this script can also be used to mix template-based distance prediction (for some domains with good templates) and template-free distance prediction (for domains without good templates)

def Usage():
	print 'python AssembleDistSubMatrix.py [-c | -g nativeDistMatrixFile | -n newName ] seqFile subMatrixFile1 subMatrixFile2 ...'
	print '	seqFile is the original sequence file for the whole chain. Each subMatrixFile is a dist sub matrix for one domain. Overlapping among submatrices is allowed. '
	print ' all matrx files shall have the same format as a .predictedDistMatrix.pkl file'
	print ' That is, each file is a tuple of six items: name, primary sequence, predicted distance prob, predicted contact prob, labelWeights, reference probabilities'
	print ' Meanwhile the primary sequence in a subMatrixFile shall be a subsequence of the seqFile'
	print '  -g: if the native distance matrix is provided, contact prediction accuracy will be calculated. A native file usually has name suffix .atomDistMatrix.pkl'
        print '  -c: if specified, output the Cb-Cb contact prob matrix in .gcnn and CASP format (i.e. .CASP.rr). '
	print '  -n: if sepcified, using this new name to avoid overwriting the old results '
	print ' The resultant distance matrix will be saved into a file named after XXXX.predictedDistMatrix.pkl where XXXX is the new name or targetName.mixed'


def ReplaceSubMatrix(baseMatrix, subMatrix, startPos=-1):
	if startPos < 0:
		print 'Please specify a correct start position in ReplaceSubMatrix()'
		exit(-1)

	baseShape = baseMatrix.shape
	subShape = subMatrix.shape

	"""
	print baseShape
	print subShape
	print startPos
	"""

	assert len(baseShape) == len(subShape)
	assert len(baseShape) >= 2

	assert (subShape[0] + startPos) <= baseShape[0]
	assert (subShape[1] + startPos) <= baseShape[1]

	if len(baseShape) > 2:
		assert subShape[2:] == baseShape[2:]

	resultMatrix = deepcopy(baseMatrix)
	resultMatrix[startPos:startPos+subShape[0], startPos:startPos+subShape[1] ] = subMatrix

	return resultMatrix


def main(argv):
	printContactMatrix = False
        nativefile = None
	newName=None

        try:
                opts, args = getopt.getopt(argv,"cg:n:",["contact=", "nativefile=", "name="])
                print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

        if len(args) < 2:
                Usage()
                exit(1)

	seqFile = args[0]
	subMatrixFiles = args[1:]


        for opt, arg in opts:
                if opt in ("-c", "--contact"):
                        printContactMatrix = True
                elif opt in ("-g", "--nativefile"):
                        nativefile = arg
                        if not os.path.isfile(nativefile):
                                print 'The specified file does not exist or is not accessible:', nativefile
                                exit(1)
		elif opt in ("-n", "--name"):
			newName = arg
                else:
                        print Usage()
                        exit(1)

	sequence = SequenceUtils.LoadFASTAFile(seqFile)
	seqLen = len(sequence)
	targetName = os.path.basename(seqFile).split('.')[0]

	## baseMatrix and subMatrix are a tuple of 6 items
	subMatrices = []
	for subMatrixFile in subMatrixFiles:
		subMatrix = DistanceUtils.LoadRawDistProbFile(subMatrixFile)
		subMatrices.append(subMatrix)

	## replace the distance matrix
	newDistMatrices = {}

	for subMatrix in subMatrices:
		subSequence = subMatrix[1]
		index = sequence.find(subSequence)
		if index<0:
			print 'ERROR: the sequence in the subMatixFile is not a substring of the whole-chain seq file'
			exit(-1)

		for response, m in subMatrix[2].iteritems():
			if not newDistMatrices.has_key(response):
				newDistMatrices[response] = np.zeros( (seqLen, seqLen, m.shape[2]), dtype=m.dtype )
			newDistMatrices[response] = ReplaceSubMatrix(newDistMatrices[response], subMatrix[2][response], index)

	## convert distance matrix to contact matrix
	newContMatrices = {}
	for response in newDistMatrices.keys():
		apt = config.Response2LabelName(response)
                labelType = config.Response2LabelType(response)
		if not labelType.startswith('Discrete'):
			print 'ERROR: unsupported labelType in AssembleSubDistMatrix.py: ', labelType
			exit(-1)

                subType = labelType[ len('Discrete'): ]
                labelOf8 = DistanceUtils.LabelsOfOneDistance(config.ContactDefinition, config.distCutoffs[subType])
                newContMatrices[apt] = ContactUtils.Distance2Contact(newDistMatrices[response], labelOf8)

	if newName is not None:
		targetName = newName

	## save the new result
	content4save = (targetName, sequence, newDistMatrices, newContMatrices, None, None)


	savefile = targetName + '.predictedDistMatrix.pkl'
        fh = open(savefile, 'wb')
        cPickle.dump(content4save, fh, protocol = cPickle.HIGHEST_PROTOCOL)
        fh.close()

	if nativefile is not None:
                print 'nativeFile=', nativefile
                acc = ContactUtils.EvaluateSingleContactPrediction(newContMatrices, nativefile)
                print '******************contact prediction accuracy*********************'
                for k, v in acc.iteritems():
                	print targetName, k, str_display(v)


        if printContactMatrix:
                for apt, m in newContMatrices.iteritems():
                        if apt == 'CbCb':
                                contactFileName = targetName + '.gcnn'
                                contactCASPFileName = targetName + '.CASP.rr'

                                np.savetxt(contactFileName, m, fmt='%.6f', delimiter=' ')
                                ContactUtils.SaveContactMatrixInCASPFormat(targetName, sequence, m, contactCASPFileName)



if __name__ == "__main__":
        main(sys.argv[1:])

