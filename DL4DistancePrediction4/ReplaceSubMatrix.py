import os
import sys
import numpy as np
import cPickle
import getopt

import config
import DistanceUtils
import ContactUtils
import OrientationUtils

from copy import deepcopy

## this script is mainly used for multi-domain protein contact/distance prediction matrix assembly.
## this script replaces a dist/contact/orientation submatrix of a multi-domain protein by another submatrix
## this script can also be used to mix template-based distance prediction (for some domains with good templates) and template-free distance prediction (for domains without good templates)

def Usage():
	print 'python ReplaceSubMatrix.py [ -s savefolder | -n name ] baseMatrixFile subMatrixFiles'
	print '	This script replaces a dist/ori submatrix of a multi-domain protein by another one'
	print '	It can also be used to mix template-based distance/ori prediction (for domains with good templates) and template-free distance prediction (for domains without good templates)'
	print '	baseMatrixFile: the whole-chain matrix file ending with .predictedDistMatrix.pkl'
	print '	subMatrixFiles: a list of matrix files separated by white space. Each file has a dist/ori matrix for one domain, which consists of up to 2 sequence segments'
	print '		Each matrix file is a tuple of six items: name, primary sequence, predicted distance/ori prob, predicted contact prob, labelWeight and labelDistribution'
	print '		The sequence in subMatrixFile shall be a subsequence of that in baseMatrixFile'
	print '		Overlapping between domains is allowed'
	#print '	-g: if the native distance matrix is provided, contact prediction accuracy will be calculated. A native file usually has name suffix .atomDistMatrix.pkl'
        #print '	-c: if specified, output the Cb-Cb contact prob matrix in .CM.txt and CASP format (i.e. .CASP.rr). '
	print '	-s: the folder for result saving, default current work directory. Please make sure this folder is different from where the input files to avoid overwritting'
	print '	-n: the base name for the resultnat file, default target name'
	print '	The resultant distance matrix will be saved into a file XXXX.predictedDistMatrix.pkl where XXXX is the base name'

## rawSeq is the original complete primary sequence of a protein. subSeq is the primary sequence of one domain.
## Each domain may have at most 2 sequence segments along the primary sequence. 
## Each segment shall have at least minLenSeg residues.
def FindIndexBySegments(rawSeq, subSeq, minLen4Seg=10):
	assert len(subSeq) >= (2*minLen4Seg)

	## starting positions of the first and second segments
	pos1 = -1
	pos2 = -1

	## length of the first and second segments
	size1 = 0
	size2 = 0

	## split subSeq into two segments and match them to rawSeq separately
	for mid in xrange(minLen4Seg, len(subSeq)-minLen4Seg):
		pos1 = rawSeq.find(subSeq[:mid])
		if pos1 < 0:
			continue
		size1 = mid

		pos2 = rawSeq[pos1+mid:].find(subSeq[mid:])
		if pos2 <0:
			continue
		else:
			pos2 = pos2 + pos1 + mid
			size2 = len(subSeq) - mid
			break

	if pos1 <0 or pos2<0:
		return None

	#print (pos1, pos2), (size1, size2)
	return (pos1, pos2), (size1, size2)

## this function shall be used when one domain has 2 segments
## starts is a tuple of two start positions, each for one segment
## sizes is a tuple of two segment lengths
## baseMatrix is the original complete distance matrix of a protein
## subMatrix is the distance matrix of a domain with 2 segments
def ReplaceSubMatrixBySegments(baseMatrix, subMatrix, starts=None, sizes=None):
	if starts is None :
		print 'ERROR: please specify a pair of correct start positions for ReplaceSubMatrixBySegments()'
		exit(1)

	if sizes is None:
		print 'ERROR: please specify a pair of correct sizes for ReplaceSubMatrixBySegments()'
		exit(1)
	
	baseShape = baseMatrix.shape
	subShape = subMatrix.shape

	assert len(baseShape) == len(subShape)
	assert len(baseShape) >= 2

	if len(baseShape) > 2:
		assert subShape[2:] == baseShape[2:]

	assert len(starts) == 2
	assert len(sizes) == 2

	assert (sizes[0]+sizes[1] == subShape[0] )
	assert (sizes[0]+sizes[1] == subShape[1] )
	start1 = starts[0]
	start2 = starts[1]
	size1 = sizes[0]	
	size2 = sizes[1]

	assert (start2 + size2 <= baseShape[0])
	assert (start2 + size2 <= baseShape[1])
	assert (start1 + size1 <= start2)

	resultMatrix = deepcopy(baseMatrix)

	resultMatrix[start1:start1+size1, start1:start1+size1 ] = subMatrix[:size1, :size1]
	resultMatrix[start2:start2+size2, start2:start2+size2 ] = subMatrix[size1:, size1:]
	resultMatrix[start2:start2+size2, start1:start1+size1 ] = subMatrix[size1:, :size1]
	resultMatrix[start1:start1+size1, start2:start2+size2 ] = subMatrix[:size1, size1:]

	return resultMatrix

## this function shall be used when one domain has only one sequence segment
def ReplaceSubMatrix(baseMatrix, subMatrix, startPos=-1):
	if startPos < 0:
		print 'ERROR: please specify a correct start position for ReplaceSubMatrix()'
		exit(1)

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
	#printContactMatrix = False
        #nativefile = None
	newName=None
	savefolder = os.getcwd()

        try:
                #opts, args = getopt.getopt(argv,"cg:n:",["contact=", "nativefile=", "name="])
                opts, args = getopt.getopt(argv,"s:n:",["savefolder=", "name="])
                #print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

        if len(args) < 2:
                Usage()
                exit(1)

	baseMatrixFile = args[0]
	subMatrixFiles = args[1:]

        for opt, arg in opts:
		"""
                if opt in ("-c", "--contact"):
                        printContactMatrix = True
                elif opt in ("-g", "--nativefile"):
                        nativefile = arg
                        if not os.path.isfile(nativefile):
                                print 'ERROR: invalid ground truth file for contact accuracy evaluation:', nativefile
                                exit(1)
		"""
		if opt in ("-s", "--savefolder"):
			savefolder = arg
		elif opt in ("-n", "--name"):
			newName = arg
                else:
                        Usage()
                        exit(1)

	baseMatrix = DistanceUtils.LoadRawDistProbFile(baseMatrixFile)
	sequence = baseMatrix[1]
	targetName = baseMatrix[0]

	## baseMatrix and subMatrix are a tuple of 6 items
	subMatrices = []
	for subMatrixFile in subMatrixFiles:
		subMatrix = DistanceUtils.LoadRawDistProbFile(subMatrixFile)

		## make sure that both matrix files are of the same type, although they may not equal
		if baseMatrix[4] is None:
			assert (subMatrix[4] is None)
		if baseMatrix[4] is not None:
			assert (subMatrix[4] is not None)

		subMatrices.append(subMatrix)

	## new distance and contact matrices with response as the keys
	newDistMatrices = {}
	newContMatrices = {}

	## replace the distance matrix
	for response, m in baseMatrix[2].iteritems():

		tmpResult = m
		for subMatrix, smfile in zip(subMatrices, subMatrixFiles):
			if not subMatrix[2].has_key(response):
				print 'WARNING: there is no response', response, ' in subMatrixFile: ', smfile
				#exit(1)

			subSequence = subMatrix[1]
	
			## try by assumming that this domain has only one seq segment
			index = sequence.find(subSequence)
			if index>=0:
				tmpResult = ReplaceSubMatrix(tmpResult, subMatrix[2][response], index)
				continue

			## try by assuming that this domain has two seq segments
			res = FindIndexBySegments(sequence, subSequence)
			if res is None:
				print 'ERROR: cannot map domain sequence to the whole chain sequence!'
				print '    domain Seq= ', subSequence
				print '    chain  Seq= ', sequence
				exit(1)

			tmpResult = ReplaceSubMatrixBySegments(tmpResult, subMatrix[2][response], starts=res[0], sizes=res[1])


		newDistMatrices[response] = tmpResult

		## derive contact matrix from distance matrix
		labelName, labelType, subType = config.ParseResponse(response)
		if not labelType.startswith('Discrete'):
			print 'ERROR: unsupported labelType by ReplaceSubDistMatrix.py: ', labelType
			exit(1)

		if labelName in config.allAtomPairNames:
                	labelOf8 = DistanceUtils.LabelsOfOneDistance(config.ContactDefinition, config.distCutoffs[subType])
                	newContMatrices[labelName] = ContactUtils.Distance2Contact(newDistMatrices[response], labelOf8)

		elif labelName in config.allOrientationNames:
			 newContMatrices[labelName] = OrientationUtils.DeriveOriContactMatrix(newDistMatrices[response], response)
		else:
			print 'ERROR: unsupported labelName in replaceSubDistMatrix(): ', labelName
			exit(1)

	"""
	if newName is not None:
		targetName = newName
	else:
		targetName = targetName + '.mixed'
	"""
	if newName is None:
		fileName = os.path.basename(baseMatrixFile).split('.')[0] + '-mixed'
	else:
		fileName = newName

	## save the new result
	content4save = (targetName, sequence, newDistMatrices, newContMatrices, baseMatrix[4], baseMatrix[5])


	savefile = os.path.join(savefolder, fileName + '.predictedDistMatrix.pkl')
        with open(savefile, 'wb') as fh:
        	cPickle.dump(content4save, fh, protocol = cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
        main(sys.argv[1:])
