import os
import sys
import numpy as np
import cPickle
import getopt

import config
import DistanceUtils
import ContactUtils
import OrientationUtils

from Common import SequenceUtils

## this script extracts a dist/contact/orientation submatrix of a multi-domain protein by domain sequence

def Usage():
	print 'python ExtractSubMatrix.py [ -s savefolder | -n name ] baseMatrixFile domainSeqFile'
	print '	This script extracts a dist/ori submatrix from baseMatrixFile by domain sequence'
	print '	baseMatrixFile: the whole-chain dist/ori matrix file, ending with .predictedDistMatrix.pkl'
	print '	domainSeqFile: a domain sequence file in FASTA format. It may consist of up to 2 sequence segments'
	print '	-s: the folder for result saving, default current work directory'
	print '	-n: the base name of the resultant file, default the base name of domainSeqFile'
	print '	The resultant matrix will be saved into a file XXXX.predictedDistMatrix.pkl where XXXX is the base name'

## rawSeq is the original complete primary sequence of a protein. subSeq is the primary sequence of one domain.
## Each domain may have at most 2 sequence segments in the whole-chain primary sequence. 
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

	return (pos1, pos2), (size1, size2)

## this function shall be used when one domain has 2 sequence segments
## starts is a tuple of two start positions, each for one segment
## sizes is a tuple of two segment lengths
## baseMatrix is the original complete distance matrix of a protein
## subMatrix is the distance matrix of a domain with 2 segments
def ExtractSubMatrixBySegments(baseMatrix, starts, sizes):
	if starts is None :
		print 'ERROR: please specify a pair of correct start positions for ExtractSubMatrixBySegments()'
		exit(1)

	if sizes is None:
		print 'ERROR: please specify a pair of correct sizes for ExtractSubMatrixBySegments()'
		exit(1)
	
	baseShape = baseMatrix.shape
	assert len(baseShape) >= 2

	assert len(starts) == 2
	assert len(sizes) == 2

	start1 = starts[0]
	start2 = starts[1]
	size1 = sizes[0]	
	size2 = sizes[1]

	assert (start2 + size2 <= baseShape[0])
	assert (start2 + size2 <= baseShape[1])
	assert (start1 + size1 <= start2)

	domainSize = size1 + size2
	domainShape = (domainSize, domainSize) + baseShape[2:]

	subMatrix = np.zeros(domainShape, dtype=baseMatrix.dtype)

	subMatrix[:size1, :size1] = baseMatrix[start1:start1+size1, start1:start1+size1 ]
	subMatrix[size1:, size1:] = baseMatrix[start2:start2+size2, start2:start2+size2 ]
	subMatrix[size1:, :size1] = baseMatrix[start2:start2+size2, start1:start1+size1 ]
	subMatrix[:size1, size1:] = baseMatrix[start1:start1+size1, start2:start2+size2 ]

	return subMatrix

## this function shall be used when one domain has only one sequence segment
def ExtractSubMatrix(baseMatrix, start, size):
	if start < 0:
		print 'ERROR: please specify a correct start position for ExtractSubMatrix()'
		exit(1)

	baseShape = baseMatrix.shape
	assert len(baseShape) >= 2

	assert (start + size ) <= baseShape[0]
	assert (start + size ) <= baseShape[1]

	subMatrix = baseMatrix[start:start+size, start:start+size]
	return subMatrix


def main(argv):
	newName=None
	savefolder = os.getcwd()

        try:
                opts, args = getopt.getopt(argv,"s:n:",["savefolder=", "name="])
        except getopt.GetoptError:
                Usage()
                exit(1)

        if len(args) != 2:
                Usage()
                exit(1)

	baseMatrixFile = args[0]
	domainSeqFile = args[1]

        for opt, arg in opts:
		if opt in ("-s", "--savefolder"):
			savefolder = arg
		elif opt in ("-n", "--name"):
			newName = arg
                else:
                        Usage()
                        exit(1)

	baseMatrix = DistanceUtils.LoadRawDistProbFile(baseMatrixFile)
	sequence = baseMatrix[1]

	domainSeq = SequenceUtils.LoadFASTAFile(domainSeqFile)
	domainName = os.path.basename(domainSeqFile).split('.')[0]
	domainIsSingleSeg = True

	## try by assumming that this domain has only one seq segment
	index = sequence.find(domainSeq)
	if index>=0:
		location = (index, len(domainSeq) )
	else:
		domainIsSingleSeg = False
		## try by assuming that this domain has two seq segments
		location = FindIndexBySegments(sequence, domainSeq)
		if location is None:
			print 'ERROR: cannot map domain sequence to the whole-chain sequence!'
			print 'domain Seq= ', domainSeq
			print 'chain  Seq= ', sequence
			exit(1)

	domDistMatrix = {}
	domContMatrix = {}

	for response, m in baseMatrix[2].iteritems():
		if domainIsSingleSeg:
			tmpResult = ExtractSubMatrix(m, location[0], location[1])
		else:
			tmpResult = ExtractSubMatrixBySegments(m, location[0], location[1])
		domDistMatrix[response] = tmpResult

	for response, m in baseMatrix[3].iteritems():
		if domainIsSingleSeg:
			tmpResult = ExtractSubMatrix(m, location[0], location[1])
		else:
			tmpResult = ExtractSubMatrixBySegments(m, location[0], location[1])
		domContMatrix[response] = tmpResult

	if newName is None:
		fileName = domainName
	else:
		fileName = newName

	## save the new result
	content4save = (domainName, domainSeq, domDistMatrix, domContMatrix, baseMatrix[4], baseMatrix[5])

	savefile = os.path.join(savefolder, fileName + '.predictedDistMatrix.pkl')
        with open(savefile, 'wb') as fh:
        	cPickle.dump(content4save, fh, protocol = cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
        main(sys.argv[1:])
