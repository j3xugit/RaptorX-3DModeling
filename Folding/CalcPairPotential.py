import numpy as np
import cPickle
import os
import sys
import getopt

from DL4DistancePrediction4 import config
from DL4DistancePrediction4 import DistanceUtils
from DL4DistancePrediction4 import OrientationUtils
from Common import PDBUtils

def Usage():
        print 'python CalcPairPotential.py -i PDB_file -e potentialFile_PKL -a labelNames -s minSeqSep -m maxCstDist -d'
        print '  -i: an input PDB file *.pdb '
	print '  -e: a protein-specific distance-based statistical potential in PKL format '
        print '  -a: label names (default CbCb), e.g., CbCb+CaCa+NO+TwoROri'
	print '  -s: min sequence separation (default 5) for the distance potential between two atoms'
	print '  -m: max allowed distance (default 20) such that all the distance larger than this value has energy 0'
	print '  -d: if specified, output energy for each atom pairs (default False)'


## atomDistMatrix is a 2D matrix, each entry is a real-valued distance and -1 represents an invalid distance
## potential is a 3D matrix, potential[i, j] is a vector representing the potential for two atoms of residues i and j
## apts is the set of atom pair types to be scored

def Score(pairwiseMatrix, potential, labelNames, outputDetails=False, minSeqSep=6, maxCstDist=None):
	
	totalScore = 0.0
	scores = dict()

	for response, pot in potential.iteritems():
		labelName, labelType, subType = config.ParseResponse(response)
		if labelName not in set(labelNames):
                	continue
		if not pairwiseMatrix.has_key(labelName):
			print 'WARNING: the atomDistMatrix does not have distance information for atom pair:', labelName
			continue
		if not labelType.startswith('Discrete'):
			print 'ERROR: unsupported labelType: ', labelType
			exit(1)


		pm = pairwiseMatrix[labelName]
		assert pm.shape == (pot.shape[0], pot.shape[1]), "the size of the pairwise potential not compatible with the matrix"

		if labelName in config.allAtomPairNames:
			## discretize the distance matrix, an invalid entry -1 will have the largest label number
			labelMatrix, _, _  = DistanceUtils.DiscretizeDistMatrix(pm, config.distCutoffs[subType], invalidDistanceSeparated=False )
		elif labelName in config.allOrientationNames:
			labelMatrix, _, _  = OrientationUtils.DiscretizeOrientationMatrix(pm, config.distCutoffs[subType], distMatrix=pairwiseMatrix['CbCb'], invalidEntrySeparated=False )
			

		size = pot.shape
		m = np.mgrid[0:size[0], 0:size[1]]
		scoreMatrix = pot[m[0], m[1], labelMatrix ]

		if labelName in config.allAtomPairNames and maxCstDist is not None:
			label4maxDist = DistanceUtils.LabelsOfOneDistance(maxCstDist, config.distCutoffs[subType])
			np.putmask(scoreMatrix, labelMatrix > label4maxDist, 0)

		scores[response] = np.sum(np.triu(scoreMatrix, minSeqSep))
		totalScore += scores[response]

		if outputDetails:
			## note that if the potential matrix is not symmetric, we have to do something more here
			indices = np.triu_indices(size[0], k=minSeqSep, m=size[1])
			scores = scoreMatrix[indices]
			labels = labelMatrix[indices]
			for i, j, s, label in zip(indices[0], indices[1], scores, labels):
				outinfo = [ str(i+1), str(j+1), apt, str(label), "{:.4f}".format(s) ] + [ "{:.3f}".format(v) for v in pot[i, j] ]
				outstr = ' '.join(outinfo) 
				print outstr

	return totalScore, scores


def main(argv):

	apts = ['CbCb']
	minSeqSep = 5
	maxCstDist = 20.0

	if len(argv) < 4:
		Usage()
		exit(1)

	try:
                opts, args = getopt.getopt(argv,"i:a:e:s:m:d",["input=", "labelNames=", "potential=", "minSeqSep=", "maxCstDist=", "details="])
                print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)


        if len(opts) < 2:
                Usage()
                exit(1)

	inputFile = None
	potentialFile = None
	outDetails = False

        for opt, arg in opts:
                if opt in ("-i", "--input"):
                        inputFile = arg

                elif opt in ("-a", "--labelNames"):
			apts = config.ParseLabelNames(arg)

		elif opt in ("-s", "--minSeqSep"):
			minSeqSep = np.int32(arg)
			if minSeqSep < 1:
				print 'ERROR: minSeqSep shall be at least 1'
				exit(1)

		elif opt in ("-m", "--maxCstDist"):
			maxCstDist = np.float32(arg)
			if maxCstDist < 0:
				print 'ERROR: maxCstDist shall be positive'
				exit(1)

                elif opt in ("-e", "--potential"):
                        potentialFile = arg
		elif opt in ("-d", "--details"):
			outDetails = True

                else:
                        print Usage()
                        exit(1)

	if inputFile is None:
                print 'ERROR: please provide an input file'
                exit(1)
        if not os.path.isfile(inputFile):
                print 'ERROR: the input file does not exist: ', inputFile
                exit(1)
	if not inputFile.endswith('.pdb'):
		print 'WARNING: the input file does not end with pdb. Please make sure that it is a PDB file'

	if potentialFile is None:
		print 'ERROR: please provide a protein-specific pairwise statistical potential'
		exit(1)
	if not os.path.isfile(potentialFile):
		print 'ERROR: the protein-specific pairwise statistical potential file does not exist: ', potentialFile
		exit(1)

	### load potential 
	with open(potentialFile, 'rb') as fh:
		name, sequence, potential = cPickle.load(fh)[:3]

	result, pdbseq, numMisMatches, numMatches = PDBUtils.ExtractCoordinatesNDSSPBySeq(sequence, inputFile)
        if numMisMatches > 5:
                print 'ERROR: too many mismatches between query sequence and ATOM record in ', inputFile
                exit(1)

	coordInfo, dssp = result
	coordinates, numInvalidAtoms = coordInfo
        if numInvalidAtoms.has_key('CA') and numInvalidAtoms['CA']>10:
                print 'ERROR: too many Ca atoms do not have valid 3D coordinates in ', inputFile
                exit(1)
        if numInvalidAtoms.has_key('CB') and numInvalidAtoms['CB']>10:
                print 'ERROR: too many Cb atoms do not have valid 3D coordinates in ', inputFile
                exit(1)

	atomDistMatrix = PDBUtils.CalcDistMatrix(coordinates)
	assert atomDistMatrix['seq4matrix'] == sequence

	if config.HasOrientationNames(apts):
		oriMatrix = PDBUtils.CalcTwoROriMatrix(coordinates)
        	atomOrientationMatrix.update(oriMatrix)

        	oriMatrix = PDBUtils.CalcCaOriMatrix(coordinates)
        	atomOrientationMatrix.update(oriMatrix)

		assert atomOrientationMatrix['seq4matrix'] == sequence
		atomDistMatrix.update(atomOrientationMatrix)

	s, scores = Score(atomDistMatrix, potential, labelNames=apts, outputDetails=outDetails, minSeqSep=minSeqSep, maxCstDist=maxCstDist)

	print 'totalEnergy=', s, scores, inputFile

if __name__ == "__main__":
        main(sys.argv[1:])

