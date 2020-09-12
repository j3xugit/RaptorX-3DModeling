import numpy as np
import cPickle
import os
import sys
import getopt

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import config
import DistanceUtils

def Usage():
        print 'python CalcDistPotential.py -i atomDistMatrix_PKL | PDB_file -e potentialFile_PKL -a atomPairType -s minSeqSep -m maxCstDist -d'
        print '  -i: an input PDB file *.pdb or an input atom distance matrix file *.atomDistMatrix.pkl '
	print '  -e: a protein-specific distance-based statistical potential in PKL format '
        print '  -a: atom pair type (default CbCb)'
	print '  -s: min sequence separation (default 6) between two atoms for energy calcuation'
	print '  -m: max allowed distance (default 20) such that all the distance larger than this value has energy 0'
	print '  -d: if specified, output energy for each atom pairs (default False)'


## atomDistMatrix is a 2D matrix, each entry is a real-valued distance and -1 represents an invalid distance
## potential is a 3D matrix, potential[i, j] is a vector representing the potential for two atoms of residues i and j
## apts is the set of atom pair types to be scored

def Score(atomDistMatrix, potential, labelNames, outputDetails=False, minSeqSep=6, maxCstDist=None):
	
	totalScore = 0.0
	for response, pot in potential.iteritems():
		labelName, labelType, subType = config.ParseResponse(response)
		if labelName not in set(labelNames):
                	continue
		if not atomDistMatrix.has_key(labelName):
			print 'WARNING: the atomDistMatrix does not have distance information for atom pair:', labelName
			continue
		if not labelType.startswith('Discrete'):
			print 'unsupported labelType: ', labelType
			exit(1)


		distm = atomDistMatrix[labelName]
		assert distm.shape == (pot.shape[0], pot.shape[1]), "the size of the distance-based statitical potential not compatible with the distance matrix"

		## discretize the distance matrix, an invalid entry -1 will have the largest label number
		labelMatrix, _, _  = DistanceUtils.DiscretizeDistMatrix(distm, config.distCutoffs[subType], invalidDistanceSeparated=False )

		size = pot.shape
		m = np.mgrid[0:size[0], 0:size[1]]
		scoreMatrix = pot[m[0], m[1], labelMatrix ]

		if maxCstDist is not None:
			label4maxDist = DistanceUtils.LabelsOfOneDistance(maxCstDist, config.distCutoffs[subType])
			np.putmask(scoreMatrix, labelMatrix > label4maxDist, 0)

		totalScore = np.sum(np.triu(scoreMatrix, minSeqSep))

		if outputDetails:
			## note that if the potential matrix is not symmetric, we have to do something more here
			indices = np.triu_indices(size[0], k=minSeqSep, m=size[1])
			scores = scoreMatrix[indices]
			labels = labelMatrix[indices]
			for i, j, s, label in zip(indices[0], indices[1], scores, labels):
				outinfo = [ str(i+1), str(j+1), apt, str(label), "{:.4f}".format(s) ] + [ "{:.3f}".format(v) for v in pot[i, j] ]
				outstr = ' '.join(outinfo) 
				print outstr

	return totalScore


def main(argv):

	apts = ['CbCb']
	minSeqSep = 6
	maxCstDist = 20.0

	if len(argv) < 4:
		Usage()
		exit(1)

	try:
                opts, args = getopt.getopt(argv,"i:a:e:s:m:d",["input=", "atomPairType=", "potential=", "minSeqSep=", "maxCstDist=", "details="])
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

                elif opt in ("-a", "--atomPairType"):
			apts = config.ParseAtomPairTypes(arg)
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
                print 'Please provide an atom dist matrix file'
                exit(1)
        if not os.path.isfile(inputFile):
                print 'The atom dist matrix file does not exist: ', inputFile
                exit(1)

	if potentialFile is None:
		print 'Please provide a protein-specific distance-based statistical potential'
		exit(1)
	if not os.path.isfile(potentialFile):
		print 'The protein-specific distance-based statistical potential file does not exist: ', potentialFile
		exit(1)

	### load potential 
	with open(potentialFile, 'rb') as fh:
		name, sequence, potential = cPickle.load(fh)[:3]


	if inputFile.endswith('.pdb'):
		## calculate atom dist matrix from PDB file
		from CalcAtomDistMatrixFromSeqPDB import CalcAtomDistMatrix
		atomDistMatrix = CalcAtomDistMatrix(sequence, inputFile, apts=apts)
		
	else:
		### load dist matrix 
		with open(inputFile, 'rb') as fh:
			atomDistMatrix = cPickle.load(fh)

	## make sure that the sequence information in both potential and atomDistMatrix is consistent
	assert atomDistMatrix['seq4matrix'] == sequence

	s = Score(atomDistMatrix, potential, labelNames=apts, outputDetails=outDetails, minSeqSep=minSeqSep, maxCstDist=maxCstDist)

	print 'totalDistEnergy=', s, apts, inputFile

if __name__ == "__main__":
        main(sys.argv[1:])

