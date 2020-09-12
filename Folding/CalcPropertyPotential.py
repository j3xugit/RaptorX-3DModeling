import numpy as np
import scipy
import cPickle
import sys
import os
import getopt
import copy
import math

from DL4PropertyPrediction import config, PropertyUtils

from DL4PropertyPrediction.Utils.CalcPhiPsiFromPDB import CalcPhiPsi

def Usage():
        print 'python CalcPropertyPotential.py -i PDB_file -e predidctedProperty_PKL [-a propertyType] [-f funcType] [-d]'
	print '  -i: the input PDB file for a 3D model '
        print '  -e: specify the propety (predicted by DL) in a PKL file with name like target.predictedProperty.pkl where target is a protein name'
        print '         This file contains a tuple of at least 3 items: name, primary sequence and predicted parameters for property distribution'
        print '  -a: property type, default PhiPsi, other properties not implemented yet'
	print '  -f: functional type for the potential, currently only CharmmPeriodic is supported (see Rosetta definition)'
	print '  -d: if specified, output score at each residue. Otherwise, only output a single score for the whole protein'

## the score function is defined as 0.5*k*(1-cos(n*(x-x0)))
## where x0 is the expected value, n is the period and k is the inverse of variance
def CharmmScore(x, x0, k, n):
	return 0.5*k*(1 - math.cos(n*(x-x0) ) )

def main(argv):
	eps = np.finfo(np.float32).eps

	propertyType = 'PhiPsi'.upper()
	funcType ='Charmm'.upper()

	allPropertyTypes = [ propertyType ]
	allFuncTypes = [ funcType ]

	inputFile = None
	predPropertyFile = None
	targetName = None
	outputDetails = False

	if len(argv) < 4:
		Usage()
		exit(1)

	try:
                opts, args = getopt.getopt(argv,"i:e:a:f:w:d",["input=", "predProperty=", "propertyType=", "funcType=", "weight=", "details="])
        except getopt.GetoptError:
                Usage()
                exit(1)

        if len(opts) < 2:
                Usage()
                exit(1)

        for opt, arg in opts:
                if opt in ("-i", "--input"):
                        inputFile = arg
		elif opt in ("-e", "--predProperty"):
			predPropertyFile = arg

                elif opt in ("-a", "--propertyType"):
			if arg.upper() not in allPropertyTypes:
				print 'ERROR: currently only support the following property types: ', allPropertyTypes
				exit(1)
                        propertyType = arg.upper()

		elif opt in ("-f", "--funcType"):
			if arg.upper() not in allFuncTypes:
				print 'ERROR: currently only support the following func types: ', allFuncTypes
				exit(1)
			funcType = arg.upper()

		elif opt in ("-d", "--details"):
			outputDetails = True

                else:
                        Usage()
                        exit(1)

	assert propertyType == 'PhiPsi'.upper()
	assert funcType == 'Charmm'.upper()

	if inputFile is None:
                print 'ERROR: Please provide a PDB file as input'
                exit(1)
        if not os.path.isfile(inputFile):
                print 'ERROR: the input PDB file does not exist: ', inputFile
                exit(1)

	if predPropertyFile is None or not os.path.isfile(predPropertyFile):
		print 'ERROR: please provide a valid file for predicted property'
		exit(1)

        targetName = os.path.basename(inputFile).split('.')[0]

	content = PropertyUtils.LoadPredictedProperties(predPropertyFile)
	assert len(content) >=3 
	name, sequence, predProperty = content[:3]

	if not predProperty.has_key('PhiPsi_vonMise2d4'):
		print 'ERROR: the property file does not contain predicted Phi and Psi information: ', inputFile
                exit(1)

	modelPhiPsi, _, _ = CalcPhiPsi(inputFile, querySeq=sequence)
	score = 0
	scores = []

	PhiPsiList = predProperty['PhiPsi_vonMise2d4']
	for i, PhiPsi, mPhiPsi in zip(range(len(sequence)), PhiPsiList, modelPhiPsi):

		phi_score = None
		psi_score = None
		## for phi
		if mPhiPsi[0] is not None:
			phi_mean = PhiPsi[0]
			phi_k  = 2./(eps + PhiPsi[2]) 
			n_periodic = 1.0
			phi_score = CharmmScore(mPhiPsi[0], phi_mean, phi_k, n_periodic)
			score += phi_score

                ## for psi
                if mPhiPsi[1] is not None:
			psi_mean = PhiPsi[1]
			psi_k = 2./(eps + PhiPsi[3])
			n_periodic = 1.0
			psi_score = CharmmScore(mPhiPsi[1], psi_mean, psi_k, n_periodic)
			score += psi_score

		scores.append((phi_score, psi_score))

	print 'totalPropertyScore = ', score
	if outputDetails:
		for res, s in zip(sequence, scores):
			print res, s

if __name__ == "__main__":
        main(sys.argv[1:])

