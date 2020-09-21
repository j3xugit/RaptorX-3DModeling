import numpy as np
import scipy
import cPickle
import sys
import os
import getopt
import copy

from DL4PropertyPrediction import config
from DL4PropertyPrediction import PropertyUtils

from Common.SequenceUtils import LoadFASTAFile

def Usage():
        print 'python GenPropertyPotential4Rosetta.py [-a propertyType | -f funcType | -w energyWeight | -s savefile | -q querySeqFile | -o ] predidctedProperty_PKL'
        print '  This script converts predicted property (e.g., phi/psi angles) to Rosetta potential'
        print '         predictedProperty_PKL: a PKL file contains a tuple of at least 3 items: name, primary sequence and predicted parameters for property distribution'
	print '		Ideally, this file shall also contain prediction of property represented as a 1-d array with element indicating the probability of disorder'
	print ''
        print '  -a: property type(default PhiPsi), others not implemented yet'
	print '  -f: functional type for the potential, currently only CHARMM, AMBERPERIODIC (default), CIRCULARHARMONIC and HARMONIC supported (see Rosetta definition)'
	print '  -w: the weight factor used for this energy (default 1.0)'
	print '	 -q: the protein sequence file. when provided, check the sequence consistency between this file and predidctedProperty_PKL'
	print '  -o: if specified, disorder information will not be used; otherwise, will be used as the weight of the potential (default)'
	print '  -s: a text file for result saving; when empty, will create a text file under current work directory with name derived from the method and some parameters'

eps = np.finfo(np.float32).eps

def GeneratePhiPsiPotential(sequence, PhiPsiList, funcType='AMBERPERIODIC', weight0=1, predDisorder=None):

	constraints = []

	for i, PhiPsi in zip(range(len(sequence)), PhiPsiList):

		if predDisorder is not None:
			weight = weight0 * (1-predDisorder[i])
		else:
			weight = weight0

		## the three elements in PhiPsi are the predicted two means (phi and psi) and two variances
		## for phi
		if i > 0:
			if funcType == 'AMBERPERIODIC':
				phi_mean = '%.4f' % (PhiPsi[0] + np.pi)
			else:
				phi_mean = '%.4f' % (PhiPsi[0])


			if funcType == 'CHARMM':
				phi_sig = '%.4f' % (weight * 2./(eps + PhiPsi[2]) )
			elif funcType == 'AMBERPERIODIC':
				phi_sig = '%.4f' % (weight /(eps + PhiPsi[2]) )
			else:
				phi_sig = '%.4f' % np.sqrt(PhiPsi[2]/(eps + weight) )


			resNum1, resNum2, resNum3, resNum4 = str(i), str(i+1), str(i+1), str(i+1)
			atomName1, atomName2, atomName3, atomName4 = 'C', 'N', 'CA', 'C'
			n_periodic = '1.0'

			if funcType in [ 'HARMONIC', 'CIRCULARHAMONIC']:
				line = ' '.join(['Dihedral', atomName1, resNum1, atomName2, resNum2, atomName3, resNum3, atomName4, resNum4, funcType, phi_mean, phi_sig])
			else:
				line = ' '.join(['Dihedral', atomName1, resNum1, atomName2, resNum2, atomName3, resNum3, atomName4, resNum4, funcType, phi_mean, n_periodic, phi_sig])
			constraints.append(line)

                ## for psi
                if i < len(sequence)-1 :
			if funcType == 'AMBERPERIODIC':
				psi_mean = '%.4f' % (PhiPsi[1] + np.pi)
			else:
				psi_mean = '%.4f' % (PhiPsi[1])

			if funcType == 'CHARMM':
				psi_sig = '%.4f' % (weight * 2./(eps + PhiPsi[3]) )
			elif funcType == 'AMBERPERIODIC':
				psi_sig = '%.4f' % (weight /(eps + PhiPsi[3]) )
			else:
				psi_sig = '%.4f' % np.sqrt(PhiPsi[3]/(eps + weight) )

			resNum1, resNum2, resNum3, resNum4 = str(i+1), str(i+1), str(i+1), str(i+2)
			atomName1, atomName2, atomName3, atomName4 = 'N', 'CA', 'C', 'N'
			n_periodic = '1.0'
			
			if funcType in [ 'HARMONIC', 'CIRCULARHAMONIC']:
				line = ' '.join(['Dihedral', atomName1, resNum1, atomName2, resNum2, atomName3, resNum3, atomName4, resNum4, funcType, psi_mean, psi_sig])
			else:
				line = ' '.join(['Dihedral', atomName1, resNum1, atomName2, resNum2, atomName3, resNum3, atomName4, resNum4, funcType, psi_mean, n_periodic, psi_sig])
			constraints.append(line)

	return constraints

def main(argv):

	propertyType = 'PhiPsi'.upper()
	funcType ='AMBERPERIODIC'

	allPropertyTypes = [ propertyType ]
	allFuncTypes = [ 'CHARMM', 'AMBERPERIODIC', 'CIRCULARHARMONIC', 'HARMONIC' ]

	inputFile = None
	targetName = None
	weight = 1.0
	wStr = 'w1'

	querySeqFile = None
	querySeq = None

	UseDisorderInfo = True
	#savefolder = os.getcwd()
	savefile = ''

	if len(argv) < 1:
		Usage()
		exit(1)
	try:
                opts, args = getopt.getopt(argv,"a:f:w:s:q:o",["propertyType=", "funcType=", "weight=", "savefile=", "querySeqFile=", "noDisorder="])
        except getopt.GetoptError:
                Usage()
                exit(1)

        if len(args) != 1:
                Usage()
                exit(1)

	inputFile = args[0]

        for opt, arg in opts:
                if opt in ("-a", "--propertyType"):
			if arg.upper() not in allPropertyTypes:
				print 'ERROR: currently only support the following property types: ', allPropertyTypes
				exit(1)
                        propertyType = arg.upper()

		elif opt in ("-f", "--funcType"):
			if arg.upper() not in allFuncTypes:
				print 'ERROR: currently only support the following func types: ', allFuncTypes
				exit(1)
			funcType = arg.upper()

		elif opt in ("-w", "--weight"):
			weight = np.float32(arg)
			wStr = 'w' + arg
			if weight < 0:
				print 'ERROR: the energy weight shall be >=0'
				exit(1)

		elif opt in ("-q", "--querySeqFile"):
			querySeqFile = arg

		elif opt in ("-s", "--savefile"):
			savefile = arg

		elif opt in ("-o", "--noDisorder"):
			UseDisorderInfo = False

                else:
                        Usage()
                        exit(1)

	assert propertyType == 'PhiPsi'.upper()

	if inputFile is None:
                print 'ERROR: please provide an input file for predicted property'
                exit(1)
        if not os.path.isfile(inputFile):
                print 'ERROR: the input file does not exist: ', inputFile
                exit(1)

	if querySeqFile is not None and os.path.isfile(querySeqFile):
		querySeq = LoadFASTAFile(querySeqFile)

        targetName = os.path.basename(inputFile).split('.')[0]

	content = PropertyUtils.LoadPredictedProperties(inputFile)
	assert len(content) >=3 
	name, sequence, predProperty = content[:3]

	if querySeq is not None and querySeq!=sequence:
		print 'ERROR: inconsistent sequences in the two files:', querySeqFile, inputFile
		exit(1)

	if not predProperty.has_key('PhiPsi_vonMise2d4'):
		print 'ERROR: the property file does not have predicted Phi/Psi: ', inputFile
                exit(1)

	PhiPsiList = predProperty['PhiPsi_vonMise2d4']
	predDisorder = None
	if UseDisorderInfo and predProperty.has_key('disorder'):
		predDisorder = predProperty['disorder']

	constraints = GeneratePhiPsiPotential(sequence, PhiPsiList, funcType=funcType, weight0=weight, predDisorder=predDisorder)

	if len(constraints)<1:
		print 'ERROR: cannot generate any constraints for Phi and Psi from ', inputFile
		exit(1)
	
	if savefile == '':
		savefile = targetName + '.PhiPsi4' + funcType + '.' + wStr + '.txt'

	with open(savefile, 'w') as f:
        	f.write('\n'.join(constraints))

if __name__ == "__main__":
        main(sys.argv[1:])

