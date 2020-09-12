import numpy as np
import sys
import os
import getopt
import cPickle

import config
import PropertyUtils

from utils import LoadFASTAFile

def Usage():
	print 'python GenAngleRestraints4CNS.py [-k energy_constant ] seq_file_in_FASTA predictedProperties_in_PKL\n'
	print ' -k: energy constant for backbone dihedral torsion angle restraint, default 5.0'


def PhiPsi2Restraints(PhiPsi, sequence, scale):
	size = PhiPsi.shape
	if len(sequence) != size[0] :
        	print 'ERROR: the sequence length does not match the PhiPsi size'
                exit(-1)

	lines = []
	for i in range(len(sequence) ):
		## for phi
		if i > 0:
			phi = '%.1f' % (PhiPsi[i][0] / np.pi * 180.)
			phiRange = '%.1f' % (PhiPsi[i][2] / np.pi * 180.)
			line = 'assign (resid ' + str(i) + ' and name c) (resid ' + str(i+1) + ' and name n) (resid ' + str(i+1) + ' and name ca) (resid ' + str(i+1) + ' and name c) ' + str(scale) +' ' + str(phi) + ' ' + str(phiRange) + ' 2 ' + ' ! phi for residue ' + str(i+1)
			lines.append(line)

		## for psi
		if i < len(sequence)-1 :
			psi = '%.1f' % (PhiPsi[i][1] / np.pi * 180.)
			psiRange = '%.1f' % (PhiPsi[i][3] / np.pi * 180.)
			line = 'assign (resid ' + str(i+1) + ' and name n) (resid ' + str(i+1) + ' and name ca) (resid ' + str(i+1) + ' and name c) (resid ' + str(i+2) + ' and name n) ' + str(scale) + ' ' + str(psi) + ' ' + str(psiRange) + ' 2 ' + ' ! psi for residue ' + str(i+1)
			lines.append(line)

	if len(lines) < 1:
		print 'ERROR: cannot generate any angle-based restraints!'
		exit(-1)

	return lines


def main(argv):

	if len(argv) < 2:
		Usage()
		exit(-1)

	scale = 5.0

	try:
		opts, args = getopt.getopt(sys.argv[1:],"k:")
	except getopt.GetoptError:
		Usage()
		exit(-1)

	for opt, arg in opts:
		if opt in ("-k", "--energyconstant"):
			scale = np.float32(arg)
		else:
			print Usage()
			exit(-1)

	if len(args)<2:
		print Usage()
		exit(-1)

	

	seqFile = args[0]
	predictedPropertyFile = args[1]


	if not os.path.isfile(seqFile):
		print 'the sequence file does not exist: ', seqFile
		exit(-1)

	if not os.path.isfile(predictedPropertyFile):
		print ' the predicted property file does not exist: ', predictedPropertyFile
		exit(-1)


	sequence = LoadFASTAFile(seqFile)
	content = PropertyUtils.LoadPredictedProperties(predictedPropertyFile)

        name, sequence2, c1, c2 = content
	assert  (sequence == sequence2)
	PhiPsi = c1['PhiPsi_vonMise2d4']

	restraints = PhiPsi2Restraints(PhiPsi, sequence, scale)

	for res in restraints:
		print res

if __name__ == "__main__":
        main(sys.argv[1:])

