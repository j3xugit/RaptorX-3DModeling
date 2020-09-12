import numpy as np
import sys
import os
import getopt
import cPickle
import re

import config
import PropertyUtils

from utils import LoadFASTAFile

def Usage():
	print 'python GenPropertyRestraints4CNS.py [-k energy_constant | -d distRestraint_file | -a angleRestraint_file | -h hbondRestraint_file  ] seq_file_in_FASTA predictedProperties_in_PKL\n'
	print '   This script generates both angle and secondary structure restraints for a sequence'
	print '   -k: energy constant for backbone dihedral torsion angle restraints, default 5.0'
	print '   -d: file for saving distance restraints derived from secondary structure, default ssnoe.tbl'
	print '   -a: file for saving angle restraints derived from predicted Phi and Psi, default dihedral.tbl'
	print '   -h: file for saving hbonds derived from predicted helix, default hbond.tbl'

## start and end: the starting and ending positions of a helix segment, counting from 0
## in CNS, the residue index counts from 1, i.e., we shall add 1 to each position in a segment
def Helix2DistRestraints(start, end, sequence):
	lines = []
	for j in range(start+1, end-4+1):

		## for C atoms
		if j+5 < end+1:
			command = "assign (resid %3d and name  C) (resid %3d and name  C) 8.16 0.10 0.10 !helix" % (j,j+5)
			lines.append(command)
                command = "assign (resid %3d and name  C) (resid %3d and name  C) 4.87 0.10 0.10 !helix" % (j,j+3)
		lines.append(command)
                command = "assign (resid %3d and name  C) (resid %3d and name  C) 6.09 0.10 0.10 !helix" % (j,j+4)
		lines.append(command) 

		## for CA atoms
		if j+5 < end+1:
                	command = "assign (resid %3d and name CA) (resid %3d and name CA) 8.63 0.10 0.10 !helix" % (j,j+5)
			lines.append(command)
                command = "assign (resid %3d and name CA) (resid %3d and name CA) 5.13 0.10 0.10 !helix" % (j,j+3)
		lines.append(command)
                command = "assign (resid %3d and name CA) (resid %3d and name CA) 6.16 0.10 0.10 !helix" % (j,j+4)
		lines.append(command)

		## for N atoms
		if j+5 < end+1:
                	command = "assign (resid %3d and name  N) (resid %3d and name  N) 8.07 0.10 0.10 !helix" % (j,j+5)
			lines.append(command)
                command = "assign (resid %3d and name  N) (resid %3d and name  N) 4.84 0.10 0.10 !helix" % (j,j+3)
		lines.append(command)
                command = "assign (resid %3d and name  N) (resid %3d and name  N) 6.10 0.10 0.10 !helix" % (j,j+4)
		lines.append(command)

		## for O atoms
		if j+5 < end+1:
                	command = "assign (resid %3d and name  O) (resid %3d and name  O) 8.40 0.10 0.10 !helix" % (j,j+5)
			lines.append(command)
                command = "assign (resid %3d and name  O) (resid %3d and name  O) 4.99 0.10 0.10 !helix" % (j,j+3)
		lines.append(command)
                command = "assign (resid %3d and name  O) (resid %3d and name  O) 6.12 0.10 0.10 !helix" % (j,j+4)
		lines.append(command)

	return lines


## start and end the starting and ending positions of a beta strand
def Beta2DistRestraints(start, end, sequence):
	lines = []
	for j in range(start+1, end):
		command = "assign (resid %3d and name  O) (resid %3d and name  O) 4.57 0.10 0.10 !unpaired E residue" % (j, j+1)
		lines.append(command)

	return lines

def Helix2HBond(start, end, sequence):
	lines = []
	for j in range(start+1, end-4+1):
		## deal with a special amino acid
		if sequence[j+3] == 'P':
			command = "assign (resid %3d and name O) (resid %3d and name N) 2.99  0.1  0.1 !helix" % (j, j+4)
		else:
			command = "assign (resid %3d and name O) (resid %3d and name N) 1.99  0.1  0.1 !helix" % (j, j+4)

		lines.append(command)

	return lines

## derive NOE restraints and HBond restraints from predicted secondary structure
## in principle, ss3str can also be a ss8 string
def SS2Restraints(ss3str, sequence):
	SSSegments = []
	for m in re.finditer('(H+|E+)', ss3str):
		SSSegments.append( (m.start(), m.end(), m.group(0)[0]) )

	## each segment starts at m.start(), but ends at m.end()-1. That is, m.end()-1 is the index of the last position in a segment

	NOEs = []
	HBonds = []
	for seg in SSSegments:
		len = (seg[1] - seg[0])
		if seg[2] == 'H' and len>=5 :
			NOEs.extend(Helix2DistRestraints(seg[0], seg[1], sequence) )
			HBonds.extend(Helix2HBond(seg[0], seg[1], sequence) )
		elif seg[2] == 'E' and len>=2 :	
			NOEs.extend(Beta2DistRestraints(seg[0], seg[1], sequence) )
		elif seg[2] not in set(['H', 'E']):
			print 'ERROR: a secondary structure segment has type other than H or E: ', seg
			exit(1)

	return NOEs, HBonds

## derive angle restraints from predicted Phi, Psi angles
def PhiPsi2Restraints(PhiPsi, sequence, scale):
	size = PhiPsi.shape
	if len(sequence) != size[0] :
        	print 'ERROR: the sequence length does not match the PhiPsi size'
                exit(1)

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
		exit(1)

	return lines


def main(argv):

	if len(argv) < 2:
		Usage()
		exit(1)

	scale = 5.0

	try:
		opts, args = getopt.getopt(sys.argv[1:],"k:d:a:h:")
	except getopt.GetoptError:
		Usage()
		exit(1)

	NOEFile = 'ssnoe.tbl'
	PhiPsiFile = 'dihedral.tbl'
	HBondFile = 'hbond.tbl'

	for opt, arg in opts:
		if opt in ("-k", "--energyconstant"):
			scale = np.float32(arg)
		elif opt in ("-d", "--NOEFile"):
			NOEFile = arg
		elif opt in ("-a", "--AngleFile"):
			PhiPsiFile = arg
		elif opt in ("-h", "--HBondFile"):
			HBondFile = arg
		else:
			Usage()
			exit(1)

	if len(args)<2:
		print Usage()
		exit(1)

	

	seqFile = args[0]
	predictedPropertyFile = args[1]


	if not os.path.isfile(seqFile):
		print 'the sequence file does not exist: ', seqFile
		exit(1)

	if not os.path.isfile(predictedPropertyFile):
		print ' the predicted property file does not exist: ', predictedPropertyFile
		exit(1)


	sequence = LoadFASTAFile(seqFile)
	content = PropertyUtils.LoadPredictedProperties(predictedPropertyFile)

        name, sequence2, c1, c2 = content
	assert  (sequence == sequence2)

	if c2.has_key('SS3_Discrete3C'):
		ss3str = c2['SS3_Discrete3C']
	elif c2.has_key('SS8_Discrete8C'):
		ss3str, _ = PropertyUtils.SS8Prob2SS3(c1['SS8_Discrete8C'])
	else:
		print 'ERROR: the property file does not contain predicted secondary structure information: ', predictedPropertyFile
		exit(1)

	NOEs, HBonds = SS2Restraints(ss3str, sequence)
	

	if not c1.has_key('PhiPsi_vonMise2d4'):
		print 'ERROR: the property file does not contain predicted Phi, Psi information: ', predictedPropertyFile
		exit(1)

	PhiPsi = c1['PhiPsi_vonMise2d4']
	AngleRestraints = PhiPsi2Restraints(PhiPsi, sequence, scale)

	f=open(NOEFile, 'w')
	f.write('\n'.join(NOEs))
	f.close()

	f=open(HBondFile, 'w')
	f.write('\n'.join(HBonds))
	f.close()

	f=open(PhiPsiFile, 'w')
	f.write('\n'.join(AngleRestraints))
	f.close()


if __name__ == "__main__":
        main(sys.argv[1:])

