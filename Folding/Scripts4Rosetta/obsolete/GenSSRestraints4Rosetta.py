import numpy as np
import sys
import os
import getopt
import cPickle
import re

from DL4PropertyPrediction import config, PropertyUtils
from Common.LoadFASTA import LoadFASTAFile

def Usage():
	print 'python GenSSRestraints4Rosetta.py [-k standard deviation ] seqfile_in_FASTA predictedProperties_in_PKL\n'
	print '   This script generates secondary structure constraints using Rosetta HARMONIC function'
	print '   -k: standard deviation for harmonic function, default 0.2'
	print '   results will be saved to a filenamed after targetName.SS.kStd.txt'

def ConstraintStr(atomType, res1, res2, meanDist, k=0.2):
	resStr = "AtomPair %s %4d %s %4d HARMONIC %f %f" % (atomType, res1, atomType, res2, meanDist, k)
	return resStr

def ConstraintStr2(atomType1, atomType2, res1, res2, meanDist, k=0.2):
	resStr = "AtomPair %s %4d %s %4d HARMONIC %f %f" % (atomType1, res1, atomType2, res2, meanDist, k)
	return resStr

## start and end: the starting and ending positions of a helix segment, counting from 0
## in Rosetta, the residue index counts from 1, i.e., we shall add 1 to each position in a segment
def Helix2DistRestraints(start, end, sequence, k=0.2):

	lines = []
	for j in range(start+1, end-4+1):

		## for C atoms
		if j+5 < end+1:
			lines.append(ConstraintStr('C', j, j+5, 8.16, k) )
		lines.append(ConstraintStr('C', j, j+3, 4.87, k) )
		lines.append(ConstraintStr('C', j, j+4, 6.09, k) )

		## for CA atoms
		if j+5 < end+1:
			lines.append(ConstraintStr('CA', j, j+5, 8.63, k) )
		lines.append(ConstraintStr('CA', j, j+3, 5.13, k) )
		lines.append(ConstraintStr('CA', j, j+4, 6.16, k) )

		## for N atoms
		if j+5 < end+1:
			lines.append(ConstraintStr('N', j, j+5, 8.07, k) )
		lines.append(ConstraintStr('N', j, j+3, 4.84, k) )
		lines.append(ConstraintStr('N', j, j+4, 6.10, k) )

		## for O atoms
		if j+5 < end+1:
			lines.append(ConstraintStr('O', j, j+5, 8.40, k) )
		lines.append(ConstraintStr('O', j, j+3, 4.99, k) )
		lines.append(ConstraintStr('O', j, j+4, 6.12, k) )

	return lines


## start and end: the starting and ending positions of a beta strand
def Beta2DistRestraints(start, end, sequence, k=0.2):
	lines = []
	for j in range(start+1, end):
		command = ConstraintStr('O', j, j+1, 4.57, k) 
		lines.append(command)

	return lines

def Helix2HBond(start, end, sequence, k=0.2):
	lines = []
	for j in range(start+1, end-4+1):
		## deal with a special amino acid
		if sequence[j+3] == 'P':
			command = ConstraintStr2('O', 'N', j, j+4, 2.99, k)
		else:
			command = ConstraintStr2('O', 'N', j, j+4, 1.99, k)

		lines.append(command)

	return lines

## derive NOE restraints and HBond restraints from predicted secondary structure
## in principle, ss3str can also be a ss8 string
def SS2Restraints(ss3str, sequence, k=0.2):
	SSSegments = []
	for m in re.finditer('(H+|E+)', ss3str):
		SSSegments.append( (m.start(), m.end(), m.group(0)[0]) )

	## each segment starts at m.start(), but ends at m.end()-1. That is, m.end()-1 is the index of the last position in a segment

	NOEs = []
	HBonds = []
	for seg in SSSegments:
		len = (seg[1] - seg[0])
		if seg[2] == 'H' and len>=5 :
			NOEs.extend(Helix2DistRestraints(seg[0], seg[1], sequence, k=k) )
			HBonds.extend(Helix2HBond(seg[0], seg[1], sequence, k=k) )
		elif seg[2] == 'E' and len>=2 :	
			NOEs.extend(Beta2DistRestraints(seg[0], seg[1], sequence, k=k) )
		elif seg[2] not in set(['H', 'E']):
			print 'ERROR: a secondary structure segment has type other than H or E: ', seg
			exit(1)

	return NOEs, HBonds


def main(argv):

	if len(argv) < 2:
		Usage()
		exit(1)
	k = 0.2

	try:
		opts, args = getopt.getopt(sys.argv[1:],"k:")
	except getopt.GetoptError:
		Usage()
		exit(1)

	for opt, arg in opts:
		if opt in ("-k", "--energyconstant"):
			k = np.float32(arg)
		else:
			Usage()
			exit(1)

	if len(args)<2:
		print Usage()
		exit(1)

	seqFile = args[0]
	predictedPropertyFile = args[1]

	if not os.path.isfile(seqFile):
		print 'ERROR: the sequence file does not exist: ', seqFile
		exit(1)

	if not os.path.isfile(predictedPropertyFile):
		print 'ERROR: the predicted property file does not exist: ', predictedPropertyFile
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

	NOEs, HBonds = SS2Restraints(ss3str, sequence, k=k)

	savefile = name + '.SS' + '.k' + str(k) + '.txt'
	f=open(savefile, 'w')
	f.write('\n'.join(NOEs + HBonds))
	f.close()


if __name__ == "__main__":
        main(sys.argv[1:])

