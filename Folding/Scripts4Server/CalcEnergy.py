import numpy as np
import sys
import os
import cPickle

from Bio.PDB import *

pdbfile = sys.argv[1]
name = os.path.basename(pdbfile).split('.')[0]

potentialFile = sys.argv[2]
fh = open(potentialFile, 'rb')
potMatrix = cPickle.load(fh)
fh.close()

parser = PDBParser()

structure = parser.get_structure(name, pdbfile)

residues = structure.get_residues()

residueList = [ r for r in residues if is_aa(r) ]

numResidues = len(residueList)


size = potMatrix.shape

print 'number of residues in the PDB file: ', numResidues
print 'protein sequence length derived from the potential matrix: ', size[0]


if size[0] != numResidues:
	print 'Error: the number of residues in the PDB file is inconsistent with the size of the potential matrix'
	exit(-1)

energy = 0

for res1, i in zip(residueList, range(numResidues) ):

	if res1.get_resname() == 'GLY':
		a1 = res1['CA']
	else:
		a1 = res1['CB']

	for res2, j in zip(residueList[i+6: ], range(i+6, numResidues) ):
		if res2.get_resname() == 'GLY':
			a2 = res2['CA']
		else:
			a2 = res2['CB']

		distance = a1 - a2

		bin = 0
		if distance > 5:
			bin = np.int32(distance - 4)
		if bin > 11:
			bin = 11 

		if bin < 8:
			energy += potMatrix[i, j, bin]


print 'the energy of the protein structure for ' + name + ' is ', energy
