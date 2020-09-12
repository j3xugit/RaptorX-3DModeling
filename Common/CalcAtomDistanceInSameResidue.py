import os
import sys
import numpy as np
from PDBUtils import StructureParser
from SelectAtoms import SelectCG
from Bio.PDB.Polypeptide import three_to_one

def Usage():
	print 'python CalcAtomDistanceInSameResidue.py pdbfileList pdbfileFolder '
	print '	this script calculates the N-O distance in the same residue'

if len(sys.argv) < 3:
	Usage()
	exit(1)

pdbfileList = sys.argv[1]
pdbfileFolder = sys.argv[2]
NOdists = dict()
CaCgdists = dict()

with open(pdbfileList, 'r') as fh:
	proteins = [ line.strip() for line in list(fh) ]

for protein in proteins:
	pdbfile = os.path.join(pdbfileFolder, protein + '.pdb')
	if not os.path.isfile(pdbfile):
		print 'WARNING: pdbfile does not exist: ', pdbfile
		continue

	parser = StructureParser(pdbfile)
	structure = parser.get_structure('NoName', pdbfile)
	model = structure[0]

	for chain in model:
		for residue in chain:
			try:
				resName = residue.get_resname()
				"""
				n = residue['N']
				o = residue['O']
				NOdist = n-o
				if not NOdists.has_key(resName):
					NOdists[resName] = [NOdist ]
				else:
					NOdists[resName].append(NOdist)
				"""
				ca = residue['CA']
				AA = three_to_one(resName)
				cg = residue[SelectCG(AA)]
				AGdist = ca - cg
				if not CaCgdists.has_key(resName):
					CaCgdists[resName] = [AGdist ]
				else:
					CaCgdists[resName].append(AGdist)
				

			except:
				print 'WARNING: missing CA or CG atoms '

"""
finalNOdists = dict()
for res, dists in NOdists.iteritems():
	finalNOdists[res] = np.mean(dists)
	AA = three_to_one(res)
	finalNOdists[AA] = finalNOdists[res]

print finalNOdists
"""

finalAGdists = dict()
for res, dists in CaCgdists.iteritems():
	finalAGdists[res] = np.mean(dists)
	AA = three_to_one(res)
	finalAGdists[AA] = finalAGdists[res]
print finalAGdists
