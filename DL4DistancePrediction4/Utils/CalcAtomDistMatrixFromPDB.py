import numpy as np
import sys
import os
import cPickle

from Bio.PDB import *
#from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import three_to_one

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import config
from SelectAtoms import SelectAtoms4Distance
from PDBUtils import ExtractPDBSeq

##this script calculate dist matrix for a protein from its PDB file.  

def Usage():
	print 'python CalcAtomDistMatrixFromPDB.py pdb_file [AtomPairTypes]'
	print '	      AtomPairTypes: the atom pairs to be used, default All'
	print '	      This script calculates atomic distance matrix from a pdb file and also outputs a sequence file'

def CalcAtomDistMatrix(pdbfile, apts=config.allAtomPairNames):

	name = os.path.basename(pdbfile).split('.')[0]
	parser = PDBParser()
	structure = parser.get_structure(name, pdbfile)
	residues = structure.get_residues()

	pdbseq, residueList = ExtractPDBSeq(residues)
	#print pdbseq

	seqLen = len(pdbseq)
	atomDistMatrix = dict()
	atomDistMatrix['name'] = name
        atomDistMatrix['pdbseq'] = pdbseq
        atomDistMatrix['seq4matrix'] = pdbseq

        for apt in apts:
                atomDistMatrix[apt] = -np.ones( (seqLen, seqLen), np.float16)

	numInvalidAtoms = dict()
        for atom in ['Ca', 'Cb', 'Cg', 'N', 'O']:
                numInvalidAtoms[atom]=0

	for i in range(seqLen):
		res1 = residueList[i]
		ca1, n1, o1, cb1, cg1 = SelectAtoms4Distance(res1)
	
		if ca1 is None:
                        numInvalidAtoms['Ca'] += 1
                if cb1 is None:
                        numInvalidAtoms['Cb'] += 1
                if cg1 is None:
                        numInvalidAtoms['Cg'] += 1
                if n1 is None:
                        numInvalidAtoms['N'] += 1
                if o1 is None:
                        numInvalidAtoms['O'] += 1
	
		for k in range(seqLen):
			res2 = residueList[k]
			ca2, n2, o2, cb2, cg2 = SelectAtoms4Distance(res2)

			for apt in apts:

                                ## Ca-Ca distance matrix, symmetric
                                if apt == 'CaCa' and (ca1 is not None) and (ca2 is not None):
                                        atomDistMatrix[apt][i, k] = ca1 - ca2

                                ## Cb-Cb, symmetric
                                if apt == 'CbCb' and (cb1 is not None) and (cb2 is not None):
                                        atomDistMatrix[apt][i, k] = cb1 - cb2

                                ## Cg-Cg, symmetric
                                if apt == 'CgCg' and (cg1 is not None) and (cg2 is not None):
                                        atomDistMatrix[apt][i, k] = cg1 - cg2

                                ## Ca-Cg, asymmetric
                                if apt == 'CaCg' and (ca1 is not None) and (cg2 is not None):
                                        atomDistMatrix[apt][i, k] = ca1 - cg2

                                ## N-O, asymmetric
                                if apt == 'NO' and (N is not None) and (O is not None):
                                        atomDistMatrix[apt][i, k]  = n1 - o2


	if atomDistMatrix.has_key('CaCa'):
		np.fill_diagonal(atomDistMatrix['CaCa'], 0)
	
	if atomDistMatrix.has_key('CbCb'):
		np.fill_diagonal(atomDistMatrix['CbCb'], 0)

	if atomDistMatrix.has_key('CgCg'):
		np.fill_diagonal(atomDistMatrix['CgCg'], 0)

	atomDistMatrix['numInvalidAtoms'] = numInvalidAtoms
        print 'numInvalidAtoms: ', numInvalidAtoms

	if np.any( np.fromiter(numInvalidAtoms.itervalues(), dtype=np.int32) >5 ):
                print 'ERROR: there are too many invalid atoms in ', pdbfile

	return atomDistMatrix


if len(sys.argv)<2:
	Usage()
	exit(1)	

pdbfile = sys.argv[1]

apts = config.allAtomPairNames
if len(sys.argv)>=3:
	apts = config.ParseAtomPairNames(sys.argv[2])

atomDistMatrix = CalcAtomDistMatrix(pdbfile, apts)
name = os.path.basename(pdbfile).split('.')[0]
savefile = name + '.atomDistMatrix.pkl'
fh = open(savefile, 'wb')
cPickle.dump(atomDistMatrix, fh, protocol=cPickle.HIGHEST_PROTOCOL)
fh.close()

savefile = name + '.fasta'
fh = open(savefile, 'w')
content = '\n'.join([ '>' + name, atomDistMatrix['pdbseq'] ])
fh.write(content)
fh.close()

