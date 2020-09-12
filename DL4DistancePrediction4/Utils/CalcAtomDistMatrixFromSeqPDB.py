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
#sys.path.append('../')

import config
from PDBUtils import MapSeq2PDB, ExtractPDBSeq, ExtractSeqFromPDBFile
from SelectAtoms import SelectAtoms4Distance

##this script calculate dist matrix for a protein. Two arguments are given: seq file and PDB file. 
## One optional argument: UseChainName. If this option is specified, then get a chain name from the protein name. Default: not use chain name
## seq file may be inconsistent with PDB file, but we will calculate the distance matrix based upon seq file by extracting coordinates from PDB file

def Usage():
	print 'python CalcAtomDistMatrixFromSeqPDB.py seq_file_in_FASTA pdb_file [AtomPairTypes]'
	print '	  this script calculates the atomic distance matrix from sequence and PDB files'
	print '   AtomPairTypesied: the atom pairs to be used, default All'


## apts is a list or set, containing allowed atom pair types
def CalcAtomDistMatrix(sequence, pdbfile, name='NoName', apts=config.allAtomPairNames, maxMisMatches=5):
	seqLen = len(sequence)
	pdbseq, residueList = ExtractSeqFromPDBFile(pdbfile)
	mapping_seq2pdb = MapSeq2PDB(sequence, pdbseq, residueList)

	##calculate the distance matrix
	atomDistMatrix = dict()
	atomDistMatrix['pdbseq'] = pdbseq
	atomDistMatrix['seq4matrix'] = str(sequence)

	for apt in apts:
		atomDistMatrix[apt] = -np.ones( (seqLen, seqLen), np.float16)

	numInvalidAtoms = dict()
	for atom in ['Ca', 'Cb', 'Cg', 'N', 'O']:
		numInvalidAtoms[atom]=0

	for i, j in zip(range(seqLen), mapping_seq2pdb):
		if j < 0:
			continue

		res1 = residueList[j]
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

		for k, l in zip(range(seqLen), mapping_seq2pdb):
			if l <0 :
				continue

			res2 = residueList[l]
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
				if apt == 'NO' and (n1 is not None) and (o2 is not None):
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

def main(argv):
	if len(argv)<2:
		Usage()
		exit(1)	

	seqfile = argv[0]
	pdbfile = argv[1]

	print seqfile, pdbfile

	apts = config.allAtomPairNames

	if len(argv)>=3:
		apts = config.ParseAtomPairNames(argv[2])

	name = os.path.basename(seqfile).split('.')[0]

	from Bio import SeqIO
	record = SeqIO.read(seqfile, "fasta")
	sequence = record.seq

	#atomDistMatrix = CalcAtomDistMatrix(sequence, pdbfile, apts=apts, maxMisMatches=10)
	atomDistMatrix = CalcAtomDistMatrix(sequence, pdbfile, apts=apts)

	atomDistMatrix['name'] = name
	savefile = name + '.atomDistMatrix.pkl'
	fh = open(savefile, 'wb')
	cPickle.dump(atomDistMatrix, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	fh.close()


if __name__ == '__main__':
        main(sys.argv[1:])


