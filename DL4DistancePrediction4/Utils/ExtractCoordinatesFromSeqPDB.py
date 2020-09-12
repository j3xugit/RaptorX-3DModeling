import sys
import os
import numpy as np
import cPickle

from Bio.PDB import *
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import three_to_one

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

from PDBUtils import ExtractCoordinatesBySeq
from SequenceUtils import LoadFASTAFile

## This script extract 3D coordinates for Cb, Ca, Cg, N amd O atoms for a template from a PDB file. 
## Two arguments are given: seq file and PDB file.
## The seq file may be inconsistent with the PDB file. In this case, seq will be mapped to PDB and then corresponding coordinates will be extracted

## When an amino acid (Glycine) does not have Cb atom, Ca will be used. If you do not want this, please set bUseAlternativeAtoms=False in ExtractCoordinates
## When some amino acids do not have CG atom, Cb or Ca will be used. If do not want this, please set bUseAlternativeAtoms=False

def Usage():
	print 'python ExtractCoordinatesFromSeqPDB.py seq_file pdb_file [numAllowedMisMatches]'
	print '	numAllowedMisMatches: the max number (default 5) of allowed mismatches in the aligned positions of query sequence and pdb sequence '


def main(argv):
	if len(argv) < 2:
		Usage()
		exit(1)

	seqfile = argv[0]
	pdbfile = argv[1]
	maxMisMatches = 5

	if len(argv) >= 3:
		maxMisMatches = np.int32(argv[2])

	sequence = LoadFASTAFile(seqfile)

	name = os.path.basename(seqfile).split('.')[0]
	atomCoordinates, pdbseq, numMisMatches = ExtractCoordinatesBySeq(sequence, pdbfile, maxMisMatches=maxMisMatches)
	if atomCoordinates is None:
		print 'ERROR: failed to extract atom coordinates from pdb file: ', pdbfile
		exit(1)

	result = dict()
	result['name'] = name
	result['sequence'] = sequence
	result['pdbseq'] = pdbseq
	result['numMisMatches'] = numMisMatches
	result['coordinates'] = atomCoordinates

	savefile = name + '.atomCoordinates.pkl'
	with open(savefile, 'wb') as fh:
		cPickle.dump(result, fh, protocol=cPickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
        main(sys.argv[1:])

