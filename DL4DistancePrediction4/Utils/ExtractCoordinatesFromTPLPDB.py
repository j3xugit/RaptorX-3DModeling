import numpy as np
import sys
import os
import cPickle

from Bio.PDB import *
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import three_to_one

sys.path.append(os.path.join(os.environ['ModelingHome'], 'Common'))
from LoadTPLTGT import load_tpl as LoadTPL

sys.path.append('../')
import config
from PDBUtils import ExtractSeqFromPDBFile
from SelectAtoms import SelectCG, SelectCB

## This script extract 3D coordinates for Cb, Ca, Cg, N amd O atoms for a template from a PDB file. 
## Two arguments are given: tpl file and PDB file.
## The tpl file may be inconsistent with the PDB file, but we will extract coordinates from PDB file based upon the tpl file

def Usage():
	print 'python ExtractCoordinatesFromPDB.py tpl_file pdb_file'

def MapTPL2PDB((DSSPsequence, SEQRESsequence):
	i=-1
	j=-1
	numMisMatches = 0

	mapping_seq2pdb = np.ones( (len(SEQRESsequence)), np.int32) * (-1)

	for a, b in zip(DSSPsequence, SEQRESsequence):
        	if a!='-' and b!='-':
                	i += 1
                	j += 1
                	mapping_seq2pdb[j] = i
                	numMisMatches += (a!=b)

        	elif a!='-':
                	i += 1
        	elif b!='-':
                	j += 1
        	else:
                	print 'ERROR:abnormal in the alignment between SEQRESsequence and DSSPsequence.'
                	exit(1)

	if numMisMatches > 0:
        	print 'ERROR:Some non-identical residues found between SEQRESsequence and DSSPsequence in ', tplfile
        	exit(1)


def ExtractCoordinatesFromTPLPDB(tplfile, pdbfile, atoms=['CB', 'CA', 'N', 'O']):
	tpl = LoadTPL(tplfile)
	SEQRESsequence = tpl['sequence']
	DSSPsequence = tpl['DSSPsequence']

	name = os.path.basename(tplfile).split('.')[0]
	if len(name)<5:
		print 'ERROR: the template name is incorrect. It must be composed of PDB ID and chain letter'
		exit(1)

	pdbseq, residueList = ExtractSeqFromPDBFile(pdbfile)
	#print pdbseq

	### check if DSSPsequence is equivalent to pdbseq
	validDSSPseq = DSSPsequence.replace('-', '')

	if validDSSPseq != pdbseq:
        	print 'ERROR: Inconsistency between DSSPsequence in ', tplfile, ' and pdbseq in ', pdbfile
        	exit(1)

	## here we simply use the sequence mapping in the TPL file
	mapping_seq2pdb = MapTPL2PDB(DSSPsequence, SEQRESsequence)

	## extract the coordinates for Ca, Cb, Cg, N and O
	atomCoordinates = []
	for i, j in zip(range(len(SEQRESsequence)), mapping_seq2pdb):
		coordinates = dict()
		if j < 0:
			assert (tpl['missing'][i]==1)
			for atom in atoms:
				coordinates[atom] = None
				coordinates['valid'] = False
			atomCoordinates.append(coordinates)
			continue

		if tpl['missing'][i] != 0:
			print 'Error: not a missing residue at index ', i
		exit(1)

		res = residueList[j]
		#print res
		for atom in atoms:
                        if atom.upper() == 'CG':
                                a = SelectCG(res.get_resname())
                                coordinates[atom] = res[a].get_vector()
                                continue

                        if atom.upper() == 'CB':
                                a = SelectCB(res.get_resname())
                                coordinates[atom] = res[a].get_vector()
                                continue

                        if res.has_id(atom):
                                coordinates[atom]= res[atom].get_vector()

                coordinates['valid'] = True
                atomCoordinates.append(coordinates)

		CaCoordinate = None
		for ca in ['CA', 'Ca', 'ca', 'cA']:
			if coordinates.has_key(ca):
				CaCoordinate = coordinates[ca]
				break
		if CaCoordinate is None:
			print 'ERROR: Ca atom is not included in atoms in ExtractCoordinatesFromTPLPDB'
			exit(1)

		if np.linalg.norm(Vector(tpl['Ca'][i]) - CaCoordinate) > 0.1 :
			print 'ERROR: maybe inconsistent Ca coordinates between the tpl file and PDB file for protein ', name
			print i
			print res1
			print tpl['Ca'][i]
			print coordinates['Ca']
			exit(1)

	return atomCoordinates

def main(argv):
	if len(argv)<2:
		Usage()
		exit(1)	

	tplfile = argv[0]
	pdbfile = argv[1]
	atomCoordinates = ExtractCoordinatesFromTPLPDB(tplfile, pdbfile)

	savefile = name + '.atomCoordinates.pkl'
	with open(savefile, 'wb') as fh:
		cPickle.dump(atomCoordinates, fh, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
        main(sys.argv[1:])

