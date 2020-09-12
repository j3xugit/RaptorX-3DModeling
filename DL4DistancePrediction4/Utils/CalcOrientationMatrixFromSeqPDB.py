import numpy as np
import sys
import os
import pickle

from Bio.PDB import *
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import three_to_one

from DL4DistancePrediction4.config import allOrientationNames as AllLabelNames
from Common.PDBUtils import ExtractSeqFromPDBFile
from Common.PDBUtils import MapSeq2PDB
from Common.SequenceUtils import LoadFASTAFile

'''
This script calculates 2 types of orientations, consisting of 4 dihedral and 3 angles
Meanwhile, there are 2 dihedral and 1 angle between two residues and 
2 dihedral and 2 angles between 4 Ca atoms: Ca(i), Ca(i+1), Ca(j), Ca(j+1)

Two arguments are needed: seq file (optional), pdb file. 
seq file may be inconsistent with pdb file, so we shall align sequence to PDB first
'''

NoCB = -np.float16(540)
ValueOfSelf = 0
NoValidCoordinates = -np.float16(720)

## a simple check. Exclude GLY from consideration
def HasAtoms4TwoROri(residue1, residue2):
	if residue1.get_resname() == 'GLY' or residue2.get_resname() == 'GLY':
		return True

	if residue1.has_id('CA') and residue1.has_id('CB') and residue1.has_id('N') and residue2.has_id('CA') and residue2.has_id('CB'):
		return True
	else:
		return False

def CalcOrientationOf2Residues(residue1, residue2):
    	ca1_cb1_cb2_ca2 = NoCB
    	n1_ca1_cb1_cb2 = NoCB
    	ca1_cb1_cb2 = NoCB
    
    	if 'CB' in residue1 and 'CB' in residue2:
		## both Cb atoms exist
		if residue1.has_id('CA') and residue2.has_id('CA'):
        		ca1_cb1_cb2_ca2 = calc_dihedral(residue1['CA'].get_vector(), residue1['CB'].get_vector() ,residue2['CB'].get_vector(), residue2['CA'].get_vector())*180/np.pi
		else:
			ca1_cb1_cb2_ca2 = NoValidCoordinates

		if residue1.has_id('N') and residue1.has_id('CA'):
        		n1_ca1_cb1_cb2 = calc_dihedral(residue1['N'].get_vector(), residue1['CA'].get_vector() ,residue1['CB'].get_vector(), residue2['CB'].get_vector())*180/np.pi
		else:
			n1_ca1_cb1_cb2 = NoValidCoordinates

		if residue1.has_id('CA'):
        		ca1_cb1_cb2 = calc_angle(residue1['CA'].get_vector() ,residue1['CB'].get_vector(), residue2['CB'].get_vector())*180/np.pi
		else:
			ca1_cb1_cb2 = NoValidCoordinates

    	elif 'CB' in residue1:
		## the 2nd Cb atom does not exist, we use Ca instead
		if residue1.has_id('N') and residue1.has_id('CA') and residue2.has_id('CA'):
        		n1_ca1_cb1_cb2 = calc_dihedral(residue1['N'].get_vector(), residue1['CA'].get_vector() ,residue1['CB'].get_vector(), residue2['CA'].get_vector())*180/np.pi
		else:
			n1_ca1_cb1_cb2 = NoValidCoordinates

		if residue1.has_id('CA') and residue2.has_id('CA'):
        		ca1_cb1_cb2 = calc_angle(residue1['CA'].get_vector() ,residue1['CB'].get_vector(), residue2['CA'].get_vector())*180/np.pi
		else:
			ca1_cb1_cb2 = NoValidCoordinates

        
    	return {'Ca1Cb1Cb2Ca2':ca1_cb1_cb2_ca2, 'N1Ca1Cb1Cb2':n1_ca1_cb1_cb2, 'Ca1Cb1Cb2':ca1_cb1_cb2 }


## calculate the angle and dihedral of 4 Ca atoms, i1 and i2 represent two adjacent sequence positions and j1 and j2 represent another two adjacent positions
def CalcOrientationOf4CAs(CAi1, CAi2, CAj1, CAj2):
	dihedral = calc_dihedral(CAi1.get_vector(), CAi2.get_vector() ,CAj1.get_vector(), CAj2.get_vector())*180/np.pi
        angle = calc_angle(CAi1.get_vector() , CAi2.get_vector(), CAj1.get_vector())*180/np.pi
	dihedral2 = calc_dihedral(CAi1.get_vector(), CAi2.get_vector() ,CAj2.get_vector(), CAj1.get_vector())*180/np.pi
        angle2 = calc_angle(CAi1.get_vector() , CAi2.get_vector(), CAj2.get_vector())*180/np.pi

	return {'Ca1Ca2Ca3Ca4': dihedral, 'Ca1Ca2Ca3': angle, 'Ca1Ca2Ca4Ca3': dihedral2, 'Ca1Ca2Ca4': angle2}

## if sequence is not None, then it will be used
## otherwise, if seq_file is not None, then load sequence from this file
## otherwise, the sequence info in the pdb file will be used
def CalcOrientationMatrix(pdbfile, sequence=None, seq_file=None):
    
	## get pdb residues and pdb sequence 
    	structure_id = os.path.splitext(os.path.basename(pdbfile))[0]
    	pdbseq, residueList = ExtractSeqFromPDBFile(pdbfile, name=structure_id)
    
    	## get sequence-pdbseq map
    	mapping_seq2pdb = None
    	if sequence is None and seq_file is None:
        	sequence = pdbseq
        	mapping_seq2pdb = np.arange(len(pdbseq))
    	elif sequence is None:
        	sequence = LoadFASTAFile(seq_file)        
       	mapping_seq2pdb = MapSeq2PDB(sequence, pdbseq, residueList)
    	if mapping_seq2pdb is None:
        	print('ERROR: cannot map the sequence information to a pdb file: ' + pdbfile)
        	exit(1)
    
    	##calculate the orientation matrices
    	OrientationMatrix = {}
    	OrientationMatrix['pdbseq'] = pdbseq
    	OrientationMatrix['seq4matrix'] = str(sequence)
    
    	seqLen = len(sequence)
    	for apt in AllLabelNames:
        	OrientationMatrix[apt] = np.ones( (seqLen, seqLen), np.float16) * NoValidCoordinates

	numInvalidAtoms=dict()
	numInvalidAtoms['TwoROri']=0
	numInvalidAtoms['FourCaOri']=0    

    	for i, j in zip(range(seqLen), mapping_seq2pdb):
        	if j<0:
            		continue
        	residue_i = residueList[j]
        	for k, l in zip(range(seqLen), mapping_seq2pdb):
            		if l<0 :
                		continue
            		if i==k:
                		for apt in AllLabelNames:
                			OrientationMatrix[apt][i, k] = ValueOfSelf
				continue

            		residue_k = residueList[l]
			if not HasAtoms4TwoROri(residue_i, residue_k):
				numInvalidAtoms['TwoROri'] += 1

            		pairOri = CalcOrientationOf2Residues(residue_i, residue_k)
            		for apt, v in pairOri.iteritems():
            			OrientationMatrix[apt][i, k] = v

	    		if i == seqLen-1 or k==seqLen-1:
				continue

	    		j1 = mapping_seq2pdb[i+1]
	    		l1 = mapping_seq2pdb[k+1]
	    		if j1<0 or l1<0:
				continue
	    		residue_i1 = residueList[j1]
	    		residue_k1 = residueList[l1]

			if not residue_i.has_id('CA') or not residue_i1.has_id('CA') or not residue_k.has_id('CA') or not residue_k1.has_id('CA'):
				numInvalidAtoms['FourCaOri'] += 1
			else:
	    			CAi = residue_i['CA']
	    			CAi1 = residue_i1['CA']
	    			CAk = residue_k['CA']
	    			CAk1 = 	residue_k1['CA']

	    			CaOri = CalcOrientationOf4CAs(CAi, CAi1, CAk, CAk1)
	    			for apt, v in CaOri.iteritems():
	    				OrientationMatrix[apt][i, k] = v

    	for apt in AllLabelNames:
		OrientationMatrix[apt] = (OrientationMatrix[apt]).astype(np.float16)

	OrientationMatrix['numInvalidAtoms'] = numInvalidAtoms
	print('numInvalidAtoms: ' + str(numInvalidAtoms) )

	if np.any( np.fromiter(numInvalidAtoms.itervalues(), dtype=np.int32) >5 ):
                print('ERROR: there are too many invalid atoms in '+pdbfile)
	
    	return OrientationMatrix


def Usage(script_name):
    	print('Usage: '+script_name+' pdb_file [seq_file]')
    	print('  1) When seq_file not provided, the sequence consisting of all residues with valid 3D coordinates is used')
    	print("  2) A valid orientation has value ranging from -180 to 180")
    	print("  3) The value is set to 540 when Cb atom does not exist")
    	print("  4) The value along the diagonal is set to 0")
    	print("  5) The value is set to 720 when a residue does not have valid 3D coordinates")

if __name__ == '__main__':
	#argv = [sys.argv[0],'example/1ctfA.pdb', 'example/1ctfA.fasta']
    	argv = sys.argv
    	if len(argv)<2 or len(argv)>3:
        	Usage(argv[0])
        	exit(1)

    	pdb_file = argv[1]

    	seq_file = None
    	if len(argv) ==3:
        	seq_file = argv[2]

    	print(argv[0],pdb_file, seq_file if seq_file else '')
    	if not os.path.exists(pdb_file):
        	print(pdb_file+' not found.')
        	exit(1)
    	if seq_file is not None and not os.path.exists(seq_file):
        	print(seq_file+' not found.')
        	exit(1)
    	oriMatrix = CalcOrientationMatrix(pdb_file, seq_file)

    	name = os.path.splitext(os.path.basename(pdb_file))[0]
    	oriMatrix['name'] = name
    	savefile = name + '.atomOrientationMatrix.pkl'
    	with open(savefile, 'wb') as pkl_file:
        	pickle.dump(oriMatrix, pkl_file, protocol=2)

