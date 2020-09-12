import os
import sys
import cPickle

from SequenceUtils import WriteFASTAFile
import PDBUtils

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print 'python CalcGroundTruthFromPDBChain.py pdbFile [chainName]'
		print '	this script calculates protein properties, distance and orientation matrices from a PDB file and a chainName'
		print '	chainName is used to determine which chain shall be selected from PDB file, default A'
		print ' save the result to a cPickle file and the corresponding sequence to a FASTA file'
		print ' the result files are named after PDBID_chainName.native.pkl and PDBID_chainName.fasta'
		exit(1)

	pdbfile = sys.argv[1]
	chainName = 'A'
	if len(sys.argv)>=3:
		chainName = sys.argv[2]

	name = os.path.basename(pdbfile).split('.')[0] + chainName
	sequence = PDBUtils.GetSEQRESOfOneChain(pdbfile, chainName)
	if sequence is None:
		print 'ERROR: cannot find SEQRES for chain ', chainName, ' in pdbfile: ', pdbfile
		exit(1)

	pdbseq, residueList, chain = PDBUtils.GetOneChain(pdbfile, chainName)
	if pdbseq is None:
		print 'ERROR: cannot find chain ', chainName, ' in the pdbfile: ', pdbfile
		exit(1)

	seq2pdb_mapping, numMisMatches, numMatches = PDBUtils.MapSeq2ResidueList(sequence, pdbseq, residueList)
	if seq2pdb_mapping is None:
		print 'ERROR: cannot map from sequence to pdb residues: ', pdbfile
		exit(1)

	if numMisMatches > 5:
		print 'ERROR: two many mismatches between SEQRES and ATOM record in ', pdbfile
		exit(1)

	if numMatches < min(30, 0.5*len(sequence)):
		print 'ERROR: more than half of SEQRES not covered by ATOM record in ', pdbfile
		exit(1)

	protein = dict()

	protein['name']=name
	protein['sequence'] = sequence	
	protein['numMisMatches'] = numMisMatches
	protein['numMatches'] = numMatches
	protein['pdbseq'] = pdbseq

	dssp = PDBUtils.ExtractDSSPByMapping(pdbfile, sequence, pdbseq, seq2pdb_mapping, residueList, chain)

	protein.update(dssp)

	coordinates, numInvalidAtoms = PDBUtils.ExtractCoordinatesByMapping(sequence, seq2pdb_mapping, residueList)
	if numInvalidAtoms.has_key('CA') and numInvalidAtoms['CA']>10:
		print 'ERROR: too many Ca atoms do not have valid 3D coordinates in chain ', chainName, ' of ', pdbfile
		exit(1)
	if numInvalidAtoms.has_key('CB') and numInvalidAtoms['CB']>10:
		print 'ERROR: too many Cb atoms do not have valid 3D coordinates in chain ', chainName, ' of ', pdbfile
		exit(1)
	
	protein['seq4matrix'] = sequence	
	protein['numInvalidAtoms'] = numInvalidAtoms

	protein['missing'] = [ c is None or (c['CA'] is None and c['CB'] is None) for c in coordinates ]

	distMatrix = PDBUtils.CalcDistMatrix(coordinates)
	protein['atomDistMatrix'] = PDBUtils.PostProcessDistMatrix(sequence, distMatrix)

	oriMatrix = PDBUtils.CalcTwoROriMatrix(coordinates)
	protein['atomOrientationMatrix'] = oriMatrix

	oriMatrix = PDBUtils.CalcCaOriMatrix(coordinates)
	protein['atomOrientationMatrix'].update(oriMatrix)

	savefile = protein['name'] + '.native.pkl'
	with open(savefile, 'w') as fh:
		cPickle.dump(protein, fh, protocol=cPickle.HIGHEST_PROTOCOL)

	## also write a FASTA file
	savefile = protein['name'] + '.fasta'
	WriteFASTAFile(name, sequence, savefile)
