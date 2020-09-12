import os
import sys
import numpy as np
import scipy
from scipy import spatial

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.vectors import calc_dihedral, calc_angle, rotaxis
#from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import three_to_one, is_aa
from Bio.SeqIO.PdbIO import PdbSeqresIterator
import Bio.PDB.Dice
from Bio.PDB.DSSP import DSSP
from Bio.PDB.PDBExceptions import PDBException

from Bio import pairwise2, SeqIO
from Bio.SubsMat.MatrixInfo import blosum80

from SelectAtoms import SelectCG, SelectCB
from SequenceUtils import ValidAA3Letters
from SSUtils import SS8Letter2SS3Letter

## the distance of N and O atoms in the same residue
NOdists = {'X': 3.06, 'ILE': 3.174216, 'GLN': 3.3115633, 'GLY': 3.2121522, 'GLU': 3.3253117, 'L': 3.2752273, 'CYS': 3.1867874, 'HIS': 3.2300916, 'G': 3.2121522, 'LYS': 3.2851896, 'PRO': 3.1523004, 'ASN': 3.2947414, 'A': 3.302093, 'C': 3.1867874, 'E': 3.3253117, 'D': 3.2847655, 'VAL': 3.1252666, 'F': 3.2107856, 'I': 3.174216, 'H': 3.2300916, 'K': 3.2851896, 'M': 3.2734911, 'THR': 3.1898365, 'N': 3.2947414, 'Q': 3.3115633, 'P': 3.1523004, 'S': 3.2022767, 'ASP': 3.2847655, 'T': 3.1898365, 'W': 3.1940844, 'V': 3.1252666, 'Y': 3.1847978, 'TRP': 3.1940844, 'SER': 3.2022767, 'PHE': 3.2107856, 'ALA': 3.302093, 'MET': 3.2734911, 'R': 3.267973, 'LEU': 3.2752273, 'ARG': 3.267973, 'TYR': 3.1847978}

# the distance of Ca and Cg atoms in the same residue. for some amino acids without Cg atoms, Cb or Ca atoms are used for Cg
CaCgdists = {'X': 2.2, 'ILE': 2.529928, 'GLN': 2.5625532, 'GLY': 0.0, 'GLU': 2.562029, 'L': 2.6051538, 'CYS': 2.8013322, 'HIS': 2.5359476, 'G': 0.0, 'LYS': 2.564922, 'PRO': 2.3897746, 'ASN': 2.539332, 'A': 1.5232254, 'C': 2.8013322, 'E': 2.562029, 'D': 2.5434172, 'VAL': 2.527417, 'F': 2.5461287, 'I': 2.529928, 'H': 2.5359476, 'K': 2.564922, 'M': 2.559007, 'THR': 2.5252295, 'N': 2.539332, 'Q': 2.5625532, 'P': 2.3897746, 'S': 2.4277365, 'ASP': 2.5434172, 'T': 2.5252295, 'W': 2.543401, 'V': 2.527417, 'Y': 2.552031, 'TRP': 2.543401, 'SER': 2.4277365, 'PHE': 2.5461287, 'ALA': 1.5232254, 'MET': 2.559007, 'R': 2.5622983, 'LEU': 2.6051538, 'ARG': 2.5622983, 'TYR': 2.552031}

# the distance of two sequentially adjacent Ca atoms
AdjacentCaCadist = 3.8

ValueOfSelf = 0
InvalidDistance = -1
NoCBDegree = np.float16(-540)
InvalidDegree = np.float16(-720)

def StructureParser(pdbfile):
	if pdbfile.endswith('.pdb'):
                parser = PDBParser(QUIET=True)
        elif pdbfile.endswith('.cif'):
                parser = MMCIFParser(QUIET=True)
        else:
                print 'ERROR: a protein structure file shall end with either .pdb or .cif'
                exit(1)
	return parser

def BuildVirtualCB(residue):
	if residue.get_resname().upper() != 'GLY':
		print 'ERROR: you shall not build a virtual Cb atom for non-Glycine residue'
		exit(1)

	# get atom coordinates as vectors
	try:
		n = residue['N'].get_vector()
		c = residue['C'].get_vector()
		ca = residue['CA'].get_vector()
	except:
		print 'WARNING: failed to obtain all of N, C and CA atoms in building a virtual Cb atom for residue ', residue.get_full_id()
		return None

	# center at origin
	n = n - ca
	c = c - ca

	# find rotation matrix that rotates n -120 degrees along the ca-c vector
	rot = rotaxis(-np.pi*120.0/180.0, c)

	# apply rotation to ca-n vector
	cb_at_origin = n.left_multiply(rot)

	# put on top of ca atom
	cb = cb_at_origin + ca
	return cb

## get the SEQRES of one specific chain. returnStr determines if a python string shall be returned (default)
def GetSEQRESOfOneChain(pdbfile, chainName, returnStr=True):
        with open(pdbfile) as fh:           
                for record in PdbSeqresIterator(fh):
                        if record.annotations['chain'] == chainName:
				if returnStr:
					return str(record.seq)
				else:
                                	return record.seq 
        return None

## extract a list of residues with valid 3D coordinates excluding non-standard amino acids
## returns the amino acid sequence as well as a list of residues with standard amino acids
def ExtractPDBSeq(residues):
	residueList = [ r for r in residues if is_aa(r, standard=True) and (r.get_resname().upper() in ValidAA3Letters) ]
	#print residueList
    	pdbseq = ''.join( [ three_to_one(r.get_resname()) for r in residueList ] )

    	return pdbseq, residueList

## get one specific chain including amino acid seq and residues with valid 3D coordinates in the pdbfile
def GetOneChain(pdbfile, chainName):
	parser = StructureParser(pdbfile)
        structure = parser.get_structure('NoName', pdbfile)
	model = structure[0]
	
	for chain in model:
		if chain.get_id() != chainName:
			continue
		residues = chain.get_residues()
		pdbseq, residueList = ExtractPDBSeq(residues)
		return pdbseq, residueList, chain

	return None, None, None	

## this function extracts one PDB chain and saves it to a file
def SaveOneChain2File(pdbfile, chainName, savefile=None):
	savefilesuffix = '.pdb'
	if pdbfile.endswith('.cif'):
		savefilesuffix = '.cif'

	parser = StructureParser(pdbfile)
	structure = parser.get_structure('NoName', pdbfile)
        model = structure[0]

	existing = False
        for chain in model:
                if chain.get_id() == chainName:
			existing = True
                        break
	if not existing:
		print 'ERROR: the specified chain ', chainName, ' does not exist in pdbfile: ', pdbfile
		return

	start = np.iinfo(np.int32).min
	end = np.iinfo(np.int32).max

	if savefile is None:
		file4save = os.path.basename(pdbfile).split('.')[0] + chainName + savefilesuffix
	else:
		file4save = savefile
	Bio.PDB.Dice.extract(structure, chainName, start, end, file4save)

## extract sequences and residue lists for each chain
## return pdbseqs, residue lists and also the chain objects
def ExtractPDBSeqByChain(structure):
	model=structure[0]
	pdbseqs = []
	residueLists = []
	chains = []
	for chain in model:
		residues = chain.get_residues()
		pdbseq, residueList = ExtractPDBSeq(residues)
		pdbseqs.append(pdbseq)
		residueLists.append(residueList)
		chains.append(chain)

	return pdbseqs, residueLists, chains

## extract sequences and residue lists from pdbfile for all the chains
def ExtractSeqFromPDBFile(pdbfile, name=None):
	#print 'In ExtractSeqFromPDBFile: ', pdbfile
	parser = StructureParser(pdbfile)
	if name is None:
		tmpName = os.path.basename(pdbfile)
	else:
		tmpName = name
	structure = parser.get_structure(tmpName, pdbfile)
	return ExtractPDBSeqByChain(structure)

## calculate the index mapping from first seq to the 2nd one, using the provided alignment.
## return a list with the same length as the first sequence (not the first enetry in alignment since alignment may contains gaps)
def Alignment2Mapping(alignment):
	S1, S2 = alignment[:2]

	## convert an aligned seq to a binary vector with 1 indicates aligned and 0 gap
	y = np.array([ 1 if a!='-' else 0 for a in S2 ])

	## get the position of each residue in the original sequence, starting from 0.
	ycs = np.cumsum(y) - 1
	np.putmask(ycs, y==0, -1)

	## map from the 1st seq to the 2nd one. set -1 for an unaligned residue in the 1st sequence
	mapping = [ y0 for a, y0 in zip(S1, ycs) if a!='-' ]

	return mapping

## calculate the number of mismatches and matches excluding residues denoted as 'X'
## regardless of the orders of the two sequences in the alignment, this function returns the same values
def CalcNumMisMatches(alignment):
	S1, S2 = alignment[:2]
	numMisMatches = np.sum( [ a!=b for a, b in zip(S1, S2) if a!='-' and b!='-' and a!='X' and b!='X' ] )
	numMatches = np.sum( [ a==b for a, b in zip(S1, S2) if a!='-' and a!='X' ] )

	return numMisMatches, numMatches
	
## map one query sequence to a list of PDB residues by sequence alignment
## pdbseq and residueList are generated by ExtractPDBSeq or ExtractPDBSeqByChain from a PDB file
## returns seq2pdb mapping, numMisMatches and numMatches
def MapSeq2ResidueList(sequence, pdbseq, residueList):

	## here we align PDB residues to query sequence instead of query to PDB residues
    	alignments = pairwise2.align.localds(pdbseq, sequence, blosum80, -5, -0.2)
	#print('#num alignments:', len(alignments))
	if not bool(alignments):
		return None, None, None

    	##find the alignment with the minimum difference
    	diffs = []
    	for alignment in alignments:
		mapping_pdb2seq = Alignment2Mapping(alignment)
        	diff = 0
		for current_map, prev_map, current_residue, prev_residue in zip( mapping_pdb2seq[1:], mapping_pdb2seq[:-1], residueList[1:], residueList[:-1] ):
			## in principle, every PDB residue with valid 3D coordinates shall appear in the query sequence.
			## otherwise, apply a big penalty
			if current_map<0:
				diff +=10
				continue

			if prev_map <0:
				continue

			## calculate the difference of sequence separation in both the PDB seq and the query seq
			## the smaller, the better
			current_id = current_residue.get_id()[1]
			prev_id = prev_residue.get_id()[1]
			id_diff = max(1, current_id - prev_id)
			map_diff = current_map - prev_map
			diff += abs( id_diff - map_diff)

		numMisMatches, numMatches = CalcNumMisMatches(alignment)
        	diffs.append(diff - numMatches)

    	diffs = np.array(diffs)
    	#print('best alignment:', diffs.argmin())
    	alignment = alignments[diffs.argmin()]

    	##print alignments
    	##print(alignment[1])
    	##print(alignment[0])

	numMisMatches, numMatches = CalcNumMisMatches(alignment)
	#print 'numMisMatches: ', numMisMatches
	#print 'numMatches: ', numMatches

	## map from the query seq to pdb
	mapping_seq2pdb = Alignment2Mapping( (alignment[1], alignment[0]) )

    	return mapping_seq2pdb, numMisMatches, numMatches

## the last two parameters are the maximum allowed mismatches and the minimum ratio of matches on the query sequence
## if a good match is identified, 
## returns seq2pdb mapping, the pdb residue list, the pdb seq, the pdb chain, the number of mismtaches and matches
def MapSeq2PDB(sequence, pdbfile, maxMisMatches=5, minMatchRatio=0.5):

	if not os.path.isfile(pdbfile):
		print 'ERROR: the pdb file does not exist: ', pdbfile
		exit(1)

	## extract PDB sequences by chains
	pdbseqs, residueLists, chains = ExtractSeqFromPDBFile(pdbfile)

	bestPDBSeq = None
	bestMapping = None
	bestResidueList = None
	bestChain = None
	minMisMatches = np.iinfo(np.int32).max
	maxMatches = np.iinfo(np.int32).min

	for pdbseq, residueList, chain in zip(pdbseqs, residueLists, chains):
		seq2pdb_mapping, numMisMatches, numMatches = MapSeq2ResidueList(sequence, pdbseq, residueList)
		if seq2pdb_mapping is None:
			continue
		if maxMatches < numMatches:
		#if numMisMatches < minMisMatches:
			bestMapping = seq2pdb_mapping
			minMisMatches = numMisMatches
			maxMatches = numMatches
			bestResidueList = residueList
			bestPDBSeq = pdbseq
			bestChain = chain

	if minMisMatches > maxMisMatches:
        	print 'ERROR: there are ', minMisMatches, ' mismatches between the query sequence and PDB file: ', pdbfile
		return None, None, None, None, None, None

	if maxMatches < min(30, minMatchRatio * len(sequence) ):
        	print 'ERROR: there are only ', maxMatches, ' matches on query sequence, less than ', minMatchRatio, ' of its length from PDB file: ', pdbfile
		return None, None, None, None, None, None

	return bestMapping, bestResidueList, bestPDBSeq, bestChain, minMisMatches, maxMatches


def ExtractCoordinatesByMapping(sequence, seq2pdb_mapping, residueList, atoms=['CB', 'CA', 'N', 'O', 'CG'], bUseAlternativeAtoms=True):
        neededAtoms = [ a.upper() for a in atoms ]

	numInvalidAtoms = dict()
        for atom in atoms:
                numInvalidAtoms[atom]=0

        atomCoordinates = []
        for i, j in zip(range(len(sequence)), seq2pdb_mapping):
                coordinates = dict()
                if j < 0:
                        for atom in neededAtoms:
                                coordinates[atom] = None
                        #coordinates['valid'] = False
                        atomCoordinates.append(coordinates)
                        continue

                res = residueList[j]
                AAname = three_to_one(res.get_resname() )
		if AAname == 'G':
			coordinates['vCB'] = BuildVirtualCB(res)

                for atom in neededAtoms:
                        coordinates[atom]= None
                        if atom == 'CG':
                                a = SelectCG(AAname, bUseAlternativeCG=bUseAlternativeAtoms)
                        elif atom == 'CB':
				a = SelectCB(AAname, bUseAlternativeCB=bUseAlternativeAtoms)
			else:
                                a = atom

			for atom2 in res:
				if atom2.get_id().upper() == a:
                               		coordinates[atom]= atom2.get_vector()
					break

			if coordinates[atom] is None:
				numInvalidAtoms[atom] += 1

                atomCoordinates.append(coordinates)

        return (atomCoordinates, numInvalidAtoms)

def ExtractCoordinatesBySeq(sequence, pdbfile, atoms=['CB', 'CA', 'N', 'O', 'CG'], maxMisMatches=5, minMatchRatio=0.5, bUseAlternativeAtoms=True):

        seq2pdb_mapping, residueList, pdbseq, chain, numMisMatches, numMatches = MapSeq2PDB(sequence, pdbfile, maxMisMatches=maxMisMatches, minMatchRatio=minMatchRatio)
	if seq2pdb_mapping is None:
		return None, None, None, None

	atomCoordinates, numInvalidAtoms = ExtractCoordinatesByMapping(sequence, seq2pdb_mapping, residueList, atoms=atoms, bUseAlternativeAtoms=bUseAlternativeAtoms)
	return (atomCoordinates, numInvalidAtoms), pdbseq, numMisMatches, numMatches

##calculate the distance matrix 
## coordinates is a list of None and Dict(). Each Dict has 3D coordinates for some atoms. The 3D coordinate of one atom is represented as Vector.  
##The results are saved in a dict(), in which each value is a matrix of float16
def CalcDistMatrix(coordinates, apts=['CbCb', 'CaCa', 'NO', 'CaCg', 'CgCg']):
        distMatrix = dict()

        for apt in apts:
                if len(apt) == 2:
                        atom1, atom2 = apt[0].upper(), apt[1].upper()
                elif len(apt) == 4:
                        atom1, atom2 = apt[:2].upper(), apt[2:].upper()
		else:
			print 'ERROR: unsupported atom pair types: ', apt
			exit(1)

                X = [ list(c[atom1]) if (c is not None and c[atom1] is not None) else [0,0,0] for c in coordinates ]
                Xvalid = [ 0 if (c is None or c[atom1] is None) else 1 for c in coordinates ]

                Y = [ list(c[atom2]) if (c is not None and c[atom2] is not None) else [0,0,0] for c in coordinates ]
                Yvalid = [ 0 if (c is None or c[atom2] is None) else 1 for c in coordinates ]

                dist = scipy.spatial.distance_matrix(X, Y).astype(np.float16)
                XYvalid = np.outer(Xvalid, Yvalid)
                np.putmask(dist, XYvalid==0, InvalidDistance)

		## set the self distance to 0
		if atom1 == atom2:
			np.fill_diagonal(dist, 0)

                distMatrix[apt] = dist

        return distMatrix

## this function assigns correct distance to two atoms of the same residue or two sequentially adajcent Ca atoms
def PostProcessDistMatrix(sequence, distMatrix, eps=0.0001):
	if distMatrix.has_key('CaCa'):
		## set the Ca-Ca distance of the same residue to 0
		np.fill_diagonal(distMatrix['CaCa'], eps)

		## set the Ca-Ca distance of two sequentially adajcent residues
		seqLen = distMatrix['CaCa'].shape[0]
		rng = np.arange(seqLen-1)
		for i in rng:
			if distMatrix['CaCa'][i, i+1] == InvalidDistance:
				distMatrix['CaCa'][i, i+1] = AdjacentCaCadist
				distMatrix['CaCa'][i+1, i] = AdjacentCaCadist

	## set the Cb-Cb distance of the same residue to 0
	if distMatrix.has_key('CbCb'):
		np.fill_diagonal(distMatrix['CbCb'], eps)

	## set the Cg-Cg distance of the same residue to 0
	if distMatrix.has_key('CgCg'):
		np.fill_diagonal(distMatrix['CgCg'], eps)

	## set the distance of N and O atoms of the same residue
	if distMatrix.has_key('NO'):
		seqLen = distMatrix['NO'].shape[0]
		assert seqLen == len(sequence)
		for AA, i in zip(sequence, range(seqLen)):
			if distMatrix['NO'][i, i] == InvalidDistance:
				distMatrix['NO'][i, i] = NOdists[AA.upper()]	

	## set the distance of Ca and Cg atoms of the same residue
	if distMatrix.has_key('CaCg'):
		seqLen = distMatrix['CaCg'].shape[0]
		assert seqLen == len(sequence)
		for AA, i in zip(sequence, range(seqLen)):
			if distMatrix['CaCg'][i, i] == InvalidDistance:
				distMatrix['CaCg'][i, i] = CaCgdists[AA.upper()]	

	return distMatrix

## calculate the inter-residue orientation using N, Ca, and Cb atoms
## coordinates shall contain a sequence of atom 3D coordinates
## the orientation angle is represented as Radian instead of Degree
def CalcTwoROriMatrix(coordinates):
        seqLen = len(coordinates)
        oriMatrix = dict()

	apts=['Ca1Cb1Cb2Ca2', 'N1Ca1Cb1Cb2', 'Ca1Cb1Cb2']
        for apt in apts:
                oriMatrix[apt] = np.full( (seqLen, seqLen), InvalidDegree, dtype=np.float16)

	numInvalidTwoROri = 0

	def ValidCB(c):
		if c['CB'] is None:
			return False

		## when CB is copied from CA, we do not use it to calculate some angles
		if c['CA'] is not None and np.linalg.norm(c['CB']-c['CA'])<0.1:
			return False
		return True

        for i in range(seqLen):
                ci = coordinates[i]
                if ci is None:
                        continue

		## we cannot replace ci CB by its CA atom since CA itself is needed for the three angles
		if ValidCB(ci):
			cicb = ci['CB']
		elif ci.has_key('vCB') and (ci['vCB'] is not None):
			cicb = ci['vCB']
		else:
			continue

                for j in range(seqLen):
                        if i == j:
                                continue

                        cj = coordinates[j]
                        if cj is None:
                                continue

                        if ci['CA'] is not None and cj['CA'] is not None:
				if ValidCB(cj):
                                        oriMatrix['Ca1Cb1Cb2Ca2'][i, j] = calc_dihedral(ci['CA'], cicb, cj['CB'], cj['CA']) * 180/np.pi
				elif cj.has_key('vCB') and cj['vCB'] is not None:
                                        oriMatrix['Ca1Cb1Cb2Ca2'][i, j] = calc_dihedral(ci['CA'], cicb, cj['vCB'], cj['CA']) * 180/np.pi
				## otherwise, assign InvalidDegree to Ca1Cb1Cb2Ca2 since we cannot replace cj CB by its CA atom

                        if ci['N'] is not None and ci['CA'] is not None:
				if ValidCB(cj):
                                        oriMatrix['N1Ca1Cb1Cb2'][i, j] = calc_dihedral(ci['N'], ci['CA'], cicb, cj['CB']) * 180/np.pi
				else:
					## replace cj CB by cj CA if the latter exists, otherwise check if vCB exists 
					if cj['CA'] is not None:
                                        	oriMatrix['N1Ca1Cb1Cb2'][i, j] = calc_dihedral(ci['N'], ci['CA'], cicb, cj['CA']) * 180/np.pi
					elif cj.has_key('vCB') and cj['vCB'] is not None:
                                        	oriMatrix['N1Ca1Cb1Cb2'][i, j] = calc_dihedral(ci['N'], ci['CA'], cicb, cj['vCB']) * 180/np.pi

                        if ci['CA'] is not None:
				if ValidCB(cj):
                                	oriMatrix['Ca1Cb1Cb2'][i, j] = calc_angle(ci['CA'], cicb, cj['CB']) * 180/np.pi
				else:
					## replace cj CB by cj CA if the latter exists, otherwise check if vCB exists
					if cj['CA'] is not None:
                                		oriMatrix['Ca1Cb1Cb2'][i, j] = calc_angle(ci['CA'], cicb, cj['CA']) * 180/np.pi
					elif cj.has_key('vCB') and cj['vCB'] is not None:
                                		oriMatrix['Ca1Cb1Cb2'][i, j] = calc_angle(ci['CA'], cicb, cj['vCB']) * 180/np.pi


	for apt in apts:
		np.fill_diagonal(oriMatrix[apt], ValueOfSelf)

	#oriMatrix['numInvalidTwoROri'] = numInvalidTwoROri

	return oriMatrix


##te the angle and dihedral of 4 Ca atoms (each is a vector), i1 and i2 represent two adjacent sequence positions and j1 and j2 represent another two adjacent positions
def CalcOrientationOf4CAs(CAi1, CAi2, CAj1, CAj2):
        angle = calc_angle(CAi1, CAi2, CAj1)*180/np.pi
	if CAj2 is None:
		angle2 = InvalidDegree
		dihedral = InvalidDegree
		dihedral2 = InvalidDegree
	else:
        	angle2 = calc_angle(CAi1, CAi2, CAj2)*180/np.pi
        	dihedral = calc_dihedral(CAi1, CAi2, CAj1, CAj2)*180/np.pi
        	dihedral2 = calc_dihedral(CAi1, CAi2, CAj2, CAj1)*180/np.pi

        return {'Ca1Ca2Ca3Ca4': dihedral, 'Ca1Ca2Ca3': angle, 'Ca1Ca2Ca4Ca3': dihedral2, 'Ca1Ca2Ca4': angle2}

## calculate the orientation defined by 4 Ca atoms. the first two Cas are adjacent and the 2nd two Cas are also adjacent.
## coordinates shall contain a sequence of atom 3D coordinates
def CalcCaOriMatrix(coordinates):
        seqLen = len(coordinates)
        oriMatrix = dict()

	apts = ['Ca1Ca2Ca3Ca4', 'Ca1Ca2Ca3', 'Ca1Ca2Ca4Ca3', 'Ca1Ca2Ca4']
        for apt in apts:
                oriMatrix[apt] = np.full( (seqLen, seqLen), InvalidDegree, dtype=np.float16)

	numInvalidCaOrientation = 0

        for i in range(seqLen-1):
                if coordinates[i] is None or coordinates[i+1] is None:
                        continue

                CAi = coordinates[i]['CA']
                CAi1 = coordinates[i+1]['CA']

		if CAi is None or CAi1 is None:
			continue

                for j in range(seqLen-1):
                        if i == j:
                                continue

                        if coordinates[j] is None or coordinates[j+1] is None:
                                continue

                        CAj = coordinates[j]['CA']
			CAj1 = coordinates[j+1]['CA']

			if CAj is None or CAj1 is None:
				continue

			orientation = CalcOrientationOf4CAs(CAi, CAi1, CAj, CAj1)
			for apt, v in orientation.iteritems():
				oriMatrix[apt][i,j] = v

		## deal with the case when j = seqLen-1
		if coordinates[seqLen-1] is not None:
			CAj = coordinates[seqLen-1]['CA']
			if CAj is not None:
				orientation = CalcOrientationOf4CAs(CAi, CAi1, CAj, None)
				for apt, v in orientation.iteritems():
					oriMatrix[apt][i, seqLen-1] = v

	for apt in apts:
		np.fill_diagonal(oriMatrix[apt], ValueOfSelf)

	#oriMatrix['numInvalidCaOrientation'] = numInvalidCaOrientation

        return oriMatrix


def ExtractDSSPByMapping(pdbfile, sequence=None, pdbseq=None, seq2pdb_mapping=None, residueList=None, chain=None):

	parser = StructureParser(pdbfile)
	structure = parser.get_structure('NoName', pdbfile)
	model = structure[0]

	try:
		dssp = DSSP(model, pdbfile, dssp='dsspcmbi')

	except PDBException as e:
		print 'EXCEPTION: ', str(e), ' from file: ', pdbfile
		exit(1)

	if sequence is None:
		return dssp

	SS8 = []
	ASA = []
	Phi = []
	Psi = []

	dsspFails = 0

	chainID = chain.get_id()
	for i, j in zip(range(len(sequence)), seq2pdb_mapping):
		if j < 0:
			SS8.append( None )
			ASA.append( None )
			Phi.append( None )
			Psi.append( None )
			continue
		
		resID = residueList[j].get_id()

		try:
			dsspInfo = dssp[ (chainID, resID) ]
			SS8.append(dsspInfo[2])
			ASA.append(dsspInfo[3])
			Phi.append(dsspInfo[4])
			Psi.append(dsspInfo[5])
		except:
			print 'EXCEPTION: there is no ', (chainID, resID), ' in the DSSP of ', pdbfile
			SS8.append( None )
                        ASA.append( None )
                        Phi.append( None )
                        Psi.append( None )
			dsspFails += 1

	result = dict()
	result['missing'] = [ 1 if s is None else 0 for s in SS8 ]
	result['SS8'] = [ s if (s is not None and s!='-' and s!='C' ) else 'L' for s in SS8 ]
	result['SS_str'] = ''.join( result['SS8'] )
	result['SS8_str'] = result['SS_str']

	## need to double check the conversion
	result['SS3'] = [ SS8Letter2SS3Letter[s] for s in result['SS8'] ]
	result['SS3_str'] = ''.join(result['SS3'])

	## to be compatible with old TPL
	result['pACC'] = [ acc*100 if acc is not None else 45 for acc in ASA ]
	result['ACC' ] = [ ]
	for pacc in result['pACC']:
		if pacc<=10:
			result['ACC'].append(0)
		elif pacc<=43:
			result['ACC'].append(1)
		else:
			result['ACC'].append(2)
	result['ACC'] = np.array(result['ACC'], dtype=np.int16)
	result['pACC'] = np.array(result['pACC'], dtype=np.int16)

	result['Phi'] = Phi
	result['Psi'] = Psi
	result['numDSSPFails'] = dsspFails

	## here we add gaps to be compatible with old TPL content
	result['DSSPsequence'] = ''.join([ pdbseq[i] if i>=0 else '-' for i in seq2pdb_mapping])

	return result

def ExtractDSSPBySeq(pdbfile, sequence=None):

	if sequence is None:
		return ExtractDSSPByMapping(pdbfile)
		
	seq2pdb_mapping, residueList, pdbseq, chain, numMisMatches, numMatches = MapSeq2PDB(sequence, pdbfile, maxMisMatches=maxMisMatches)
        if seq2pdb_mapping is None:
                return None, None, None, None

	result = ExtractDSSPByMapping(pdbfile, sequence, pdbseq, seq2pdb_mapping, residueList, chain)

	return result, pdbseq, numMisMatches, numMatches

def ExtractCoordinatesNDSSPBySeq(sequence, pdbfile, atoms=['CB', 'CA', 'N', 'O', 'CG'], maxMisMatches=5, minMatchRatio=0.5, bUseAlternativeAtoms=True):
	assert sequence is not None

	seq2pdb_mapping, residueList, pdbseq, chain, numMisMatches, numMatches = MapSeq2PDB(sequence, pdbfile, maxMisMatches=maxMisMatches, minMatchRatio=minMatchRatio)
        if seq2pdb_mapping is None:
                return None, None, None, None

	dssp = ExtractDSSPByMapping(pdbfile, sequence, pdbseq, seq2pdb_mapping, residueList, chain)

	##numInvalidAtoms is the number of atoms in PDB files without valid coordinates, but mapped by query sequence
	coordinates = ExtractCoordinatesByMapping(sequence, seq2pdb_mapping, residueList, atoms=atoms, bUseAlternativeAtoms=bUseAlternativeAtoms)

	return (coordinates, dssp), pdbseq, numMisMatches, numMatches
