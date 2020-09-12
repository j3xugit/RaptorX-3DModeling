import os
import sys
import numpy as np
import cPickle

#from Bio.PDB.QCPSuperimposer import QCPSuperimposer
from Bio.SVDSuperimposer import SVDSuperimposer

from AlignmentUtils import ReadAlignment, ExtractSeqFromAlignment
from Common import PDBUtils

def Usage():
	print 'python EvaluateAlignmentQuality.py alnFile queryNativeFile templateNativeFile [savefolder]'
	print '	this script evaluates the quality of an alignment in terms of GDT and TMscore'
	print '	alnFile: the pairwise alignment file, in which template is usually placed before query'
	print '	queryNativeFile: the file containing the Ca coordinates of query, e.g., a PDB file'
	print ' templateNativeFile: the file containing the Ca coordinates of template, e.g., a PDB file'
	print '	the script will print out TM, GDT, uGDT, RMSD and also save all quality info to a cPickle file ending with .quality.pkl'

## superimpose two lists of atoms with coordinates. They shall have the same length
## atoms1 and atoms2 are an array of n*3 where n is the number of atoms
## in addition to RMSD, this function returns the per-position distance deviation
def Superimpose(atoms1, atoms2):
	assert len(atoms1) == len(atoms2)
	#aligner = QCPSuperimposer()
	aligner = SVDSuperimposer()
	aligner.set(atoms1, atoms2)
	aligner.run()
	RMSD = aligner.get_rms()

	## calculate the distance deviation at each position
	atoms2_transformed = aligner.get_transformed()
	diff = atoms1 - atoms2_transformed
	diff2 = np.power(diff, 2)
	deviations = np.sqrt(np.sum(diff2, axis=1))

	return RMSD, deviations	
	
## calc unnormalized GDT
def CalcUGDT(deviations):
	count, _ = np.histogram(deviations, bins=[0, 0.5, 1, 2, 4, 8])
	## get accumulative count
	count = np.cumsum(count)

	return sum(count[1:5])/4., sum(count[0:4])/4.

## calc TMscore
def CalcTM(deviations, normLen):
	if normLen<=15:
		d0=0.5
	else:
		d0=1.24*( (normLen-15.0)**(1./3.) )-1.8

	tm = [ 1/(1+d*d/d0/d0) for d in deviations ]
	tm = sum(tm)/normLen

	return tm
		
## this function returns the atoms in the aligned positions and also the residue indices of these atoms 
def MatchAtoms(alignment, tplAtoms, queryAtoms):
	tplSeqInAlignment, querySeqInAlignment = alignment

	## i is the residue index of template sequence, starting from 0
	## j is the residue index of query, starting from 0
	i = 0
	j = 0

	tplMatchedAtoms = []
	queryMatchedAtoms = []

	## queryMapping and tplMapping map the matched atoms back to residue index
	queryMapping = []
	tplMapping = []
	for t, q in zip(tplSeqInAlignment, querySeqInAlignment):
		if t!='-' and q!='-' and (tplAtoms[i] is not None) and (queryAtoms[j] is not None): 
			tplMatchedAtoms.append(tplAtoms[i])
			tplMapping.append(i)
			queryMatchedAtoms.append(queryAtoms[j])
			queryMapping.append(j)

		if t!='-':
			i += 1
		if q!='-':
			j += 1

	return np.array(tplMatchedAtoms), np.array(queryMatchedAtoms), tplMapping, queryMapping


## expand deviations of matched atoms to the whole sequence
## 999 is used for the deviation of unaligned positions
def ExpandDeviations(deviations, mapping, seqLen):
	fullDeviations = [ 999 ] * seqLen
	for d, index in zip(deviations, mapping):
		fullDeviations[index] = d
	return fullDeviations

if len(sys.argv)<4:
	Usage()
	exit(1)

alnFile = sys.argv[1]
queryNativeFile = sys.argv[2]
templateNativeFile = sys.argv[3]

savefolder=os.getcwd()
if len(sys.argv) >=5:
	savefolder = sys.argv[4]

alignment = ReadAlignment(alnFile)
tplSeq, querySeq = ExtractSeqFromAlignment(alignment)

## get the template Ca coordinates by template sequence
tplInfo, _, tplNumMisMatches, _ = PDBUtils.ExtractCoordinatesBySeq(tplSeq, templateNativeFile, atoms=['CA'] )
if tplNumMisMatches > 5:
	print 'ERROR: too many residue mismatches between template sequence and its PDB file'
	exit(1)
tplCaAtoms, tplNumInvalidAtoms = tplInfo

## get the query Ca coordinates by query sequence
queryInfo, _, queryNumMisMatches, _ = PDBUtils.ExtractCoordinatesBySeq(querySeq, queryNativeFile, atoms=['CA'] )
if queryNumMisMatches > 5:
	print 'ERROR: too many residue mismatches between query sequence and its PDB file'
	exit(1)
queryCaAtoms, queryNumInvalidAtoms = queryInfo

## data type conversion
tplAtoms = [ list(a['CA']) if a['CA'] is not None else None for a in tplCaAtoms ]
queryAtoms = [ list(a['CA']) if a['CA'] is not None else None for a in queryCaAtoms ]

assert len(tplAtoms) == len(tplSeq)
assert len(queryAtoms) == len(querySeq)

## superimpose tplCaAtoms and queryCaAtoms based upon the alignment
tplMatchedAtoms, queryMatchedAtoms, tplMapping, queryMapping = MatchAtoms(alignment, tplAtoms, queryAtoms)
if len(tplMatchedAtoms) < 5:
	print 'ERROR: there are fewer than 5 aligned residues with valid coordinates in the alignment: ', alnFile
	exit(1)

RMSD, deviations = Superimpose(tplMatchedAtoms, queryMatchedAtoms)

uGDT, uGHA = CalcUGDT(deviations)
GDT = uGDT * 100./len(querySeq)
GHA = uGHA * 100./len(querySeq)
TM = CalcTM(deviations, len(querySeq) )

print 'Quality of', os.path.basename(alnFile), ':', TM, GDT, uGDT, RMSD

fullDeviations = np.array( ExpandDeviations(deviations, queryMapping, len(querySeq)) ).astype(np.float16)
#print fullDeviations

result = dict()
result['alnfile'] = alnFile
result['TM'] = TM
result['GDT'] = GDT
result['uGDT'] = uGDT
result['GHA'] = GHA
result['uGHA'] = uGHA
result['RMSD'] = RMSD
result['deviations'] = fullDeviations

savefile = os.path.basename(alnFile).split('.')[0] + '.quality.pkl'
savefile = os.path.join(savefolder, savefile)
with open(savefile, 'w') as fh:
	cPickle.dump(result, fh, protocol=cPickle.HIGHEST_PROTOCOL)
