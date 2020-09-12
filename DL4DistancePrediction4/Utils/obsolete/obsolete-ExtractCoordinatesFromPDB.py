import numpy as np
import sys
import os
import cPickle

from Bio.PDB import *
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import three_to_one

sys.path.append(os.path.join(os.environ['ModelingHome'], 'Common'))
from LoadTPLTGT import load_tpl as LoadTPL

## This script extract 3D coordinates for Cb, Ca, Cg, N amd O atoms for a template from a PDB file. 
## Two arguments are given: tpl file and PDB file.
## The tpl file may be inconsistent with the PDB file, but we will extract coordinates from PDB file based upon the tpl file

def Usage():
	print 'python ExtractCoordinatesFromPDB.py tpl_file  pdb_file'

if len(sys.argv)<3:
	Usage()
	exit(-1)	

tplfile = sys.argv[1]
pdbfile = sys.argv[2]
name = os.path.basename(tplfile).split('.')[0]

parser = PDBParser()
structure = parser.get_structure(name, pdbfile)

pdbseq = ''

"""
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
	#print(pp.get_sequence())
	pdbseq += pp.get_sequence()
"""

if len(name)!=5:
	print 'the template name is incorrect. It must be composed of PDB ID and chain letter'
	exit(-1)

chain=name[4]
residues = structure[0][chain].get_residues()
residueList = [ r for r in residues if is_aa(r) ]
#numResidues = len(residueList)
pdbseq = ''.join( [ three_to_one(r.get_resname()) for r in residueList ] )

#print pdbseq

#from Bio import SeqIO
tpl = LoadTPL(tplfile)
tplseq = tpl['sequence']

#record = SeqIO.read(seqfile, "fasta")
##print(record.seq)

##align two sequences

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum80
##alignments = pairwise2.align.localds(pdbseq, tplseq, blosum80, -5, -1)
alignments = pairwise2.align.localms(pdbseq, tplseq, 3, -1, -0.5, -0.0 )

#print '#alignments:', len(alignments)

##find the alignment with the minimum residue number difference
diffs = []
for alignment in alignments:
	print alignment[0]
	print alignment[1]

	mapping_pdb2seq = np.ones( (len(pdbseq)), np.int32) * (-1)

	i = -1
	j = -1

	CaDiff = 0

	for a, b in zip(alignment[0], alignment[1]):
		if a!='-' and b!='-':
			i += 1
			j += 1
			mapping_pdb2seq[i] = j

			tplCa = Vector(tpl['Ca'][j] )
			pdbCa = residueList[i]['CA'].get_vector()
			CaDiff += np.linalg.norm( tplCa - pdbCa )

		elif a!='-':
			i += 1
		elif b!='-':
			j += 1
		else:
			print 'abnormal in the alignment'

	"""
	##calculate the residue number difference

	diff = 0
	for i in range(1, len(pdbseq)):
		if mapping_pdb2seq[i] == -1:
			diff += 10
			continue
		prev_mapping = mapping_pdb2seq[i-1]
		if prev_mapping == -1:
			continue
		res_id = residueList[i].get_id()[1]
		prev_res_id = residueList[i-1].get_id()[1]
		id_diff = max(1, res_id - prev_res_id)
		mapping_diff = mapping_pdb2seq[i] - prev_mapping
		diff += abs(id_diff - mapping_diff)

	print 'residue number diff: ', diff
	diffs.append(diff)
	"""

	print 'Ca diff: ', CaDiff
	diffs.append(CaDiff)

diffs = np.array(diffs)

#print '#alignments:', len(alignments)
#print 'best alignment:', diffs.argmin()

alignment = alignments[diffs.argmin()]
		
##print alignments
print alignment[0]
print alignment[1]

pdbseq_aligned = alignment[0]
fastaseq_aligned = alignment[1]

i=-1
j=-1
numMisMatches = 0

mapping_seq2pdb = np.ones( (len(tplseq)), np.int32) * (-1)

for a, b in zip(pdbseq_aligned, fastaseq_aligned):
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
		print 'abnormal in the alignment'
		

assert j == len(tplseq) - 1

print 'numMisMatches: ', numMisMatches

if numMisMatches > 5:
	print 'ERROR: there are more than 5 mismatches among the aligned positions for the seq file ', seqfile
	exit(-1)

## extract the coordinates for Ca, Cb, Cg, N and O

"""
CaCaMatrix = np.ones( (len(record.seq), len(record.seq)), np.float16) * (-1)
CbCbMatrix = np.ones( (len(record.seq), len(record.seq)), np.float16) * (-1)
CgCgMatrix = np.ones( (len(record.seq), len(record.seq)), np.float16) * (-1)
CaCgMatrix = np.ones( (len(record.seq), len(record.seq)), np.float16) * (-1)
NOMatrix   = np.ones( (len(record.seq), len(record.seq)), np.float16) * (-1)
"""

atomCoordinates = []


for i, j in zip(range(len(tplseq)), mapping_seq2pdb):
	if j < 0:

		assert (tpl['missing'][i]==1)

		coordinates = dict()
		fake = Vector(x=-999., y=-999., z=-999.)
		coordinates['Ca'] = fake
		coordinates['Cg'] = fake
		coordinates['Cb'] = fake
		coordinates['N'] = fake
		coordinates['O'] = fake
		coordinates['valid'] = False

		atomCoordinates.append(coordinates)
		continue


	#print 'j=', j

	#assert tpl['missing'][i] == 0
	if tpl['missing'][i] != 0:
		print 'Error: not a missing residue at index ', i
		exit(-1)

	res1 = residueList[j]

	ca1 = res1['CA']
	n1 = res1['N']
	o1 = res1['O']

	if res1.has_id('CB'):
		cb1 = res1['CB']
	else:
		cb1 = res1['CA']

	if res1.has_id('CG'):
		cg1 = res1['CG']
	elif res1.has_id('CG1'):
		cg1 = res1['CG1']
	elif res1.has_id('CG2'):
		cg1 = res1['CG2']
	else:
		cg1 = cb1


	"""
	print record.seq[i], res1.get_resname()
	print n1
	print o1
	print ca1
	print cb1
	print cg1
	"""
	coordinates = dict()
	coordinates['Ca'] = ca1.get_vector()
	coordinates['Cg'] = cg1.get_vector()
	coordinates['Cb'] = cb1.get_vector()
	coordinates['N'] = n1.get_vector()
	coordinates['O'] = o1.get_vector()
	coordinates['valid'] = True

	assert np.linalg.norm(Vector(tpl['Ca'][i]) - coordinates['Ca']) < 0.1

	atomCoordinates.append(coordinates)


savefile = name + '.atomCoordinates.pkl'
fh = open(savefile, 'wb')
cPickle.dump(atomCoordinates, fh, protocol=cPickle.HIGHEST_PROTOCOL)
fh.close()

"""
np.savetxt(name + '.CaCg.txt', CaCgMatrix, fmt='%.2f')
np.savetxt(name + '.CaCa.txt', CaCaMatrix, fmt='%.2f')
np.savetxt(name + '.CbCb.txt', CbCbMatrix, fmt='%.2f')
"""
