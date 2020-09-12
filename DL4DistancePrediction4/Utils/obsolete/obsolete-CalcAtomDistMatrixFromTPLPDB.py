import numpy as np
import sys
import os
import cPickle

sys.path.append('../')
import config

from Bio.PDB import *
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import three_to_one

##this script calculate dist matrix for a protein. Two arguments are given: tpl file and PDB file.
## we calculate the distance matrix based upon SEQRES record in the tpl file by mapping it to DSSP sequence
## the PDB file contains coordinates for residues in DSSP sequence

def Usage():
	print 'python CalcAtomDistMatrixFromPDB.py tpl_file  pdb_file'

if len(sys.argv)<3:
	Usage()
	exit(-1)	

tplfile = sys.argv[1]

f = open(tplfile, 'r')
tplContent = [ line.strip() for line in list(f) ]
f.close()

if not tplContent[4].startswith('SEQRES'):
        print 'incorrect tpl file format at line: ', tplContent[4]
        exit(-1)

if not tplContent[5].startswith('DSSP'):
        print 'incorrect tpl file format at line: ', tplContent[5]
        exit(-1)
SEQRESsequence = tplContent[4].split('= ')[-1]
DSSPsequence = tplContent[5].split('= ')[-1]

name = os.path.basename(tplfile).split('.')[0]
if len(name)!=5:
	print 'the protein target name is incorrect. It must be composed of PDB ID and chain letter'
	exit(-1)

pdbfile = sys.argv[2]
parser = PDBParser()
structure = parser.get_structure(name, pdbfile)

chain=name[4]
residues = structure[0][chain].get_residues()
residueList = [ r for r in residues if is_aa(r) ]
#numResidues = len(residueList)
pdbseq = ''.join( [ three_to_one(r.get_resname()) for r in residueList ] )

#print pdbseq

### check if DSSPsequence is equivalent to pdbseq
validDSSPseq = DSSPsequence.replace('-', '')

if validDSSPseq != pdbseq:
	print 'Inconsistency between DSSPsequence in ', tplfile, ' and pdbseq in ', pdbfile
	print 'pdbseq: ', pdbseq
	print 'DSPseq: ', validDSSPseq
	diffs = [i for i in xrange(min(len(pdbseq), len(validDSSPseq) ) ) if pdbseq[i] != validDSSPseq[i] ]
	print 'inconsistent positions: ', diffs

	exit(-1)

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
		print 'abnormal in the alignment between SEQRESsequence and DSSPsequence.'
		exit(-1)

if numMisMatches > 0:
	print 'Some non-identical residues found between SEQRESsequence and DSSPsequence in ', tplfile
	exit(-1)


"""
from Bio import SeqIO
record = SeqIO.read(seqfile, "fasta")
##print(record.seq)

##align two sequences

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum80
alignments = pairwise2.align.localds(pdbseq, record.seq, blosum80, -5, -1)
##alignments = pairwise2.align.localxx(pdbseq, record.seq)

print '#alignments:', len(alignments)

##find the alignment with the minimum residue number difference
diffs = []
for alignment in alignments:
	print alignment[0]
	print alignment[1]

	mapping_pdb2seq = np.ones( (len(pdbseq)), np.int32) * (-1)

	i = -1
	j = -1

	for a, b in zip(alignment[0], alignment[1]):
		if a!='-' and b!='-':
			i += 1
			j += 1
			mapping_pdb2seq[i] = j

		elif a!='-':
			i += 1
		elif b!='-':
			j += 1
		else:
			print 'abnormal in the alignment'

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

	diffs.append(diff)

diffs = np.array(diffs)

print '#alignments:', len(alignments)
print 'best alignment:', diffs.argmin()

alignment = alignments[diffs.argmin()]
		
##print alignments
print alignment[0]
print alignment[1]

pdbseq_aligned = alignment[0]
fastaseq_aligned = alignment[1]


i=-1
j=-1
numMisMatches = 0

mapping_seq2pdb = np.ones( (len(record.seq)), np.int32) * (-1)

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
		

assert j == len(record.seq) - 1

print 'numMisMatches: ', numMisMatches

if numMisMatches > 5:
	print 'ERROR: there are more than 5 mismatches among the aligned positions for the seq file ', seqfile
	exit(-1)
"""

##calculate the distance matrix
fullSeqLen = len(SEQRESsequence)

CaCaMatrix = np.ones( (fullSeqLen, fullSeqLen), np.float16) * (-1)
CbCbMatrix = np.ones( (fullSeqLen, fullSeqLen), np.float16) * (-1)
CgCgMatrix = np.ones( (fullSeqLen, fullSeqLen), np.float16) * (-1)
CaCgMatrix = np.ones( (fullSeqLen, fullSeqLen), np.float16) * (-1)
NOMatrix   = np.ones( (fullSeqLen, fullSeqLen), np.float16) * (-1)


for i, j in zip(range(fullSeqLen), mapping_seq2pdb):
	if j < 0:
		continue

	#print 'j=', j

	res1 = residueList[j]
	#print res1

	ca1 = res1['CA']
	n1 = res1['N']
	o1 = res1['O']


	if res1.get_resname().upper() == 'GLY':
		cb1 = res1['CA']
	else:
		cb1 = res1['CB']

	cg1 = (config.SelectCG( three_to_one( res1.get_resname() ) ) ).upper()

	#print 'cg1=', cg1

	if res1.has_id(cg1):
		cg1 = res1[ cg1 ]
	else:
		cg1 = None

	"""
	if res1.has_id('CG'):
		cg1 = res1['CG']
	elif res1.has_id('CG1'):
		cg1 = res1['CG1']
	elif res1.has_id('CG2'):
		cg1 = res1['CG2']
	else:
		cg1 = cb1
	"""


	"""
	print record.seq[i], res1.get_resname()
	print n1
	print o1
	print ca1
	print cb1
	print cg1
	"""

	for k, l in zip(range(fullSeqLen), mapping_seq2pdb):
		if l <0 :
			continue

		res2 = residueList[l]

		#print 'l=', l
		#print res2

		ca2 = res2['CA']
		n2 = res2['N']
		o2 = res2['O']

		if res2.get_resname().upper() == 'GLY':
			cb2 = res2['CA']
		else:
			cb2 = res2['CB']

		cg2 = (config.SelectCG ( three_to_one( res2.get_resname() )  ) ).upper()

		if res2.has_id(cg2):
			cg2 = res2[ cg2 ]
		else:
			cg2 = None

		"""
		if res2.has_id('CG'):
			cg2 = res2['CG']
		elif res2.has_id('CG1'):
			cg2 = res2['CG1']
		elif res2.has_id('CG2'):
			cg2 = res2['CG2']
		else:
			cg2 = cb2
		"""

		## Ca-Ca distance matrix, symmetric
		CaCaMatrix[i, k] = ca1 - ca2

		## Cb-Cb, symmetric
		CbCbMatrix[i, k] = cb1 - cb2

		## Cg-Cg, symmetric
		if (cg1 is not None) and (cg2 is not None):
			CgCgMatrix[i, k] = cg1 - cg2

		## Ca-Cg, asymmetric
		if cg2 is not None:
			CaCgMatrix[i, k] = ca1 - cg2

		## N-O, asymmetric
		NOMatrix[i, k]  = n1 - o2

		"""
		if i==320 and k==44:
			print j, l, cb1.get_vector(), cb2.get_vector(), CbCbMatrix[i, k]
		"""

np.fill_diagonal(CaCaMatrix, 0)
np.fill_diagonal(CbCbMatrix, 0)
np.fill_diagonal(CgCgMatrix, 0)

atomDistMatrix = dict()
atomDistMatrix['CaCa'] = CaCaMatrix
atomDistMatrix['CbCb'] = CbCbMatrix
atomDistMatrix['CgCg'] = CgCgMatrix
atomDistMatrix['CaCg'] = CaCgMatrix
atomDistMatrix['NO']   = NOMatrix
atomDistMatrix['name'] = name
atomDistMatrix['pdbseq'] = pdbseq
atomDistMatrix['seq4matrix'] = SEQRESsequence

savefile = name + '.atomDistMatrix.pkl'
fh = open(savefile, 'wb')
cPickle.dump(atomDistMatrix, fh, protocol=cPickle.HIGHEST_PROTOCOL)
fh.close()

