import numpy as np
import os
import sys
import cPickle

import SequenceUtils
from SequenceUtils import ValidAA1Letters, Valid1AALetters, AALetter2OrderOf1LetterCode, AALetter2OrderOf3LetterCode, AA1LetterOrder23LetterOrder

from SSUtils import SS8Letter2Code, SS8Letter2SS3Code, SS8Letter2SS3Letter 

"""
This script reads in a tpl or tgt file. The tpl and tgt files contain both obsolete and useful data, so we have to be very careful in using them.
Each protein is stored as a python dict().
To use the position-specfic frequency matrix, please use the keyword PSFM.
To use the position-specific scoring matrix, please use the keyword PSSM.
Both PSFM and PSSM actually encode information derived from the profile HMM built by HHpred or HHblits,
so there is no need to directly use the keys containing 'hmm'

Other keys such as psp and psm are obosolete. Please DO NOT use them.
Especially, DO NOT use a mix of psp(psm) and PSFM(PSSM) in calculating the similarity of two aligned positions since their columns are arranged in a different order.
PSFM and PSSM columns are arranged by the alphabetical order of amino acids in their 1-letter code
psp and psm columns are arranged by the alphabetical order of amino acids in their 3-letter code

The predicted secondary structure information (key SS2) in tpl file shall not be used.
"""

## the rows and cols of gonnet matrix are arranged in the alpahbetical oder of amino acids in the 3-letter code
## the below matrix is a mutation probability matrix?
gonnet = [ 
[ 1.7378,  0.870964,0.933254,0.933254, 1.12202,  0.954993, 1,        1.12202,  0.831764, 0.831764,  0.758578, 0.912011, 0.851138, 0.588844, 1.07152,  1.28825,  1.14815,   0.436516,  0.60256,  1.02329],
[ 0.870964,2.95121, 1.07152, 0.933254, 0.60256,  1.41254,  1.09648,  0.794328, 1.14815,  0.57544,   0.60256,  1.86209,  0.676083, 0.47863,  0.812831, 0.954993, 0.954993,  0.691831,  0.660693, 0.630957],
[ 0.933254,1.07152, 2.39883, 1.65959,  0.660693, 1.1749,   1.23027,  1.09648,  1.31826,  0.524807,  0.501187, 1.20226,  0.60256,  0.489779, 0.812831, 1.23027,  1.12202,   0.436516,  0.724436, 0.60256],
[ 0.933254,0.933254,1.65959, 2.95121,  0.47863,  1.23027,  1.86209,  1.02329,  1.09648,  0.416869,  0.398107, 1.12202,  0.501187, 0.354813, 0.851138, 1.12202,  1,         0.301995,  0.524807, 0.512861],
[ 1.12202, 0.60256, 0.660693,0.47863, 14.1254,   0.57544,  0.501187, 0.630957, 0.74131,  0.776247,  0.707946, 0.524807, 0.812831, 0.831764, 0.489779, 1.02329,  0.891251,  0.794328,  0.891251, 1], 
[ 0.954993,1.41254, 1.1749,  1.23027,  0.57544,  1.86209,  1.47911,  0.794328, 1.31826,  0.645654,  0.691831, 1.41254,  0.794328, 0.549541, 0.954993, 1.04713,  1,         0.537032,  0.676083, 0.707946],
[ 1,       1.09648, 1.23027, 1.86209,  0.501187, 1.47911,  2.29087,  0.831764, 1.09648,  0.537032,  0.524807, 1.31826,  0.630957, 0.40738,  0.891251, 1.04713,  0.977237,  0.371535,  0.537032, 0.645654],
[ 1.12202, 0.794328,1.09648, 1.02329,  0.630957, 0.794328, 0.831764, 4.57088,  0.724436, 0.354813,  0.363078, 0.776247, 0.446684, 0.301995, 0.691831, 1.09648,  0.776247,  0.398107,  0.398107, 0.467735],
[ 0.831764,1.14815, 1.31826, 1.09648,  0.74131,  1.31826,  1.09648,  0.724436, 3.98107,  0.60256,   0.645654, 1.14815,  0.74131,  0.977237, 0.776247, 0.954993, 0.933254,  0.831764,  1.65959,  0.630957],
[ 0.831764,0.57544, 0.524807,0.416869, 0.776247, 0.645654, 0.537032, 0.354813, 0.60256,  2.51189,   1.90546,  0.616595, 1.77828,  1.25893,  0.549541, 0.660693, 0.870964,  0.660693,  0.851138, 2.04174],
[ 0.758578,0.60256, 0.501187,0.398107, 0.707946, 0.691831, 0.524807, 0.363078, 0.645654, 1.90546,   2.51189,  0.616595, 1.90546,  1.58489,  0.588844, 0.616595, 0.74131,   0.851138,  1,        1.51356],
[ 0.912011,1.86209, 1.20226, 1.12202,  0.524807, 1.41254,  1.31826,  0.776247, 1.14815,  0.616595,  0.616595, 2.0893,   0.724436, 0.467735, 0.870964, 1.02329,  1.02329,   0.446684,  0.616595, 0.676083],
[ 0.851138,0.676083,0.60256, 0.501187, 0.812831, 0.794328, 0.630957, 0.446684, 0.74131,  1.77828,   1.90546,  0.724436, 2.69153,  1.44544,  0.57544,  0.724436, 0.870964,  0.794328,  0.954993, 1.44544],
[ 0.588844,0.47863, 0.489779,0.354813, 0.831764, 0.549541, 0.40738,  0.301995, 0.977237, 1.25893,   1.58489,  0.467735, 1.44544,  5.01187,  0.416869, 0.524807, 0.60256,   2.29087,   3.23594,  1.02329],
[ 1.07152, 0.812831,0.812831,0.851138, 0.489779, 0.954993, 0.891251, 0.691831, 0.776247, 0.549541,  0.588844, 0.870964, 0.57544,  0.416869, 5.7544,   1.09648,  1.02329,   0.316228,  0.489779, 0.660693],
[ 1.28825, 0.954993,1.23027, 1.12202,  1.02329,  1.04713,  1.04713,  1.09648,  0.954993, 0.660693,  0.616595, 1.02329,  0.724436, 0.524807, 1.09648,  1.65959,  1.41254,   0.467735,  0.645654, 0.794328],
[ 1.14815, 0.954993,1.12202, 1,        0.891251, 1,        0.977237, 0.776247, 0.933254, 0.870964,  0.74131,  1.02329,  0.870964, 0.60256,  1.02329,  1.41254,  1.77828,   0.446684,  0.645654, 1], 
[ 0.436516,0.691831,0.436516,0.301995, 0.794328, 0.537032, 0.371535, 0.398107, 0.831764, 0.660693,  0.851138, 0.446684, 0.794328, 2.29087,  0.316228, 0.467735, 0.446684, 26.3027,    2.5704,   0.549541],
[ 0.60256, 0.660693,0.724436,0.524807, 0.891251, 0.676083, 0.537032, 0.398107, 1.65959,  0.851138,  1,        0.616595, 0.954993, 3.23594,  0.489779, 0.645654, 0.645654,  2.5704,    6.0256,   0.776247],
[ 1.02329, 0.630957,0.60256, 0.512861, 1,        0.707946, 0.645654, 0.467735, 0.630957, 2.04174,   1.51356,  0.676083, 1.44544,  1.02329,  0.660693, 0.794328, 1,         0.549541,  0.776247, 2.18776] 
]

gonnet = np.array(gonnet, np.float32)

M_M, M_I, M_D, I_M, I_I, D_M, D_D, _NEFF, I_NEFF, D_NEFF = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

## HMMNull is the background score of amino acids, in the alphabetical order of 1-letter code
HMMNull = [3706,5728,4211,4064,4839,3729,4763,4308,4069,3323,5509,4640,4464,4937,4285,4423,3815,3783,6325,4665,0]
HMMNull = np.array(HMMNull, dtype=np.float32)


def DetectMultiHIS(sequence):
	HISflag = np.zeros( (len(sequence)), np.uint8)
	##scan from left to right
	pos1 = sequence.find('HHH')
	if pos1>=0 and pos1<10 : 
		i = pos1 + 3
		while i<len(sequence) and sequence[i] == 'H':
			i = i + 1
		HISflag[:i] = 1

	pos2 = sequence.rfind('HHH')
	if pos2!=-1 and pos2>len(sequence)-10 :
		i = pos2 -1
	        while i>=0 and sequence[i] == 'H':
			i = i - 1
		HISflag[i+1:] = 1

	return HISflag

def ReadTemplateFeatures(lines, start_position, length, one_protein):
	i = start_position
	one_protein['templateFeatureHeader'] = lines[i:i+2]
	i += 2

	## missing indicates which residues have no coordinates
	one_protein['missing'] = np.zeros( (length), np.uint8)

	one_protein['SS3Coding'] = np.zeros( (length, 3), np.float32)
	one_protein['SS8Coding'] = np.zeros( (length, 8), np.float32)
	one_protein['SS_str'] = ''
	#   one_protein['core'] = np.zeros(length, np.int )
	one_protein['ACC'] = np.zeros(length, np.uint8)
	one_protein['pACC'] = np.zeros(length, np.uint16)
	one_protein['CNa'] = np.zeros(length, np.uint8)
	one_protein['CNb'] = np.zeros(length, np.uint8)
	one_protein['Ca'] = np.ones( (length, 3), np.float32) * (-999)
	one_protein['Cb'] = np.ones( (length, 3), np.float32) * (-999)

	seqStr = ''
	for l in range(length):
		fields = lines[i+l].split()
		if len(fields) < 9:
			print 'Error: wrong format in template Features section at line: ', line
			exit(-1)

		seqStr += fields[1]
		one_protein['missing'][l] = np.uint( fields[2])

		one_protein['SS3Coding'][l, SS8Letter2SS3Code[ fields[3]] ] = 1
		one_protein['SS8Coding'][l, SS8Letter2Code[ fields[3]] ] = 1
		one_protein['SS_str'] += ( 'L' if fields[3] == 'C' else fields[3] )
		one_protein['ACC'][l] = np.uint8(fields[5])

		#one_protein['core'][i] = np.int( fields[4] )

		## pACC is the percentage of ACC
		one_protein['pACC'][l] = np.uint16(fields[6])

       		##CNa and CNb are the number of contacts formed by the Ca and Cb atoms, respectively
		one_protein['CNa'][l] = np.uint8(fields[7])
		one_protein['CNb'][l] = np.uint8(fields[8])


		if len(fields) == 9:
			continue
		elif len(fields) >= 15:
       			## coordinates of Ca and Cb atoms
			one_protein['Ca'][l] = [np.float32(num) for num in fields[9:12]]
			one_protein['Cb'][l] = [np.float32(num) for num in fields[12:15]]
		else:
			print 'ERROR: wrong format at line: ', lines[ i+ l ]
			exit(1)

	assert (length == len(seqStr) )
	assert all( seqLetter == aa for seqLetter, aa in zip(seqStr, one_protein['sequence']) )

	i += length
	return i, one_protein

def ReadPSP(lines, start_position, length, one_protein):
	##From left to right, the columns of PSM are arranged in the alphabetical order of amino acids in 3-letter code

	i = start_position
	one_protein['PSPHeader'] = lines[i: i+2 ]
			
	## psp1 equivalent to psm
	one_protein['psp1'] = np.zeros( (length, 20), np.int8)

	## psp2 is the position-specific frequency matrix, represented by an integer between 0 and 100
	one_protein['psp2'] = np.zeros( (length, 20), np.uint8)

	## psp3 is the rightmost 2 columns in the PSP file
	one_protein['psp3'] = np.zeros( (length, 2), np.float32)

	i += 2
	seqStr = ''
	for l in range(length):
		fields = lines[ i+l ].split()
		seqStr += fields[1]

		cont = fields[2:]
		assert (len(cont) >=42 )

		one_protein['psp1'][l] = [np.int8(num) for num in cont[0:20] ]
		one_protein['psp2'][l] = [np.uint8(num) for num in cont[20:40] ]
		one_protein['psp3'][l] = [np.float32(num) for num in cont[40:]]


	if ( len(seqStr) != len(one_protein['sequence']) ):
		print 'ERROR: PSP sequence length inconsistent with the original sequence length for protein: ', one_protein['name']
		print 'PSP primary sequence: ', seqStr
		print 'original sequence: ', one_protein['sequence']
		exit(1)

	i += length

	return i, one_protein

def ReadPSM(lines, start_position, length, one_protein):
	##From left to right, the columns of PSM are arranged in the alphabetical order of amino acids in 3-letter code

	i = start_position
	one_protein['PSMHeader'] = lines[i:i+2]
	i += 2

	rows = []
	for l in range(length):
		one_row = [np.int16(num) for num in lines[ i + l ].split() ]
		assert ( len(one_row) == 20 )
		rows.append(one_row)

	one_protein['psm'] = np.array(rows) 
	i += length
	return i, one_protein

## a function reading HMM block for tgt and tpl files. We also use this function to read the HMM block from a profleHMM file (generated by HHpred/HHblits package)
## for the tgt and tpl file, the number of lines for the header is 6; for the profileHMM file, the header has 4 lines
def ReadHHM(lines, start_position, length, one_protein, numLines4Header=6):

	i = start_position
	one_protein['HMMHeader'] = lines[i: i+ numLines4Header]
	i += numLines4Header

	## the columns of hmm1 are in the alphabetical order of amino acids in 1-letter code, different from the above PSP and PSM matrices
	one_protein['hmm1'] = np.zeros( (length, 20), np.float32)
	one_protein['hmm2'] = np.zeros( (length, 10), np.float32)
	one_protein['hmm1_prob'] = np.zeros((length, 20), np.float32)
	one_protein['hmm1_score'] = np.zeros((length, 20), np.float32)

	seqStr = ''
	
	for l in range(length):
		##this line is for emission score. The amino acids are ordered from left to right alphabetically by their 1-letter code
		fields = lines[ i + l*3 + 0].replace("*", "99999").split()

		assert len(fields) == 23
		one_protein['hmm1'][l] = np.array([ - np.int32(num) for num in fields[2: -1] ])/1000.
		aa = fields[0]
		seqStr += aa
		"""
		if  (l + 1) != np.int32(fields[1]):
			print 'Error: inconsistent residue number in file for protein ', one_protein['name'], ' at line: ', lines[ i + l*3 + 0 ]
			exit(-1)
		"""

		##the first 7 columns of this line is for state transition
		one_protein['hmm2'][l][0:7] = [ np.exp(-np.int32(num)/1000.0*0.6931) for num in lines[i + l*3 + 1].replace("*", "99999").split()[0:7]]

		##the last 3 columns of this line is for Neff of Match, Insertion and Deletion
		one_protein['hmm2'][l][7:10] = [ np.int32(num)/1000.0 for num in lines[i + l*3 + 1].split()[7:10]]

		## _NEFF is for match, I_NEFF for insertion and D_NEFF for deletion. More comments are needed for the below code
		rm = 0.1
		one_protein['hmm2'][l][M_M] = (one_protein['hmm2'][l][_NEFF]*one_protein['hmm2'][l][M_M] + rm*0.6)/(rm + one_protein['hmm2'][l][_NEFF])
		one_protein['hmm2'][l][M_I] = (one_protein['hmm2'][l][_NEFF]*one_protein['hmm2'][l][M_I] + rm*0.2)/(rm + one_protein['hmm2'][l][_NEFF])
		one_protein['hmm2'][l][M_D] = (one_protein['hmm2'][l][_NEFF]*one_protein['hmm2'][l][M_D] + rm*0.2)/(rm + one_protein['hmm2'][l][_NEFF])

		ri = 0.1
		one_protein['hmm2'][l][I_I] = (one_protein['hmm2'][l][I_NEFF]*one_protein['hmm2'][l][I_I] + ri*0.75)/(ri + one_protein['hmm2'][l][I_NEFF])
		one_protein['hmm2'][l][I_M] = (one_protein['hmm2'][l][I_NEFF]*one_protein['hmm2'][l][I_M] + ri*0.25)/(ri + one_protein['hmm2'][l][I_NEFF])

		rd = 0.1
		one_protein['hmm2'][l][D_D] = (one_protein['hmm2'][l][D_NEFF]*one_protein['hmm2'][l][D_D] + rd*0.75)/(rd + one_protein['hmm2'][l][D_NEFF])
		one_protein['hmm2'][l][D_M] = (one_protein['hmm2'][l][D_NEFF]*one_protein['hmm2'][l][D_M] + rd*0.25)/(rd + one_protein['hmm2'][l][D_NEFF])
		

		one_protein['hmm1_prob'][l,] = pow(2.0, one_protein['hmm1'][l,])
		wssum = sum(one_protein['hmm1_prob'][l, ])

		#print 'l = ', l, 'sum= ', wssum

		## renormalize to make wssum = 1
		if wssum > 0 : 
			one_protein['hmm1_prob'][l, ] /= wssum
		else:
			one_protein['hmm1_prob'][l, AALetter2OrderOf1LetterCode[aa] ] = 1.

		"""
		## if the probability sum is not equal to 1
		if abs(wssum - 1.) > 0.1 :
			one_protein['hmm1_prob'][l, ] = 0
			one_protein['hmm1_prob'][l, AALetter2OrderOf1LetterCode[aa] ] = 1.
		"""


		## add pseudo count		
		g = np.zeros( (20), np.float32)
		for j in range(20):
			orderIn3LetterCode_j = AA1LetterOrder23LetterOrder[j]
			for k in range(20):
				orderIn3LetterCode_k = AA1LetterOrder23LetterOrder[k]
				g[j] += one_protein['hmm1_prob'][l, k] * gonnet[ orderIn3LetterCode_k, orderIn3LetterCode_j ]
			g[j] *= pow(2.0, -1.0*HMMNull[j] / 1000.0)

		#print 'l=', l, ' gsum= ', sum(g)
		## sum(g) is very close to 1, here we renormalize g to make its sum to be exactly 1
		g = g/sum(g)

		ws_tmp_neff = one_protein['hmm2'][l][_NEFF] - 1
		one_protein['hmm1'][l, ] = (ws_tmp_neff * one_protein['hmm1_prob'][l, ] + g*10) / (ws_tmp_neff+10)
		
		## recalculate the emission score and probability after pseudo count is added	
		one_protein['hmm1_prob'][l,] = one_protein['hmm1'][l, ]
		one_protein['hmm1'][l, ] = np.log2(one_protein['hmm1_prob'][l, ])
		one_protein['hmm1_score' ][l, ] = one_protein['hmm1'][l, ] + HMMNull[:20]/1000.0

		## PSFM: position-specific frequency matrix, PSSM: position-specific scoring matrix
		one_protein['PSFM'] = one_protein['hmm1_prob']
		one_protein['PSSM'] = one_protein['hmm1_score']

	#assert ( seqStr == one_protein['sequence'] )
	if len(seqStr) != len(one_protein['sequence']):
		print 'ERROR: inconsistent sequence length in HMM section and orignal sequence for protein: ', one_protein['name'] 
		exit(1)

	comparison = [ (aa=='X' or bb=='X' or aa==bb) for aa, bb in zip(seqStr, one_protein['sequence']) ]
	if not all(comparison):
		print 'ERROR: inconsistent sequence between HMM section and orignal sequence for protein: ', one_protein['name']
		print ' original seq: ', one_protein['sequence']
		print ' HMM seq: ', seqStr
		exit(1)

	return i + 3*length, one_protein

## this function reads a profile HMM file generated by HHpred/HHblits package
def load_hhm(hhmfile):
	with open(hhmfile, 'r') as fh:
		content = [ r.strip() for r in list(fh) ]
	if not bool(content):
		print 'ERROR: empty profileHMM file: ', hhmfile
		exit(1)
	if not content[0].startswith('HHsearch'):
		print 'ERROR: this file may not be a profileHMM file generated by HHpred/HHblits: ', hhmfile
		exit(1)
	if len(content) < 10:
		print 'ERROR: this profileHMM file is too short: ', hhmfile
		exit(1)

	requiredSections = ['name', 'length', 'sequence', 'NEFF',  'hmm1', 'hmm2', 'hmm1_prob', 'hmm1_score', 'PSFM', 'PSSM', 'DateCreated']	
	protein = {}

	## get sequence name
	if not content[1].startswith('NAME '):
		print 'ERROR: the protein name shall appear at the second line of profileHMM file: ', hhmfile
		exit(1)
	fields = content[1].split()
	if len(fields) < 2:
		print 'ERROR: incorrect name format in profileHMM file: ', hhmfile
		exit(1)
	protein['name'] = fields[1]

	i = 0
	while i < len(content):
		row = content[i]
		if len(row)<1:
			i += 1
			continue

		if row.startswith('DATE '):
			protein['DateCreated'] = row[6:]
			i += 1
			continue

		if row.startswith('NEFF '):
			protein['NEFF'] = np.float32(row.split()[1])
			i += 1
			continue

		if row.startswith('LENG '):
			protein['length'] = np.int32(row.split()[1])
			i += 1
			continue

		if row.startswith('>ss_pred'):
			## read the predicted secondary structure
			start = i+1
			end = i+1
			while not content[end].startswith('>'):
				end += 1
			protein['SSEseq'] = ''.join(content[start:end]).replace('C', 'L')
			if len(protein['SSEseq']) != protein['length']:
				print 'ERROR: inconsistent sequence length and predicted SS sequence length in hmmfile: ', hhmfile
				exit(1)
			i = end
			continue

		if row.startswith('>ss_conf'):
			## read the predicted secondary structure confidence score
			start = i+1
			end = i+1
			while not content[end].startswith('>'):
				end += 1

			SSEconfStr = ''.join(content[start:end])
			protein['SSEconf'] = [ np.int16(score) for score in SSEconfStr ]

			if len(protein['SSEconf']) != protein['length']:
				print 'ERROR: inconsistent sequence length and predicted SS confidence sequence length in hhmfile: ', hhmfile
				exit(1)

			i = end
			continue
			

		if row.startswith('>' + protein['name']):
			## read in the sequence in the following lines
			start = i+1
			end = i+1
			while not content[end].startswith('>') and ( not content[end].startswith('#') ):
				end += 1

			## at this point, content[end] shall start with >
			protein['sequence'] = ''.join(content[start:end])
			if len(protein['sequence']) != protein['length']:
				print 'ERROR: inconsistent sequence length in hmmfile: ', hhmfile
				exit(1)
			i = end
			continue

		if len(row) == 1 and row[0]=='#' and content[i+1].startswith('NULL') and content[i+2].startswith('HMM'):
			i, protein = ReadHHM(content, i+1, protein['length'], protein, numLines4Header=4)
			continue

		i += 1

	
	## double check to see some required sections are read in
	for section in requiredSections:
		if not protein.has_key(section):
			print 'ERROR: one section for ', section, ' is missing in the hmm file: ', hhmfile
			print 'ERROR: it is also possible that the hmm file has a format incompatible with this script.'
			exit(1)

	protein['requiredSections'] = requiredSections

	return protein
	
	
def load_tpl(tpl_file):

	if tpl_file.endswith('.tpl.pkl'):
		with open(tpl_file, 'rb') as fh:
			return cPickle.load(fh)

	### For position-sepcific frequency matrix, please use one_protein['PSFM']. For position-specific scoring matrix, please use one_protein['PSSM']
	### it is better not to use psm, psp1, psp2, psp3, hmm1, hmm2, hmm1_prob
	### Be careful if you want to use SS8Coding and SS3Coding. Please make sure that the order of 8-state and 3-state SS types is consistent with other source code.
	### SS_str is the string of 8-state secondary structure types

	requiredSections = ['name', 'length', 'chain_id', 'HISflag', 'sequence', 'DSSPsequence', 'SS_str', 'ACC', 'pACC', 'NEFF', 'DateCreated', 'psm', 'psp1', 'psp2', 'psp3', 'hmm1', 'hmm2', 'hmm1_prob', 'PSFM', 'PSSM', 'SS3Coding', 'SS8Coding', 'missing', 'CNa', 'CNb', 'Ca', 'Cb']
	one_protein = {}

	f = open(tpl_file, 'r')
	lines = [ l.strip() for l in list(f) ]
	f.close()

	length = None

	i = 0
	while ( i<len(lines) ):
		line = lines[i]
		
		if line.startswith('Version '):
			one_protein['Version'] = line.split()[-1]
			i += 1
			continue

		if line.startswith('Template Name'):
			one_protein['name'] = line.split('=')[1].strip()
			i += 1
			continue

		if line.startswith('Chain ID'):
			fields = line.split('=')
			if len(fields)<2:
				one_protein['chain_id'] = ''
			else:
				one_protein['chain_id'] = fields[1].strip()
			i += 1
			continue

		if line.startswith('Length '):
			one_protein['length'] = np.int32( line.split('=')[1] )
			length = one_protein['length']
			assert (length >= 1)
			##print 'length=', length

			i += 1
			continue
	
		if line.startswith('SEQRES sequence'):	
			one_protein['sequence'] = line.split('=')[1].strip()
			one_protein['HISflag'] = DetectMultiHIS(one_protein['sequence'])
			##print 'seq length=', len(one_protein['sequence'])

			i += 1
			continue

		if line.startswith('DSSP '):
			##DSSP sequence, which usually is a subsequence of 'sequence'. The residues with missing coordinates are represented as a hypen
			one_protein['DSSPsequence'] = line.split('=')[1].strip()

			i += 1
			continue

		if line.startswith('NEFF '):
			one_protein['NEFF'] = np.float32( line.split('=')[1] )
			assert (one_protein['NEFF'] >= 1.)

			i += 1
			continue

		if line.startswith('Date '):
			one_protein['DateCreated'] = line.split('=')[1]
			i += 1
			continue

		if line.startswith('//////////// Features'):
			i, one_protein = ReadTemplateFeatures(lines, i, length, one_protein)
			continue
		
		if line.startswith('//////////// Original PSM'):
			i, one_protein = ReadPSM(lines, i, length, one_protein)
			continue

		if line.startswith('//////////// Original PSP'):
			i, one_protein = ReadPSP(lines, i, length, one_protein)
			continue

		if line.startswith('//////////// Original SS2'):
			"""
			predicted secondary structure probability by PSIPRED or DeepCNF. We simply skip this section
			since 1) we have native SS for templates 2) the format of this section is inconsistent with the other SS prediction files

			one_protein['SS2'] = np.zeros( (length, 3), np.float32)
			for l in range(length):
				one_protein['SS2'][l] = [ np.float32(num) for num in lines[ i + l ].split()[3:] ]
			"""
			one_protein['obsoleteSS2'] = lines[i: i+2+length]
			i += (2 + length)
			continue

		if line.startswith('//////////// Original HHM'):
			i, one_protein =ReadHHM(lines, i, length, one_protein)
			continue

		i += 1
	
	## double check to see some required sections are read in
	for section in requiredSections:
		if not one_protein.has_key(section):
			print 'Error: one section for ', section, ' is missing in the tgt file: ', tgt_file
			print '	it is also possible that the tgt file has a wrong or different format compatible with this script.'
			exit(1)

	one_protein['requiredSections'] = requiredSections	

	assert ( length == len(one_protein['sequence']) )
	assert ( length == len(one_protein['DSSPsequence']) )
	assert all( letter=='X' or (letter in Valid1AALetters) for letter in one_protein['sequence'] )
	assert all( l1 == l2 for l1, l2 in zip(one_protein['sequence'], one_protein['DSSPsequence']) if l2 != '-' )

	return one_protein

def load_tgt(tgt_file):

	if tgt_file.endswith('.tgt.pkl'):
		with open(tgt_file, 'rb') as fh:
			return cPickle.load(fh)

	### For position-sepcific frequency matrix, please use one_protein['PSFM']. For position-specific scoring matrix, please use one_protein['PSSM']
	### It is better not use psm, psp1, psp2, psp3, hmm1, hmm2, hmm1_prob

	requiredSections = ['name', 'length', 'sequence', 'SSEseq', 'SSEconf', 'ACCseq', 'ACCconf', 'NEFF', 'EVD', 'DateCreated', 'psm', 'psp1', 'psp2', 'psp3', 'hmm1', 'hmm2', 'hmm1_prob', 'PSFM', 'PSSM', 'SS3', 'SS8', 'ACC', 'ACC_prob' ]
	one_protein = {}

	f = open(tgt_file)
	lines = [ l.strip() for l in list(f) ]
	f.close()
	length = 0

	i = 0
	while (i < len(lines) ):
		line = lines[i]
		if line.startswith('Sequence Name'):
			one_protein['name'] = line.split() [-1]
			i += 1
			continue
		if line.startswith('Length '):
			one_protein['length'] = np.int32( line.split()[-1] )
			length = one_protein['length']
			i += 1
			continue
		if line.startswith('Sequence ') and not line.startswith('Sequence Name'):
			one_protein['sequence'] = line.split()[-1]
			i += 1
			continue

		if line.startswith('SSEseq '):
			one_protein['SSEseq'] = line.split()[-1].replace('C', 'L')
			i += 1
			continue

		if line.startswith('SSEconf '):
			SSEconfStr = line.split()[-1]
			one_protein['SSEconf'] = [ np.int16(score) for score in SSEconfStr ]
			i += 1
			continue
		if line.startswith('ACCseq '):
			one_protein['ACCseq'] = line.split()[-1]
			i += 1
			continue
		if line.startswith('ACCconf '):
			one_protein['ACCconf'] = line.split()[-1]
			i += 1
			continue
		if line.startswith('NEFF '):
			one_protein['NEFF'] = np.float32( line.split()[-1] )
			i += 1
			continue

		if line.startswith('EVD '):
			fields = line.split()
			if len(fields) <= 2:
				one_protein['EVD'] = np.array([0, 1]).astype(np.float32)
			elif len(fields) >=4:
				one_protein['EVD'] = np.array(map(np.float32, fields[2:]) ).astype(np.float32)
			else:
				print 'ERROR: wrong format in line: ', line
				exit(1)
			i += 1
			continue

		if line.startswith('Date '):
			one_protein['DateCreated'] = ' '.join(line.split()[2:] )
			i += 1
			continue

		if line.startswith('//////////// Original PSM'):
			i, one_protein = ReadPSM(lines, i, length, one_protein)
			continue

		if line.startswith('//////////// Original PSP'):
			i, one_protein = ReadPSP(lines, i, length, one_protein)
			continue

		if line.startswith('//////////// Original DIS'):
			one_protein['DISOHeader'] = lines[ i: i+2]
			i += 2
			## read disorder prediction
			for l in range(length):
				fields = lines[ i + l ]. split()
				if len(fields) != 5:
					print 'Error: wrong format at line: ', lines[ i + l ]
					exit(-1)
				one_protein['DISO_prob'] = np.float32( fields[-1] )

			i += length
			continue

		if line.startswith('//////////// Original HHM'):

			i, one_protein = ReadHHM(lines, i, length, one_protein)
			continue

		if line.startswith('//////////// Original SS3+SS8+ACC'):
			one_protein['SSACCheader'] = lines[i:i+3]	
			assert ( one_protein['SSACCheader'][2] == 'SS3:   H     E     L  | SS8:   H     G     I     E     B     T     S     L  | ACC:  Bury   Medium  Exposed' )

			i += 3

			one_protein['SS3'] = np.zeros( (length, 3), np.float32 )
			one_protein['SS8'] = np.zeros( (length, 8), np.float32 )
			one_protein['ACC_prob'] = np.zeros( (length, 3), np.float32)

			for l in range(length):
				cont = lines[ i + l ].split()[:-2]
				one_protein['SS3'][l] = [np.float32(num) for num in cont[0:3]]
				one_protein['SS8'][l] = [np.float32(num) for num in cont[3:11]]
				one_protein['ACC_prob'][l] = [np.float32(num) for num in cont[11:]]

			one_protein['ACC'] = np.argmax( one_protein['ACC_prob'], axis=1 ).astype(np.uint8)

			i += length
			continue

		i += 1

	## double check to see some required sections are read in
	for section in requiredSections:
		if not one_protein.has_key(section):
			print 'ERROR: one section for ', section, ' is missing in the tgt file: ', tgt_file
			print '	it is also possible that the tgt file has a wrong or different format compatible with this script.'
			exit(1)

	one_protein['requiredSections'] = requiredSections

	return one_protein

## for test
if __name__ == "__main__":
	if len(sys.argv) < 2:
		print 'python LoadTPLTGT.py tpl_file or tgt_file or hhm_file'
		print '	the input file shall end with .tgt, .tpl. or .hhm'
		exit(1)

	file = sys.argv[1]

	if file.endswith('.tgt'):
		protein = load_tgt(file)
	elif file.endswith('.tpl'):
		protein = load_tpl(file)
	elif file.endswith('.hhm'):
		protein = load_hmm(file)
	else:
		print 'ERROR: the input file shall have suffix .tgt or .tpl or .hhm'
		exit(1)

	savefile = os.path.basename(file) + '.pkl'
	fh = open(savefile, 'wb')
	cPickle.dump( protein, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	fh.close()
