import os
import sys
import numpy as np
import cPickle

from SequenceUtils import AALetter2OrderOf1LetterCode, AA1LetterOrder23LetterOrder, DetectMultiHIS
"""
This script reads an .hhm file and save it as a python dict().
To use the position-specfic frequency matrix, please use the keyword PSFM.
To use the position-specific scoring matrix, please use the keyword PSSM.
PSFM and PSSM are derived from the HMM block, so there is no need to directly use the keys containing 'hmm'.
PSFM and PSSM columns are arranged by the alphabetical order of amino acids in their 1-letter code.

the following keys are also available: name, sequence, NEFF, DateCreated
"""

M_M, M_I, M_D, I_M, I_I, D_M, D_D, _NEFF, I_NEFF, D_NEFF = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

## HMMNull is the background score of amino acids, in the alphabetical order of 1-letter code
HMMNull = np.array([3706,5728,4211,4064,4839,3729,4763,4308,4069,3323,5509,4640,4464,4937,4285,4423,3815,3783,6325,4665,0], dtype=np.float32)

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

## a function reading the HMM block in .hhm,. tgt and .tpl files. The HMM block is generated by HHpred/HHblits package
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
		return i+3*length, None
		#exit(1)

	comparison = [ (aa=='X' or bb=='X' or aa==bb) for aa, bb in zip(seqStr, one_protein['sequence']) ]
	if not all(comparison):
		print 'ERROR: inconsistent sequence between HMM section and orignal sequence for protein: ', one_protein['name']
		print ' original seq: ', one_protein['sequence']
		print ' HMM seq: ', seqStr
		return i+3*length, None
		#exit(1)

	return i + 3*length, one_protein

## this function reads a profile HMM file generated by HHpred/HHblits package
def load_hhm(hhmfile):
	if hhmfile.endswith('.hhm.pkl'):
		with open(hhmfile, 'rb') as fh:
			return cPickle.load(fh)

	with open(hhmfile, 'r') as fh:
		content = [ r.strip() for r in list(fh) ]
	fh.close()

	if not bool(content):
		print 'ERROR: empty profileHMM file: ', hhmfile
		exit(1)
	if not content[0].startswith('HHsearch '):
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
	protein['annotation'] = ' '.join(fields[2:])

	i = 2
	while i < len(content):
		row = content[i]
		if len(row)<1:
			i += 1
			continue

		if row.startswith('DATE '):
			protein['DateCreated'] = row[6:]
			i += 1
			continue

		if row.startswith('LENG '):
			protein['length'] = np.int32(row.split()[1])
			i += 1
			continue

		if row.startswith('FILT '):
			protein['filterinfo'] = row[len('FILT '):]
			i += 1
			continue

		if row.startswith('NEFF '):
			protein['NEFF'] = np.float32(row.split()[1])
			i += 1
			continue

		if row.startswith('>ss_dssp '):
			## read native secondary structure
                        start = i+1
                        end = i+1
                        while not content[end].startswith('>'):
                                end += 1
			protein['nativeSS8'] = ''.join(content[start:end]).replace('C', 'L').replace('-', 'L')
			if len(protein['nativeSS8']) != protein['length']:
				print 'ERROR: inconsistent sequence length and native SS length in hmmfile: ', hhmfile
				exit(1)
			i = end
			continue

		if row.startswith('>ss_pred'):
			## read predicted secondary structure
			start = i+1
			end = i+1
			while not content[end].startswith('>'):
				end += 1
			protein['SSEseq'] = ''.join(content[start:end]).replace('C', 'L')
			if len(protein['SSEseq']) != protein['length']:
				print 'ERROR: inconsistent sequence length and predicted SS length in hmmfile: ', hhmfile
				exit(1)
			i = end
			continue

		if row.startswith('>ss_conf'):
			## read predicted secondary structure confidence score
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
			
		if row.startswith('>' + protein['name']) and (content[i-1].strip() != 'SEQ') and (not protein.has_key('sequence')):
			## read the primary sequence in the following lines
			start = i+1
			end = i+1
			while not content[end].startswith('>') and ( not content[end].startswith('#') ):
				end += 1

			## at this point, content[end] shall start with >
			protein['sequence'] = ''.join(content[start:end])
			if len(protein['sequence']) != protein['length']:
				print 'ERROR: inconsistent sequence length in hmmfile: ', hhmfile
				return None
				#exit(1)
			i = end
			protein['HISflag'] = DetectMultiHIS(protein['sequence'])
			continue

		if len(row) == 1 and row[0]=='#' and content[i+1].startswith('NULL') and content[i+2].startswith('HMM'):
			nullfields = content[i+1].split()[1:]
			HMMNull[:20] = [ np.float32(f) for f in nullfields ]

			i, protein = ReadHHM(content, i+1, protein['length'], protein, numLines4Header=4)
			if protein is None:
				return None
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

## for test only
if __name__ == "__main__":
        if len(sys.argv) < 2:
                print 'python LoadHHM.py hhm_file'
                print ' the input file shall end with .hhm'
                exit(1)

        file = sys.argv[1]

        if file.endswith('.hhm') or file.endswith('.hhm.pkl'):
                protein = load_hhm(file)
		if protein is None:
			print 'ERROR: failed to load the hhm file: ', file
			exit(1)
        else:
                print 'ERROR: the input file shall have suffix .hhm'
                exit(1)

        savefile = os.path.basename(file) + '.pkl'
        with open(savefile, 'wb') as fh:
        	cPickle.dump(protein, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	
