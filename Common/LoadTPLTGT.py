import numpy as np
import os
import sys
import cPickle

import SequenceUtils
from SequenceUtils import ValidAA1Letters, Valid1AALetters, AALetter2OrderOf1LetterCode, AALetter2OrderOf3LetterCode, AA1LetterOrder23LetterOrder, DetectMultiHIS

from SSUtils import SS8Letter2Code, SS8Letter2SS3Code, SS8Letter2SS3Letter 
from LoadHHM import ReadHHM
from LoadHHM import load_hhm as LoadHHM

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
			print 'ERROR: wrong format in template Features section at line: ', line
			exit(1)

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
			if one_protein is None:
				print 'ERROR: failed to read the HHM section in ', tpl_file
				exit(1)
			continue

		i += 1
	
	## double check to see some required sections are read in
	for section in requiredSections:
		if not one_protein.has_key(section):
			print 'ERROR: one section for ', section, ' is missing in the tgt file: ', tgt_file
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
					print 'ERROR: wrong format at line: ', lines[ i + l ]
					exit(1)
				one_protein['DISO_prob'] = np.float32( fields[-1] )

			i += length
			continue

		if line.startswith('//////////// Original HHM'):
			i, one_protein = ReadHHM(lines, i, length, one_protein)
			if one_protein is None:
				print 'ERROR: failed to load the HHM section in ', tgt_file
				exit(1)
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

def load_hhm(hhmfile):
	return LoadHHM(hhmfile)

## mainly for test
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

	if protein is None:
		print 'ERROR: failed to load ', file
		exit(1)

	savefile = os.path.basename(file) + '.pkl'
	with open(savefile, 'wb') as fh:
		cPickle.dump(protein, fh, protocol=cPickle.HIGHEST_PROTOCOL)
