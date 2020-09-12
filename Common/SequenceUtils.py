import os
import numpy as np
from Bio import SeqIO

## the 22 amino acids in the order of its 3-letter code
AA3LetterCode=['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
AA1LetterCode=['A',   'R',   'N',   'D',   'B',   'C',   'E',   'Q',   'Z',   'G',   'H',   'I',   'L',   'K',   'M',   'F',   'P',   'S',   'T',   'W',   'Y',   'V'  ]

## AAOrderBy3Letter is the alphabetical order of amino acids by its 3-letter code. This order is used in many amino acid substitution matrices.  map the two rare amino acids ASX and GLX to 20
AAOrderBy3Letter=[0,   1,     2,     3,    20,     4,     5,     6,     20,    7,     8,     9,     10,    11,    12,    13,    14,    15,    16,    17,    18,   19   ]

## AAOrderBy1Letter is the alphabetical order of amino acids by its 1-letter code. Again map the two rare amino acids B and Z to 20.
AAOrderBy1Letter=[0,   14,    11,    2,    20,     1,     3,     13,    20,    5,     6,     7,     9,      8,    10,    4,     12,    15,    16,    18,    19,   17   ]

##the 20 frequent amino acids
ValidAA1Letters=set(['A',   'R',    'N',  'D',    'C',  'E',   'Q',   'G',   'H',   'I',   'L',   'K',   'M',   'F',   'P',   'S',   'T',   'W',   'Y',   'V'])
Valid1AALetters = ValidAA1Letters
ValidAA3Letters=set(['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'])

## mapping between 1-letter code and 3-letter code
AA3LetterCode21LetterCode = {}
AA1LetterCode23LetterCode = {}

## map between 3-letter order and 1-letter order
AA3LetterOrder21LetterOrder = {}
AA1LetterOrder23LetterOrder = {}

##map from one AA (3-letter code or 1-letter code) to its order in terms of 3-letter code
AALetter2OrderOf3LetterCode = {}

## map from one AA (3-letter or 1-letter code) to its order in terms of 1-letter code
AALetter2OrderOf1LetterCode = {}

## initialize the mapping between letters and orders; all non-standard amino acids are mapped to 20
for i in range(26):
	letter = chr( ord('A') + i )
	AALetter2OrderOf3LetterCode[letter] = 20
	AALetter2OrderOf1LetterCode[letter] = 20


for l3, l1, o3, o1 in zip(AA3LetterCode, AA1LetterCode, AAOrderBy3Letter, AAOrderBy1Letter):
        AA1LetterCode23LetterCode[l1] = l3
        AA3LetterCode21LetterCode[l3] = l1

        AALetter2OrderOf3LetterCode[l1] = o3
        AALetter2OrderOf3LetterCode[l3] = o3

        AALetter2OrderOf1LetterCode[l1] = o1
        AALetter2OrderOf1LetterCode[l3] = o1

        AA3LetterOrder21LetterOrder[o3] = o1
        AA1LetterOrder23LetterOrder[o1] = o3

## default order of the amino acids
AAOrders = AALetter2OrderOf3LetterCode

## represent one AA by a binary vector of length 20. Only the 20 frequent amino acids are mapped. The other are represented as a zero vector.
AAVectors = np.zeros((26,20)).astype(np.int32)
for aa in ValidAA1Letters:
        index = ord(aa) - ord('A')
        AAVectors[index][ AAOrders[aa] ] = 1

## one-hot encoding of a protein sequence, represented as a L*20 binary matrix
def SeqOneHotEncoding(sequence):
        seq2int = (np.array(map(ord, sequence)) - ord('A') ).astype(np.int32)
        return AAVectors[seq2int]

## one-hot encoding of a protein sequence, represented as a L*21 binary matrix. gaps correspond to the last column.
def SeqOneHotEncodingWithGaps(sequence):
	m = np.zeros((len(sequence), 21), dtype=np.int32)
	seq = sequence.upper()
	for s, i in zip(seq, range(len(seq)) ):
		if s == '-':
			m[i][20] = 1
		else:
			index = ord(aa) - ord('A')
			m[i][:20] = AAVectors[index]
	return m

## convert a protein sequence into a list of orders (in terms of 3-letter code), which can be used to access of amino acid substitution matrices such as BLOSUMs
## these matrices are indexed by the alphabetical order of amino acids in their 3-letter code	
def Seq2OrderOf3LetterCode(sequence):
	indices = [ AALetter2OrderOf3LetterCode[aa] for aa in sequence.upper() ]		
	return indices

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

## read a FASTA file
def LoadFASTAFile(seqFile, returnStr=True):
	if not os.path.isfile(seqFile):
		print 'ERROR: an invalid sequence file: ', seqFile
		exit(1)
        record = SeqIO.read(seqFile, "fasta")
	if returnStr:
		return str(record.seq)

        return record.seq

## write to a FASTA file. sequence is a string of amino acid letters in upper case
def WriteFASTAFile(name, sequence, savefile):
	firstLine = '>' + name
	secondLine = sequence
	lines = '\n'.join([firstLine, secondLine])
	with open(savefile, 'w') as fh:
		fh.writelines(lines)
