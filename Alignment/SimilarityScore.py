import sys
import os
import copy
import numpy as np

from Common import SequenceUtils

"""
In this script, the rows and cols of all the mutation matrices are arranged in the alphabetical order of amino acids in their 3-letter code.
That is, A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, Z
"""

Ori_BLOSUM_62 = [
[  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -5 ],  #A
[ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -5 ],  #R
[ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3, -5 ],  #N
[ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3, -5 ],  #D
[  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -5 ],  #C
[ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2, -5 ],  #Q
[ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2, -5 ],  #E
[  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -5 ],  #G
[ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3, -5 ],  #H
[ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -5 ],  #I
[ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -5 ],  #L
[ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2, -5 ],  #K
[ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -5 ],  #M
[ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -5 ],  #F
[ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -5 ],  #P
[  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2, -5 ],  #S
[  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -5 ],  #T
[ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -5 ],  #W
[ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -5 ],  #Y
[  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -5 ],  #V
[ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 ]]  #Z
#  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z

Ori_BLOSUM_62 = np.array(Ori_BLOSUM_62, dtype=np.int32)
newBLOSUM62 = copy.deepcopy(Ori_BLOSUM_62)
newBLOSUM62[20,:] = 0
newBLOSUM62[:,20] = 0


Ori_BLOSUM_45 = [
[  5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -2, -2,  0, -5 ],  #A  
[ -2,  7,  0, -1, -3,  1,  0, -2,  0, -3, -2,  3, -1, -2, -2, -1, -1, -2, -1, -2, -5 ],  #R  
[ -1,  0,  6,  2, -2,  0,  0,  0,  1, -2, -3,  0, -2, -2, -2,  1,  0, -4, -2, -3, -5 ],  #N  
[ -2, -1,  2,  7, -3,  0,  2, -1,  0, -4, -3,  0, -3, -4, -1,  0, -1, -4, -2, -3, -5 ],  #D  
[ -1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -5 ],  #C  
[ -1,  1,  0,  0, -3,  6,  2, -2,  1, -2, -2,  1,  0, -4, -1,  0, -1, -2, -1, -3, -5 ],  #Q  
[ -1,  0,  0,  2, -3,  2,  6, -2,  0, -3, -2,  1, -2, -3,  0,  0, -1, -3, -2, -3, -5 ],  #E  
[  0, -2,  0, -1, -3, -2, -2,  7, -2, -4, -3, -2, -2, -3, -2,  0, -2, -2, -3, -3, -5 ],  #G  
[ -2,  0,  1,  0, -3,  1,  0, -2, 10, -3, -2, -1,  0, -2, -2, -1, -2, -3,  2, -3, -5 ],  #H  
[ -1, -3, -2, -4, -3, -2, -3, -4, -3,  5,  2, -3,  2,  0, -2, -2, -1, -2,  0,  3, -5 ],  #I  
[ -1, -2, -3, -3, -2, -2, -2, -3, -2,  2,  5, -3,  2,  1, -3, -3, -1, -2,  0,  1, -5 ],  #L  
[ -1,  3,  0,  0, -3,  1,  1, -2, -1, -3, -3,  5, -1, -3, -1, -1, -1, -2, -1, -2, -5 ],  #K  
[ -1, -1, -2, -3, -2,  0, -2, -2,  0,  2,  2, -1,  6,  0, -2, -2, -1, -2,  0,  1, -5 ],  #M  
[ -2, -2, -2, -4, -2, -4, -3, -3, -2,  0,  1, -3,  0,  8, -3, -2, -1,  1,  3,  0, -5 ],  #F  
[ -1, -2, -2, -1, -4, -1,  0, -2, -2, -2, -3, -1, -2, -3,  9, -1, -1, -3, -3, -3, -5 ],  #P  
[  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -3, -1, -2, -2, -1,  4,  2, -4, -2, -1, -5 ],  #S  
[  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1,  2,  5, -3, -1,  0, -5 ],  #T  
[ -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2,  1, -3, -4, -3, 15,  3, -3, -5 ],  #W  
[ -2, -1, -2, -2, -3, -1, -2, -3,  2,  0,  0, -1,  0,  3, -3, -2, -1,  3,  8, -1, -5 ],  #Y  
[  0, -2, -3, -3, -1, -3, -3, -3, -3,  3,  1, -2,  1,  0, -3, -1,  0, -3, -1,  5, -5 ],  #V  
[ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 ]] #Z
# A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z   

Ori_BLOSUM_45 = np.array(Ori_BLOSUM_45, dtype=np.int32)
newBLOSUM45 = copy.deepcopy(Ori_BLOSUM_45)
newBLOSUM45[20,:] = 0
newBLOSUM45[:,20] = 0

Ori_BLOSUM_80 = [
[  5, -2, -2, -2, -1, -1, -1,  0, -2, -2, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -5 ],  #A
[ -2,  6, -1, -2, -4,  1, -1, -3,  0, -3, -3,  2, -2, -4, -2, -1, -1, -4, -3, -3, -5 ],  #R
[ -2, -1,  6,  1, -3,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -4, -3, -4, -5 ],  #N
[ -2, -2,  1,  6, -4, -1,  1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4, -5 ],  #D
[ -1, -4, -3, -4,  9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -5 ],  #C
[ -1,  1,  0, -1, -4,  6,  2, -2,  1, -3, -3,  1,  0, -4, -2,  0, -1, -3, -2, -3, -5 ],  #Q
[ -1, -1, -1,  1, -5,  2,  6, -3,  0, -4, -4,  1, -2, -4, -2,  0, -1, -4, -3, -3, -5 ],  #E
[  0, -3, -1, -2, -4, -2, -3,  6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -5 ],  #G
[ -2,  0,  0, -2, -4,  1,  0, -3,  8, -4, -3, -1, -2, -2, -3, -1, -2, -3,  2, -4, -5 ],  #H
[ -2, -3, -4, -4, -2, -3, -4, -5, -4,  5,  1, -3,  1, -1, -4, -3, -1, -3, -2,  3, -5 ],  #I
[ -2, -3, -4, -5, -2, -3, -4, -4, -3,  1,  4, -3,  2,  0, -3, -3, -2, -2, -2,  1, -5 ],  #L
[ -1,  2,  0, -1, -4,  1,  1, -2, -1, -3, -3,  5, -2, -4, -1, -1, -1, -4, -3, -3, -5 ],  #K
[ -1, -2, -3, -4, -2,  0, -2, -4, -2,  1,  2, -2,  6,  0, -3, -2, -1, -2, -2,  1, -5 ],  #M
[ -3, -4, -4, -4, -3, -4, -4, -4, -2, -1,  0, -4,  0,  6, -4, -3, -2,  0,  3, -1, -5 ],  #F
[ -1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4,  8, -1, -2, -5, -4, -3, -5 ],  #P
[  1, -1,  0, -1, -2,  0,  0, -1, -1, -3, -3, -1, -2, -3, -1,  5,  1, -4, -2, -2, -5 ],  #S
[  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2,  1,  5, -4, -2,  0, -5 ],  #T
[ -3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2,  0, -5, -4, -4, 11,  2, -3, -5 ],  #W
[ -2, -3, -3, -4, -3, -2, -3, -4,  2, -2, -2, -3, -2,  3, -4, -2, -2,  2,  7, -2, -5 ],  #Y
[  0, -3, -4, -4, -1, -3, -3, -4, -4,  3,  1, -3,  1, -1, -3, -2,  0, -3, -2,  4, -5 ],  #V
[ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 ]] #Z
# A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z 

Ori_BLOSUM_80 = np.array(Ori_BLOSUM_80, dtype=np.int32)

newBLOSUM80 = copy.deepcopy(Ori_BLOSUM_80)
newBLOSUM80[20,:] = 0
newBLOSUM80[:,20] = 0

## both a and b are 1-letter code
def BLOSUM62(a, b):

        if ord(a) < ord('A') or ord(a) > ord('Z'): 
                print "a non standard amino acid: " + a
                return 0 # a non-standard amino acid, neither reward nor penalize
        
	orderIn3LetterCode_a = SequenceUtils.AALetter2OrderOf3LetterCode[a]

        if ord(b) < ord('A') or ord(b) > ord('Z'):
                print "a non standard amino acid: " + b
                return 0 #a non-standard amino acid
        
	orderIn3LetterCode_b = SequenceUtils.AALetter2OrderOf3LetterCode[b]

	return Ori_BLOSUM_62[orderIn3LetterCode_a, orderIn3LetterCode_b]

def BLOSUM45(a, b):
        
        if ord(a) < ord('A') or ord(a) > ord('Z'):
                print "a non standard amino acid: " + a
                return 0 # a non-standard amino acid, neither reward nor penalize
        
	orderIn3LetterCode_a = SequenceUtils.AALetter2OrderOf3LetterCode[a]
        
	if ord(b) < ord('A') or ord(b) > ord('Z'):
                print "a non standard amino acid: " + b
                return 0 # //a non-standard amino acid
        
	orderIn3LetterCode_b = SequenceUtils.AALetter2OrderOf3LetterCode[b]

	return Ori_BLOSUM_45[orderIn3LetterCode_a, orderIn3LetterCode_b]

def BLOSUM80(a, b):

        if ord(a) < ord('A') or ord(a) > ord('Z'):
                print "a non standard amino acid: " + a
                return 0 # a non-standard amino acid, neither reward nor penalize
        
	orderIn3LetterCode_a = SequenceUtils.AALetter2OrderOf3LetterCode[a]

        if ord(b) < ord('A') or ord(b) > ord('Z'):
                print "a non standard amino acid: " + b
                return 0 # a non-standard amino acid
        
	orderIn3LetterCode_b = SequenceUtils.AALetter2OrderOf3LetterCode[b]

	return Ori_BLOSUM_80[orderIn3LetterCode_a, orderIn3LetterCode_b]
######################################## The End of BLOSUM Score ##########################################################

######################### score derived from sequence profile (profile HMM built by HHpred or HHblits)#####################
# calculate profile score by PSFM * PSSM, for two directions: seq-template and template-seq
def MutationOf2Pos5(tPos, sPos, temp, seq):
	"""
        if tPos < 0 or tPos >= temp['length']:
                print "MutationOf2Pos5: out of range by tPos=" + str(tPos) + " in template=" + temp['name']
                exit(1)
        
        if sPos < 0 or sPos >= seq['length']:
                print "MutationOf2Pos5: out of range by sPos=" + str(sPos) + " in sequence=" + seq['name']
                exit(1);
	"""

        m = np.inner(seq['PSFM'][sPos], temp['PSSM'][tPos])
        m += np.inner(temp['PSFM'][tPos], seq['PSSM'][sPos])

        return m

## calculate profile score by sequence PSFM * template PSSM
def MutationOf2Pos5_ST(tPos, sPos, temp, seq):
	"""
        if tPos < 0 or tPos >= temp['length']:
                print "MutationOf2Pos5: out of range by tPos=" + str(tPos) + " in template=" + temp['name']
                exit(1)
        
        if sPos < 0 or sPos >= seq['length']:
                print "MutationOf2Pos5: out of range by sPos=" + str(sPos) + " in sequence=" + seq['name']
                exit(1);
	"""
        m = np.inner(seq['PSFM'][sPos], temp['PSSM'][tPos])
        ##m += np.inner(temp['PSFM'][tPos], seq['PSSM'][sPos])

        return m

# similarity score between sequence and profile: two-way
def MutationOf2Pos6(tPos, sPos, temp, seq):
	"""
        if tPos < 0 or tPos >= temp['length']:
                print "MutationOf2Pos6: out of range by tPos=" + str(tPos) + " in template=" + temp['name']
                exit(1)

        if sPos < 0 or sPos >= seq['length']:
                print "MutationOf2Pos6: out of range by sPos=" + str(sPos) + " in sequence=" + seq['name']
                exit(1)
	"""

        tAA = temp['sequence'][tPos]    # the template residue at tPos
        sAA = seq['sequence'][sPos]     # the sequence residue at sPos

	## x and y are the order in 1-letter code
	x = SequenceUtils.AALetter2OrderOf1LetterCode[tAA]
	y = SequenceUtils.AALetter2OrderOf1LetterCode[sAA]

	#here we silently skip the non-standard amino acids in the sequence or template
	#Maybe we shall report a WARNING???
        m = 0
        if y >= 0 and y < 20:
                m += temp['PSSM'][tPos][y]  # // score between the template profile and sequence residue
        if x < 20 and x >= 0:
                m += seq['PSSM'][sPos][x]  # //score between the sequence profile and template residue

        return m


#// similarity score between the primary sequence and the profile, one-way: primary used for target and profile for template
def MutationOf2Pos6_ST(tPos, sPos, temp, seq):
	"""
        if tPos < 0 or tPos >= temp['length']:
                print "MutationOf2Pos6: out of range by tPos=" + str(tPos) + " in template=" + temp['name']
                exit(1)
        
        if sPos < 0 or sPos >= seq['length']:
                print "MutationOf2Pos6: out of range by sPos=" + str(sPos) + " in sequence=" + seq['name']
                exit(1)
	"""

        sAA = seq['sequence'][sPos]   # the sequence residue at sPos
        y = SequenceUtils.AALetter2OrderOf1LetterCode[sAA]

#        //here we silently skip the non-standard amino acids in the sequence or template
#        //Maybe we shall report a WARNING???
        if y >= 0 and y < 20:
		m = temp['PSSM'][tPos][y] #  score between the template profile and sequence residue
	else:
		m = 0.
        return m

######################################################### The End of Profile Score ##########################################################

#------------- below are two amino acid scoring matrices for very weak similarity. One is derived from structure alignment and the other from amino acid properties  ----------------#

Ori_HDSM = [
[  2.09,  -0.50,  -0.57,  -0.73,   0.33,  -0.75,  -0.12,   0.27,  -1.42,  -0.97,  -0.39,  -0.38,  -0.04,  -0.76,  -0.53,   0.34,   0.13,  -0.66,  -1.25,   0.02], #A
[ -0.50,   2.87,   0.60,   0.13,  -1.30,   0.13,   0.99,  -0.96,   0.54,  -1.40,  -1.19,   1.42,  -0.63,  -1.40,   0.21,  -0.06,  -0.15,  -0.04,  -0.75,  -1.52], #R
[ -0.57,   0.60,   3.60,   1.78,  -2.08,   0.33,  -0.16,   0.79,   0.76,  -2.43,  -2.10,   0.83,  -2.01,  -2.25,  -1.10,   0.40,   0.30,  -2.89,  -0.36,  -2.17], #N
[ -0.73,   0.13,   1.78,   4.02,  -2.51,   0.34,   1.20,  -1.20,  -0.01,  -2.77,  -2.65,   0.66,  -2.58,  -2.19,   0.72,   0.71,  -0.75,  -1.91,  -1.21,  -2.02], #D
[  0.33,  -1.30,  -2.08,  -2.51,   6.99,  -0.83,  -1.97,  -2.11,  -1.50,   0.13,  -0.31,  -2.19,   1.04,   1.13,  -2.19,   0.31,  -0.59,  -0.76,   0.13,   0.34], #C
[ -0.75,   0.13,   0.33,   0.34,  -0.83,   2.60,   1.23,  -0.12,  -0.46,  -1.47,  -1.49,   0.92,  -0.13,  -2.31,   0.24,   1.04,   0.60,  -0.81,  -0.61,  -1.38], #Q
[ -0.12,   0.99,  -0.16,   1.20,  -1.97,   1.23,   2.97,  -0.41,  -0.62,  -1.81,  -2.11,   1.11,  -1.86,  -1.61,  -0.26,   0.31,  -0.21,  -2.70,  -1.64,  -1.84], #E
[  0.27,  -0.96,   0.79,  -1.20,  -2.11,  -0.12,  -0.41,   4.36,  -0.40,  -2.93,  -1.98,  -0.71,  -1.86,  -2.67,  -0.04,   0.29,  -0.81,  -1.21,  -1.62,  -1.96], #G
[ -1.42,   0.54,   0.76,  -0.01,  -1.50,  -0.46,  -0.62,  -0.40,   5.89,  -1.76,  -0.93,   0.31,  -1.04,  -0.22,  -1.44,  -0.74,  -0.52,  -1.48,  -0.12,  -0.35], #H
[ -0.97,  -1.40,  -2.43,  -2.77,   0.13,  -1.47,  -1.81,  -2.93,  -1.76,   2.76,   1.56,  -1.81,   0.99,   0.76,  -2.00,  -1.75,  -0.96,   0.25,   0.08,   1.94], #I
[ -0.39,  -1.19,  -2.10,  -2.65,  -0.31,  -1.49,  -2.11,  -1.98,  -0.93,   1.56,   2.43,  -1.96,   1.61,   1.23,  -1.56,  -2.30,  -0.86,  -0.14,   0.70,   0.81], #L
[ -0.38,   1.42,   0.83,   0.66,  -2.19,   0.92,   1.11,  -0.71,   0.31,  -1.81,  -1.96,   2.91,  -1.62,  -2.41,  -0.19,  -0.06,  -0.10,  -1.94,  -1.72,  -1.27], #K
[ -0.04,  -0.63,  -2.01,  -2.58,   1.04,  -0.13,  -1.86,  -1.86,  -1.04,   0.99,   1.61,  -1.62,   3.75,   0.80,  -1.09,  -1.34,  -1.58,   0.87,  -0.41,   0.61], #M
[ -0.76,  -1.40,  -2.25,  -2.19,   1.13,  -2.31,  -1.61,  -2.67,  -0.22,   0.76,   1.23,  -2.41,   0.80,   3.28,  -0.91,  -1.11,  -0.69,   2.29,   1.96,   0.51], #F
[ -0.53,   0.21,  -1.10,   0.72,  -2.19,   0.24,  -0.26,  -0.04,  -1.44,  -2.00,  -1.56,  -0.19,  -1.09,  -0.91,   5.45,  -0.29,   0.93,  -5.34,  -1.98,  -1.11], #P
[  0.34,  -0.06,   0.40,   0.71,   0.31,   1.04,   0.31,   0.29,  -0.74,  -1.75,  -2.30,  -0.06,  -1.34,  -1.11,  -0.29,   2.36,   1.20,  -1.18,  -1.56,  -1.11], #S
[  0.13,  -0.15,   0.30,  -0.75,  -0.59,   0.60,  -0.21,  -0.81,  -0.52,  -0.96,  -0.86,  -0.10,  -1.58,  -0.69,   0.93,   1.20,   2.04,  -0.57,  -0.41,   0.05], #T
[ -0.66,  -0.04,  -2.89,  -1.91,  -0.76,  -0.81,  -2.70,  -1.21,  -1.48,   0.25,  -0.14,  -1.94,   0.87,   2.29,  -5.34,  -1.18,  -0.57,   6.96,   2.15,  -1.09], #W
[ -1.25,  -0.75,  -0.36,  -1.21,   0.13,  -0.61,  -1.64,  -1.62,  -0.12,   0.08,   0.70,  -1.72,  -0.41,   1.96,  -1.98,  -1.56,  -0.41,   2.15,   3.95,   0.21], #Y
[  0.02,  -1.52,  -2.17,  -2.02,   0.34,  -1.38,  -1.84,  -1.96,  -0.35,   1.94,   0.81,  -1.27,   0.61,   0.51,  -1.11,  -1.11,   0.05,  -1.09,   0.21,   2.05]] #V
#  A       R       N       D       C       Q       E       G       H       I       L       K       M       F       P       S       T       W       Y       V       

Ori_HDSM_core = np.array(Ori_HDSM, dtype=np.float32)
Ori_HDSM = np.full((21, 21), -3.5, dtype=np.float32)
Ori_HDSM[:20, :20] = Ori_HDSM_core

def HDSM(a, b):
	orderIn3LetterCode_a = SequenceUtils.AALetter2OrderOf3LetterCode[a]
	orderIn3LetterCode_b = SequenceUtils.AALetter2OrderOf3LetterCode[b]
	return Ori_HDSM[orderIn3LetterCode_a, orderIn3LetterCode_b]

newHDSM = np.zeros((21, 21), dtype=np.float32)
newHDSM[:20, :20] = Ori_HDSM_core


Ori_CC50 = [
[1.000,0.620,0.257,0.133,0.411,0.684,0.681,0.368,0.361,0.481,0.728,0.586,0.735,0.488,0.431,0.355,0.440,0.475,0.479,0.484], #A
[0.620,1.000,0.407,0.283,0.303,0.862,0.630,0.410,0.679,0.203,0.525,0.929,0.739,0.683,0.399,0.689,0.655,0.659,0.594,0.208], #R
[0.257,0.407,1.000,0.872,0.325,0.465,0.404,0.463,0.692,0.146,0.134,0.410,0.284,0.197,0.411,0.790,0.403,0.320,0.315,0.155], #N
[0.133,0.283,0.872,1.000,0.215,0.383,0.466,0.362,0.674,0.034,0.026,0.535,0.155,0.076,0.302,0.707,0.300,0.202,0.201,0.046], #D
[0.411,0.303,0.325,0.215,1.000,0.389,0.255,0.414,0.861,0.365,0.356,0.275,0.428,0.651,0.426,0.825,0.742,0.696,0.673,0.368], #C
[0.684,0.862,0.465,0.383,0.389,1.000,0.740,0.477,0.711,0.270,0.710,0.859,0.626,0.640,0.465,0.475,0.711,0.657,0.655,0.278], #Q
[0.681,0.630,0.404,0.466,0.255,0.740,1.000,0.370,0.362,0.157,0.401,0.743,0.245,0.184,0.363,0.396,0.359,0.288,0.297,0.168], #E
[0.368,0.410,0.463,0.362,0.414,0.477,0.370,1.000,0.468,0.282,0.273,0.394,0.396,0.333,0.465,0.470,0.455,0.423,0.412,0.288], #G
[0.361,0.679,0.692,0.674,0.861,0.711,0.362,0.468,1.000,0.276,0.265,0.416,0.387,0.567,0.454,0.930,0.704,0.661,0.738,0.285], #H
[0.481,0.203,0.146,0.034,0.365,0.270,0.157,0.282,0.276,1.000,0.499,0.168,0.465,0.491,0.364,0.263,0.380,0.443,0.449,0.948], #I
[0.728,0.525,0.134,0.026,0.356,0.710,0.401,0.273,0.265,0.499,1.000,0.606,0.791,0.741,0.358,0.252,0.371,0.688,0.447,0.498], #L
[0.586,0.929,0.410,0.535,0.275,0.859,0.743,0.394,0.416,0.168,0.606,1.000,0.512,0.444,0.363,0.687,0.643,0.292,0.307,0.174], #K
[0.735,0.739,0.284,0.155,0.428,0.626,0.245,0.396,0.387,0.465,0.791,0.512,1.000,0.815,0.445,0.366,0.438,0.821,0.487,0.466], #M
[0.488,0.683,0.197,0.076,0.651,0.640,0.184,0.333,0.567,0.491,0.741,0.444,0.815,1.000,0.400,0.550,0.650,0.918,0.919,0.492], #F
[0.431,0.399,0.411,0.302,0.426,0.465,0.363,0.465,0.454,0.364,0.358,0.363,0.445,0.400,1.000,0.449,0.463,0.470,0.470,0.370], #P
[0.355,0.689,0.790,0.707,0.825,0.475,0.396,0.470,0.930,0.263,0.252,0.687,0.366,0.550,0.449,1.000,0.791,0.646,0.724,0.271], #S
[0.440,0.655,0.403,0.300,0.742,0.711,0.359,0.455,0.704,0.380,0.371,0.643,0.438,0.650,0.463,0.791,1.000,0.699,0.904,0.387], #T
[0.475,0.659,0.320,0.202,0.696,0.657,0.288,0.423,0.661,0.443,0.688,0.292,0.821,0.918,0.470,0.646,0.699,1.000,0.742,0.444], #W
[0.479,0.594,0.315,0.201,0.673,0.655,0.297,0.412,0.738,0.449,0.447,0.307,0.487,0.919,0.470,0.724,0.904,0.742,1.000,0.452], #Y
[0.484,0.208,0.155,0.046,0.368,0.278,0.168,0.288,0.285,0.948,0.498,0.174,0.466,0.492,0.370,0.271,0.387,0.444,0.452,1.000]] #V
#A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V      

Ori_CC50_core = np.array(Ori_CC50, dtype=np.float32) - 0.5
Ori_CC50 = np.full((21,21), -1, dtype=np.float32)
Ori_CC50[:20, :20] = Ori_CC50_core

def CC50(a, b):
	orderIn3LetterCode_a = SequenceUtils.AALetter2OrderOf3LetterCode[a]
	orderIn3LetterCode_b = SequenceUtils.AALetter2OrderOf3LetterCode[b]
	return Ori_CC50[orderIn3LetterCode_a, orderIn3LetterCode_b]

newCC50 = np.zeros((21, 21), dtype=np.float32)
newCC50[:20, :20] = Ori_CC50_core

#============= the below matrix is not directly used for Similarity Score ===========#
Ori_GONNET = [
[   24,   -6,   -3,   -3,    5,   -2,    0,    5,   -8,   -8,  -12,   -4,   -7,  -23,    3,   11,    6,  -36,  -22,    1],  #A
[   -6,   47,    3,   -3,  -22,   15,    4,  -10,    6,  -24,  -22,   27,  -17,  -32,   -9,   -2,   -2,  -16,  -18,  -20],  #R
[   -3,    3,   38,   22,  -18,    7,    9,    4,   12,  -28,  -30,    8,  -22,  -31,   -9,    9,    5,  -36,  -14,  -22],  #N
[   -3,   -3,   22,   47,  -32,    9,   27,    1,    4,  -38,  -40,    5,  -30,  -45,   -7,    5,    0,  -52,  -28,  -29],  #D
[    5,  -22,  -18,  -32,  115,  -24,  -30,  -20,  -13,  -11,  -15,  -28,   -9,   -8,  -31,    1,   -5,  -10,   -5,    0],  #C
[   -2,   15,    7,    9,  -24,   27,   17,  -10,   12,  -19,  -16,   15,  -10,  -26,   -2,    2,    0,  -27,  -17,  -15],  #Q
[    0,    4,    9,   27,  -30,   17,   36,   -8,    4,  -27,  -28,   12,  -20,  -39,   -5,    2,   -1,  -43,  -27,  -19],  #E
[    5,  -10,    4,    1,  -20,  -10,   -8,   66,  -14,  -45,  -44,  -11,  -35,  -52,  -16,    4,  -11,  -40,  -40,  -33],  #G
[   -8,    6,   12,    4,  -13,   12,    4,  -14,   60,  -22,  -19,    6,  -13,   -1,  -11,   -2,   -3,   -8,   22,  -20],  #H
[   -8,  -24,  -28,  -38,  -11,  -19,  -27,  -45,  -22,   40,   28,  -21,   25,   10,  -26,  -18,   -6,  -18,   -7,   31],  #I
[  -12,  -22,  -30,  -40,  -15,  -16,  -28,  -44,  -19,   28,   40,  -21,   28,   20,  -23,  -21,  -13,   -7,    0,   18],  #L
[   -4,   27,    8,    5,  -28,   15,   12,  -11,    6,  -21,  -21,   32,  -14,  -33,   -6,    1,    1,  -35,  -21,  -17],  #K
[   -7,  -17,  -22,  -30,   -9,  -10,  -20,  -35,  -13,   25,   28,  -14,   43,   16,  -24,  -14,   -6,  -10,   -2,   16],  #M
[  -23,  -32,  -31,  -45,   -8,  -26,  -39,  -52,   -1,   10,   20,  -33,   16,   70,  -38,  -28,  -22,   36,   51,    1],  #F
[    3,   -9,   -9,   -7,  -31,   -2,   -5,  -16,  -11,  -26,  -23,   -6,  -24,  -38,   76,    4,    1,  -50,  -31,  -18],  #P
[   11,   -2,    9,    5,    1,    2,    2,    4,   -2,  -18,  -21,    1,  -14,  -28,    4,   22,   15,  -33,  -19,  -10],  #S
[    6,   -2,    5,    0,   -5,    0,   -1,  -11,   -3,   -6,  -13,    1,   -6,  -22,    1,   15,   25,  -35,  -19,    0],  #T
[  -36,  -16,  -36,  -52,  -10,  -27,  -43,  -40,   -8,  -18,   -7,  -35,  -10,   36,  -50,  -33,  -35,  142,   41,  -26],  #W
[  -22,  -18,  -14,  -28,   -5,  -17,  -27,  -40,   22,   -7,    0,  -21,   -2,   51,  -31,  -19,  -19,   41,   78,  -11],  #Y
[    1,  -20,  -22,  -29,    0,  -15,  -19,  -33,  -20,   31,   18,  -17,   16,    1,  -18,  -10,    0,  -26,  -11,   34]]  #V
#   A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V      

Ori_GONNET = np.array(Ori_GONNET, dtype=np.int32)
newGONNET = np.zeros((21,21), dtype=np.int32)
newGONNET[:20, :20] = Ori_GONNET

########################################## Start of Secondary Structure Mutation Score ################################

## note that in the tgt file, secondary structure is arranged in the order of helix, sheet and loop. 
## it is important to make them consistent 

SSMutation = [
[0.941183,  -2.32536,  -0.87487],  #HELIX
[-2.11462,  1.41307,   -0.401386], #SHEET
[-0.760861, -0.540041, 0.269711]]  #LOOP
# HELIX    SHEET     LOOP

SSMutation = np.array(SSMutation, dtype=np.float32)

#8-state secondary structure mutation score
##from template to the predicted SS8 of the target
SS8Mutation = [
[0.9823,        -0.1944, -1.905, -0.5508, -1.051, -1.163], #H(I)
[-0.07923,1.139,        -0.7431, 0.2849,  0.07965,-0.1479], #G
[-1.868, -0.7317, 1.274, -0.7433, -0.2456, -0.07621], #E(B)
[-0.4469, 0.2968,-0.8554, 0.9231, 0.2446,  -0.1803], #T
[-1.064,  0.0251, -0.3282,0.3049, 0.6813,  0.2468 ], #S
[-1.327, -0.3154, -0.2324, -0.2839, 0.1512, 0.3150]] #L
#H(I)    G       E(B)    T        S        L

SS8Mutation = np.array(SS8Mutation, dtype=np.float32)

HELIX = 0
SHEET = 1
LOOP = 2

SS8Letter2SS3Code = {'H':HELIX, 'G':LOOP, 'I':LOOP, 'E':SHEET, 'B':LOOP, 'T':LOOP, 'S':LOOP, 'L':LOOP }

# calculate secondary structure mutation score
def SSMutationScore_3State(tPos, sPos, temp, seq):
        if tPos < 0 or tPos >= temp['length']:
                print "SSMutationScore_3State: out of range by tPos=" + str(tPos) + " in template=" + temp['name']
                exit(-1)

        if sPos < 0 or sPos >= seq['length']:
                print "SSOf2Pos1: out of range by sPos=" + str(sPos) + " in sequence=" + seq['name']
                exit(-1)

        ss_type = temp['SS_str'][tPos]

        if ss_type not in SS8Letter2SS3Code:
                print "SSMutationScore_3State: unknown secondary structure type at position " + str(tPos) + " in template " + temp['name']
                return 0 # shall we terminate the program or just skip it ???

	ss_type = SS8Letter2SS3Code[ss_type]

	score = np.dot(SSMutation[ss_type], seq['SS3'][sPos])
	return score


SS82SS6 = { 'H':0, 'G':1, 'I':0, 'E':2, 'B':2, 'T':3, 'S':4, 'L':5, 'C':5 }
SS82SS6_n = { 0:0, 1:1, 2:0, 3:2, 4:2, 5:3, 6:4, 7:5, 8:5 }

## it is better not to use this function
def SSMutationScore_6State(tPos, sPos, temp, seq):
        if tPos < 0 or tPos >= temp['length']:
		print "SSOf2Pos2: out of range by tPos=" + str(tPos) + " in template=" + temp['name']
		sys.exit(-1)

        if sPos < 0 or sPos >= seq['length']:
		print "SSOf2Pos2: out of range by sPos=" + str(sPos) + " in sequence=" + seq['name']
		sys.exit(-1)

	tSS6 = SS82SS6[temp['SS_str'][tPos] ]

        score = 0.0
        for j in range(8): 
                #//sum over the secondary structure mutation score by predicted probability
                sSS6 = SS82SS6_n[j]
                score += SS8Mutation[tSS6][sSS6] * seq['SS8'][sPos][j] 
                
        return score

############################################## End of Secondary Structure Score #################################

############################################### Start of Solvent Accessibility Mutation score #################################

## in the tgt file, the predicted solvent accessibility is arranged in the order of buried, medium and exposed
## in the tpl file, the native solvent accessibility is arranged in the order of buried, medium and exposed

ACCMutation = [
[0.760885,  -0.0701501,  -0.903965],  #buried
[ -0.0798508,  0.218512,  -0.0623839], #medium
[-1.14008,  -0.099655,  0.233613]] #exposed
# buried     medium      exposed 

ACCMutation = np.array(ACCMutation, dtype=np.float32)

def ACCMutationScore_3State(tPos, sPos, temp, seq):
        if tPos < 0 or tPos >= temp['length']:
                print "ACC_New_Score2: out of range by tPos=" + str(tPos) + " in template=" + temp['name']
                exit(1)
        if sPos < 0 or sPos >= seq['length']:
                print "ACC_New_Score2: out of range by sPos=" + str(sPos) + " in sequence=" + seq['name']
                exit(1)

        tACC = temp['ACC'][tPos]
        if tACC > 2 or tACC < 0: 
                print "ACC_New_Score2: unknown solvent accessibility status at position " + str(tPos) + "of template " + temp['name']
                exit(1)

        sACC = seq['ACC'][sPos]

        if sACC > 2 or sACC < 0:
                print "ACC_New_Score2: Unknown solvent accessibility status at position " + str(tPos) + "of sequence " + seq['name']
                exit(1)

	score = np.dot(ACCMutation[tACC], seq['ACC_prob'][sPos])

        return score

############################ End of Solvent Accessibility Mutation Score ##################################################

############################ Start of Singleton Score ########################################################################     
Ori_Singleton = [
[  -60,  -13,  -18,    0,   58,   94,    1,   21,   33],   #A
[   99,  -54,  -53,  129,  -34,  -58,   94,   -3,   -5],   #R
[   87,   13,    2,   88,   25,    5,    6,  -31,  -46],   #N
[  111,   21,  -38,  117,   33,    7,   34,  -20,  -55],   #D
[  -33,   41,  191,  -66,    6,  121,  -70,  -23,  124],   #C
[  107,  -32,  -97,  146,   12,  -57,  126,   43,  -29],   #Q
[   67,  -34,  -77,  101,    4,  -51,   83,   22,  -16],   #E
[   47,  104,  116,   19,   58,  101,  -50,  -25,  -54],   #G
[   51,  -18,   18,   33,  -35,   -3,    3,  -30,    4],   #H
[  -58,    5,  119,  -89,  -20,   50,  -18,   38,  135],   #I
[  -75,  -23,   91,  -41,   19,   91,  -26,   15,  120],   #L
[  184,   -7,  -91,  210,   -4,  -90,  196,   45,  -51],   #K
[  -65,  -20,   74,  -27,   20,   65,  -25,   11,   74],   #M
[  -47,    9,  139,  -64,  -19,   79,  -44,   -3,  124],   #F
[  126,   91,   46,  131,   76,   20,  -35,  -46,  -63],   #P
[   37,   31,    8,   32,    6,   -8,  -14,  -19,  -25],   #S
[   29,   19,   29,   17,  -38,  -60,    1,   -8,   -9],   #T
[  -33,  -34,  123,  -26,  -48,   65,  -32,  -15,  123],   #W
[   -9,  -26,   97,  -26,  -70,   28,   -7,  -22,   93],   #Y
[  -40,   19,  109,  -93,  -35,   10,   -4,   32,  101]]   #V
# bury  medi  expo  bury  medi  expo   bury medi  expo
#         HELIX            SHEET            LOOP

Ori_Singleton = np.array(Ori_Singleton, dtype=np.int32)
newSingleton = np.zeros((21, 9), dtype=np.int32)
newSingleton[:20, ] = Ori_Singleton

#==== TMscore > 0.35 ======#
##singleton matrix (substition matrix from 17W CLEF alignment) -> (more positive, better)
WS_Singleton = [
[ 0.4521,  0.1993,  0.1440, -0.1481, -0.4385, -0.5268, -0.0161, -0.1930, -0.2784],  #A
[-0.3654,  0.3570,  0.4287, -0.6435,  0.2095,  0.2655, -0.3786,  0.0195,  0.1137],  #R
[-0.5198, -0.0223,  0.1781, -0.6236, -0.1045,  0.0789,  0.0198,  0.2663,  0.3942],  #N
[-0.6642,  0.0034,  0.3186, -0.6616, -0.1438,  0.0117, -0.0439,  0.2922,  0.3959],  #D
[ 0.1426, -0.4187, -0.9984,  0.4265,  0.0323, -0.3947,  0.4247,  0.0587, -0.5173],  #C
[-0.2788,  0.3363,  0.5703, -0.6899,  0.0179,  0.2306, -0.4087, -0.0564,  0.1450],  #Q
[-0.4877,  0.3586,  0.7434, -0.8903, -0.0308,  0.3240, -0.5427, -0.1179,  0.2259],  #E
[-0.4438, -0.5533, -0.5185, -0.1935, -0.2888, -0.2710,  0.3172,  0.3411,  0.4983],  #G
[-0.2226,  0.0608, -0.1088, -0.1423,  0.1681,  0.0601,  0.1487,  0.1427, -0.0238],  #H
[ 0.3407, -0.2074, -0.6775,  0.6630,  0.0352, -0.3408,  0.0132, -0.3905, -0.7773],  #I
[ 0.5600,  0.0477, -0.3708,  0.2819, -0.2322, -0.5238,  0.0296, -0.2958, -0.6507],  #L
[-0.5763,  0.2980,  0.5645, -0.9008,  0.1691,  0.4039, -0.5606,  0.0123,  0.3085],  #K
[ 0.4432,  0.0740, -0.2794,  0.1453, -0.2255, -0.4119,  0.0810, -0.1618, -0.4200],  #M
[ 0.3000, -0.1884, -0.7094,  0.4458,  0.0652, -0.3039,  0.1840, -0.1604, -0.5609],  #F
[-0.7645, -0.4526, -0.2861, -0.6109, -0.1524,  0.2258,  0.3069,  0.3887,  0.5543],  #P
[-0.3276, -0.1341,  0.0339, -0.3139,  0.0630,  0.2579,  0.0899,  0.1754,  0.2426],  #S
[-0.2333, -0.1518, -0.2018, -0.0647,  0.3434,  0.5032,  0.0151,  0.1004,  0.0398],  #T
[ 0.2079, -0.0059, -0.4754,  0.2423,  0.1374, -0.1509,  0.1852, -0.0916, -0.4484],  #W
[ 0.0997, -0.0204, -0.4840,  0.2422,  0.2807, -0.1224,  0.0971, -0.0134, -0.3813],  #Y
[ 0.1356, -0.3279, -0.7181,  0.7086,  0.2294, -0.0558,  0.0044, -0.3293, -0.6776]]  #V 
# bury     medi     expo     bury     medi     expo     bury     medi     expo
#         HELIX                       SHEET                      LOOP

WS_Singleton = np.array(WS_Singleton, dtype=np.float32)
newWSSingleton= np.zeros((21, 9), dtype=np.float32)
newWSSingleton[:20, ] = WS_Singleton

# [singleton], sequence profile of the target is used here
def SingletonScore_ProfileBased(tPos, sPos, temp, seq):
        if tPos < 0 or tPos >= temp['length']:
                print "ContactCapacityOf2Pos_old: out of range by tPos=" + str(tPos) + " in template=" + temp['name']
                exit(1)
        
        if sPos < 0 or sPos >= seq['length']:
                print "ContactCapacityOf2Pos_old: out of range by sPos=" + str(sPos) + " in sequence=" + seq['name']
                exit(1)

        ss = temp['SS_str'][tPos]

        if ss not in SS8Letter2SS3Code:
                print "Unknown secondary structure type at position " + str(tPos) + " of template " + temp['name']
                exit(1)
        
        ac = temp['ACC'][tPos]
        if ac < 0 or ac > 2:
                print "Unknown solvent accessibility type at position " + str(tPos) + " of template " + temp['name']
                exit(1)
        
	ss = SS8Letter2SS3Code[ss]
	score = np.dot(seq['PSFM'][sPos], Ori_Singleton[:,ss*3+ac])
        return score

# [ws_singleton], sequence profile of the target not used here, calculate the fitness score of one amino acid in a specific enviroment described by a combination of secondary structure and ACC
def SingletonScore_WS(tPos, sPos, temp, seq):
        if tPos < 0 or tPos >= temp['length']:
                print "ContactCapacityOf2Pos_WS: out of range by tPos=" + str(tPos) + " in template=" + temp['name']
                exit(1)
        
        if sPos < 0 or sPos >= seq['length']:
                print "ContactCapacityOf2Pos_WS: out of range by sPos=" + str(sPos) + " in sequence=" + seq['name']
                exit(1)

        ss = temp['SS_str'][tPos]
        if ss not in SS8Letter2SS3Code:
                print "Unknown secondary structure type at position " + str(tPos) + " of template " + temp['name']
                exit(1)
        
        ac = temp['ACC'][tPos]
        if ac < 0 or ac > 2:
                print "Unknown solvent accessibility type at position " + str(tPos) + " of template " + temp['name']
                exit(1)

        ss = SS8Letter2SS3Code[ss]

        res = seq['sequence'][sPos]
	
	j = SequenceUtils.AALetter2OrderOf3LetterCode[res]
	
        if j >= 20 or j < 0:
		return 0

        return WS_Singleton[j][ss*3+ac]

############################# End of Singleton Score ################################################


###########################The below contact capacity score is not used ##############################
Ori_Contact = [
[ 195.9, 309.8, -59.8,   4.0, 174.6,  -9.7,  -3.2, 109.2,  22.7,  37.2,  57.8,  10.3,  15.9,  10.6,   2.7, -27.8, -23.8,   7.9, -51.4, -68.2,  13.7, -35.1,-102.0,  15.5,  -3.1,-129.0,  13.8],  #A
[ 168.5, 282.7, -65.9, -19.2,  77.2, -43.1, -58.4,   4.5, -23.1, -45.0, -33.0, -16.0, -16.3, -37.2,  13.4,  24.7, -30.0,  52.4,  83.4, -15.9,  94.2, 153.7,   4.3, 133.4, 237.9,  52.9, 156.8],  #R
[ 152.9, 195.5, -82.8, -35.9,  63.7, -35.8, -55.3,  14.2,  -9.1, -12.3, -17.5,  -2.8,  -0.8, -21.5,  19.2,  17.4, -17.0,  48.1,  43.7, -17.8,  71.7,  60.0, -32.9,  84.1, 112.3, -48.9,  83.1],  #N
[ 111.1, 185.9,-100.2, -72.2,  33.0, -49.5, -69.2, -15.8, -15.4,  -8.8, -27.5,   5.8,  31.4, -16.1,  34.9,  58.9, -13.7,  72.0,  95.3,  -0.4, 102.8, 128.8,  -1.1, 128.1, 150.0,   7.5, 133.1],  #D
[ 314.1, 568.7,  70.8, 148.4, 302.2,  87.6,  99.3, 191.6,  55.8,  48.4,  94.3,   9.0,  -0.7,  27.5, -24.4, -56.2, -13.2, -27.6, -78.9, -71.6, -34.8, -78.8,-115.2, -42.5, -53.7,-168.0, -64.1],  #C
[ 144.8, 245.8, -84.1, -42.4,  55.5, -50.1, -65.1,  -4.5, -16.1, -29.2, -28.7,  -5.8,   7.2, -25.4,  30.4,  37.5, -19.3,  63.7,  79.6, -16.5,  98.0, 106.4,  -0.8, 105.6, 167.9,   9.5, 126.5],  #Q
[  97.9, 205.4,-113.3, -74.0,  35.7, -58.6, -70.9, -25.7, -17.6, -17.8, -44.8,   9.4,  31.7, -23.8,  52.9,  68.6,  -4.0, 103.9, 123.0,  16.4, 139.2, 173.1,  28.0, 161.8, 227.3,  50.7, 168.1],  #E
[ 136.6, 204.3, -86.5,   5.1,  96.4, -31.0,  -1.2,  96.1,  12.8,  35.3,  65.3,  19.8,  26.1,  15.5,  28.9,  -9.4, -38.4,  32.7, -46.8, -67.9,  34.1, -47.6, -88.8,  32.7, -52.0,-105.8,  17.9],  #G
[ 156.3, 305.5, -57.9,   6.1,  89.0, -12.4, -39.4,  32.2, -11.0, -30.8,  -5.4, -12.9, -28.3, -29.0,  -0.2,  -8.0, -36.2,  19.9,  34.7, -35.2,  55.6,  81.4, -40.3,  67.0, 158.3, -29.2,  82.2],  #H
[ 257.9, 484.1,  15.5, 101.1, 232.0,  54.3,  39.6, 124.7,  40.4,  16.6,  53.2,  -9.1, -18.9,  11.7, -19.7, -47.8, -31.1, -16.0, -53.0, -69.8, -20.8, -51.3,-100.8, -20.4, -20.6,-131.7, -23.8],  #I
[ 287.8, 449.7,  14.6,  87.0, 242.9,  43.4,  26.3, 122.3,  27.0,   1.6,  44.7, -16.8, -32.5,  -2.5, -23.1, -48.9, -41.0, -15.3, -40.7, -67.1,  -9.9, -22.7, -91.8,  -5.4,  21.1,-115.7,  10.6],  #L
[ 110.0, 237.7,-100.8, -54.4,  42.8, -54.4, -70.0, -20.5, -25.2, -35.9, -44.9,  -7.5,   4.0, -39.2,  39.5,  64.6, -17.3,  95.2, 144.6,  19.3, 151.2, 227.1,  66.8, 200.6, 306.9, 108.8, 256.5],  #K
[ 281.4, 397.2, -32.9,  71.2, 182.6,  15.5,   5.0,  85.8,  12.0,  -7.5,  35.5,  -7.6, -30.0,  -5.4, -11.2, -36.7, -38.2,   2.8, -36.0, -67.4,   5.4, -11.8, -83.3,  21.1,  38.0, -96.8,  12.1],  #M
[ 239.4, 398.5,  18.2,  83.9, 196.1,  41.2,  24.8, 105.6,  13.5, -16.7,  37.4, -24.4, -47.6,  -8.0, -32.9, -51.1, -48.2, -19.3, -25.1, -75.3,   2.7,  21.2, -83.2,  35.2,  97.9, -72.4,  65.4],  #F
[  68.4, 258.6, -81.0, -79.6,  34.4, -25.1,   5.3,  -9.1,   5.4,  25.5, -15.2,   3.7,   1.4,  -6.7,  14.7,  13.9, -12.7,  33.1,  19.3, -21.3,  43.3,  39.8,  -8.0,  47.3,  47.4, -45.3,  54.1],  #P
[ 137.8, 219.6, -85.0, -36.9,  76.7, -38.6, -34.6,  21.0,  -1.4,   9.3,  -3.9,   6.8,  14.2, -15.7,  27.3,   2.3, -18.2,  45.6,  -0.9, -29.9,  56.1,  11.6, -46.1,  60.4,  36.5, -73.4,  55.4],  #S
[ 164.7, 297.9, -69.0,  -4.2, 113.5, -28.7, -25.3,  15.7,   2.1,  -8.4, -23.4,   0.7,   4.7, -27.8,  10.6,  -7.5, -14.7,  30.6,  -7.8, -23.4,  48.5,  -7.0, -36.4,  51.6,  14.8, -60.7,  35.4],  #T
[ 224.0, 555.6,   9.3,  68.1, 194.7,  34.6,  15.8,  81.4,  12.4, -38.9,  14.5, -24.2, -53.6, -14.3, -34.6, -37.9, -45.2,  -9.5,  -7.1, -63.0,   2.5,  68.0, -82.8,  34.3, 168.0, -54.1,  78.9],  #W
[ 246.7, 418.8,   7.2,  70.6, 179.0,  33.8,  13.5, 108.7,  11.2, -14.9,  34.8, -26.1, -50.1, -15.2, -34.9, -48.2, -44.9, -13.9, -17.3, -70.6,  15.4,  35.7, -83.2,  39.2, 102.5, -74.6,  79.6],  #Y
[ 238.3, 466.6,  -6.9,  85.9, 218.4,  37.9,  32.8, 114.0,  36.3,  30.7,  54.0,  -3.8,  -8.8,  14.4, -11.2, -44.3, -25.7, -11.7, -55.0, -63.4, -14.7, -58.9,-105.8, -19.3, -38.3,-137.9, -23.5]]  #V
# HELIX  SHEET  LOOP  HELIX  SHEET  LOOP    HELIX  SHEET  LOOP    HELIX  SHEET  LOOP   HELIX  SHEET  LOOP  HELIX  SHEET  LOOP    HELIX  SHEET  LOOP  HELIX  SHEET  LOOP   HELIX  SHEET  LOOP
#          0                    1                    2                    3                    4                    5                    6                    7                    8      -> contact number

Ori_Contact = np.array(Ori_Contact, dtype=np.float32)
########################################## End of Contact Capacity Scoring Matrix ######################################

## purely for test
if __name__ == '__main__':
	print BLOSUM62('A', 'A')
	print BLOSUM45('A', 'A')
	print BLOSUM80('A', 'A')

	from Common import LoadTPLTGT

	tmp = LoadTPLTGT.load_tpl("1pazA.tpl") 
	tgt = LoadTPLTGT.load_tgt("1pazA.tgt")

	x, y = 2, 10

	print MutationOf2Pos5(x, y, tmp, tgt) 
	print MutationOf2Pos5_ST(x, y, tmp, tgt) 
	print MutationOf2Pos6(x, y, tmp, tgt) 
	print MutationOf2Pos6_ST(x, y, tmp, tgt) 
	print SSMutationScore_3State(x, y, tmp, tgt)
	print SSMutationScore_6State(x, y, tmp, tgt)
