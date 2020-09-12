import os
import sys
import numpy as np

from Common.LoadTPLTGT import load_tpl as LoadTPL

def DiffOfTwoArrays(a1, a2):
	b1 = np.array(a1, dtype=np.float32)
	b2 = np.array(a2, dtype=np.float32)
	diff = abs(b1-b2)
	return np.mean(diff)

def HammingDiffOfTwoArrays(a1, a2):
        diff = np.array( (a1 != a2), dtype=np.int32)
        return np.mean(diff)

if len(sys.argv)<3:
	print 'python CheckConsistency .tpl .tpl.pkl'
	exit(1)

file1=sys.argv[1]
file2=sys.argv[2]

p1=LoadTPL(file1)
p2=LoadTPL(file2)

if p1['sequence'] != p2['sequence']:
	print 'inconsistent sequence'
	exit(1)

print 'SS3 diff: ', DiffOfTwoArrays(p1['SS3Coding'], p2['SS3Coding'])
print 'SS8 diff: ', DiffOfTwoArrays(p1['SS8Coding'], p2['SS8Coding'])
print 'pACC diff: ', DiffOfTwoArrays(p1['pACC'], p2['pACC'])
print 'ACC diff: ', HammingDiffOfTwoArrays(p1['ACC'], p2['ACC'])
print 'PSFM diff: ', DiffOfTwoArrays(p1['PSFM'], p2['PSFM'])
print 'PSSM diff: ', DiffOfTwoArrays(p1['PSSM'], p2['PSSM'])
print 'missing diff: ', HammingDiffOfTwoArrays(p1['missing'], p2['missing'])
