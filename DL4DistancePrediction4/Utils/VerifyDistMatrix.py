import numpy as np
import sys
import os
import cPickle

if len(sys.argv) < 2:
	print 'Usage: python VerifyDistMatrix.py target_name'
	exit(-1)

target = sys.argv[1]

distcbDir = 'Test500-distcb/'
tplDistDir = 'pdb25-tpl-atomDistMatrix/'
pdbDistDir = 'pdb25-7952-atomDistMatrix/'

distcbFile = distcbDir + target + '.distcb'
distcbMatrix = np.loadtxt(distcbFile, dtype=np.float16)

tplDistFile = tplDistDir + target + '.atomDistMatrix.pkl'
tplfh = open(tplDistFile, 'rb')
tplDistMatrix = cPickle.load(tplfh)
tplfh.close()

pdbDistFile = pdbDistDir + target + '.atomDistMatrix.pkl'
pdbfh = open(pdbDistFile, 'rb')
pdbDistMatrix = cPickle.load(pdbfh)
pdbfh.close()

length = pdbDistMatrix['CaCa'].shape[0]


## Ca-Ca matrix comparison

assert (distcbMatrix.shape == tplDistMatrix['CaCa'].shape )
assert (distcbMatrix.shape == pdbDistMatrix['CaCa'].shape )

tplDiff = np.max( np.triu(abs(distcbMatrix - tplDistMatrix['CbCb']), 1) )
pdbDiff = np.max( np.tril(abs(distcbMatrix - pdbDistMatrix['CbCb']), -1) )
##pdbDiff = np.argmax( np.tril(abs(distcbMatrix - pdbDistMatrix['CbCb']), -1) )
##pdbDiff = np.unravel_index(pdbDiff, distcbMatrix.shape)

Diff2tplpdb = np.histogram(abs(tplDistMatrix['CaCa'] - pdbDistMatrix['CaCa']) )
print Diff2tplpdb

CbDiff2tplpdb = np.max(abs(tplDistMatrix['CbCb'] - pdbDistMatrix['CbCb']) )


if tplDiff>0.1 or pdbDiff>0.1 or Diff2tplpdb>0.1 or CbDiff2tplpdb>0.1:
	##print 'WARNING: inconsistent dist matrix for target: ', target, tplDiff, pdbDiff, pdbDistMatrix['CbCb'][pdbDiff], Diff2tplpdb, CbDiff2tplpdb
	print 'WARNING: inconsistent dist matrix for target: ', target, tplDiff, pdbDiff, Diff2tplpdb, CbDiff2tplpdb
	"""
	np.savetxt(target+'.tpl.CaCa.txt', tplDistMatrix['CaCa'], fmt='%.3f')
	np.savetxt(target+'.pdb.CaCa.txt', pdbDistMatrix['CaCa'], fmt='%.3f')
	np.savetxt(target+'.tpl.CbCb.txt', tplDistMatrix['CbCb'], fmt='%.3f')
	"""
	np.savetxt(target+'.pdb.CbCb.txt', pdbDistMatrix['CbCb'], fmt='%.3f')

"""
##calculate correlation coefficient between different dist matrices

CaCbCoef = np.histogram( pdbDistMatrix['CaCa'].flatten()- pdbDistMatrix['CbCb'].flatten() )
CaCgCoef = np.histogram( pdbDistMatrix['CaCa'].flatten()- pdbDistMatrix['CgCg'].flatten() )
CaNOCoef = np.histogram( pdbDistMatrix['CaCa'].flatten()- pdbDistMatrix['NO'].flatten() )
CaCaCgCoef = np.histogram( pdbDistMatrix['CaCa'].flatten()- pdbDistMatrix['CaCg'].flatten() )

print CaCbCoef
print CaCgCoef
print CaNOCoef
print CaCaCgCoef
"""

"""
A = np.vstack([pdbDistMatrix['CaCa'].flatten(), np.ones(len(pdbDistMatrix['CaCa'].flatten()))]).T.astype(np.float32)
res=np.linalg.lstsq(A, pdbDistMatrix['CbCb'].flatten().astype(np.float32))

print res

np.savetxt(target +'.CaCg.txt', pdbDistMatrix['CaCg'], fmt='%.2f')
"""
