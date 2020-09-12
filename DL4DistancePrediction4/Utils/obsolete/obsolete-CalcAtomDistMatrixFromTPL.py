import numpy as np
import sys
import os
import cPickle

sys.path.append(os.path.join(os.environ['ModelingHome'], 'Common'))
from LoadTPLTGT import load_tpl as TPLReader

##this script reads in a tpl file and calculate the Ca-Ca, Ca-Cb and Cb-Cb dist matrices and save them as a dictionary in a PKL file.

if len(sys.argv)<2:
	print 'Usage: python CalcAtomDistMatrixFromTPL.py tplFile'
	exit(-1)

tplFile = sys.argv[1]

if not os.path.isfile(tplFile):
	print 'Cannot find the tpl file: ', tplFile
	exit(-1)

protein = TPLReader(tplFile)

if not protein.has_key('Ca') or not protein.has_key('Cb') or not protein.has_key('missing'):
	print 'This protein does not have Ca or Cb or missing information for dist matrix calculation'
	exit(-1)

## calclulate Ca-Ca distance matrix
dMatrix = []

for ca1 in zip(protein['Ca']):
	dVector = [ np.linalg.norm( ca1 - ca2) for ca2 in protein['Ca'] ]
	dMatrix.append(dVector)

CaMatrix = np.array(dMatrix).astype(np.float16)


## calculate Ca-Cb distance matrix
dMatrix = []

for ca1 in zip(protein['Ca']):
	dVector = [ np.linalg.norm( ca1 - cb2) for cb2 in protein['Cb'] ]
	dMatrix.append(dVector)

CaCbMatrix = np.array(dMatrix).astype(np.float16)

## calculate Cb-Cb distance matrix
dMatrix = []
for cb1 in zip(protein['Cb']):
	dVector = [ np.linalg.norm( cb1 - cb2) for cb2 in protein['Cb'] ]
	dMatrix.append(dVector)
CbMatrix = np.array(dMatrix).astype(np.float16)

length = protein['missing'].shape[0]
for i, missing in zip(range(length), protein['missing']):
	if missing:
		CaMatrix[i, :] = -1.
		CaMatrix[:, i] = -1.
		CaCbMatrix[i, :] = -1.
		CaCbMatrix[:, i] = -1. 
		CbMatrix[i, :] = -1.
		CbMatrix[:, i] = -1.

np.fill_diagonal(CaMatrix, 0)
np.fill_diagonal(CbMatrix, 0)

atomDistMatrix = dict()
atomDistMatrix['CaCa'] = CaMatrix
atomDistMatrix['CaCb'] = CaCbMatrix
atomDistMatrix['CbCb'] = CbMatrix

name = os.path.basename(tplFile).split('.')[0]
atomDistMatrix['name'] = name

savefile = name + '.atomDistMatrix.pkl'
fh = open(savefile, 'wb')
cPickle.dump(atomDistMatrix, fh, protocol=cPickle.HIGHEST_PROTOCOL)
fh.close()

"""
np.savetxt(name+'.CaCa.txt', CaMatrix, fmt='%.2f')
np.savetxt(name+'.CaCb.txt', CaCbMatrix, fmt='%.2f')
np.savetxt(name+'.CbCb.txt', CbMatrix, fmt='%.2f')
"""
