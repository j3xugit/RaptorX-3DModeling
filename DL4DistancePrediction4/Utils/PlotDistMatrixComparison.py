import matplotlib
matplotlib.use('Agg') # plot but not display

import matplotlib.pyplot as plt

import numpy as np
import sys
import cPickle

if len(sys.argv) < 4:
	print 'Usage: python plotDistMatrixComparison.py targetName predDistBoundFile groundTruth_PKL'
        print '\tpredDistBoundFile: the predicted distance bound file ending with .bound.pkl or bound.txt'
	print '\tgroundTurth_PKL: the ground truth file ending with .native.pkl that contains experimental distance matrix information'
	exit(1)

target = sys.argv[1]
predDistBoundFile = sys.argv[2]
nativeDistMatrixFile = sys.argv[3]

## load in distance bound file, return sequence and predicted dist matrix
def LoadDistBoundFile(boundFile):
	seq = ''
	predDistList = []
	with open(boundFile, 'r') as f:
		for line in f.readlines():
			if line.startswith("SEQ"):
				seq += line.strip().split()[1]
				continue

			if len(line)>0 and line[0].isdigit():
				fields = line.strip().split()
				if len(fields) < 3:
					print 'WARNING: incorrect format in line: ', line
					exit(1)
				i = np.int32(fields[0]) - 1
				j = np.int32(fields[1]) - 1
				d = np.float32(fields[2])
				predDistList.append((i, j, d))
				continue

	## convert dist list to matrix
	predDistMatrix = np.full( (len(seq), len(seq)), -1, dtype=np.float32)

	for e in predDistList:
		predDistMatrix[ e[0], e[1] ] = e[2]
		predDistMatrix[ e[1], e[0] ] = e[2]

	return predDistMatrix, seq

def LoadDistBoundFileInPKL(boundFile, atomPairType='CbCb'):
	f = open(boundFile, 'r')
	data, name, seq = cPickle.load(f)
	f.close()

	if not data.has_key(atomPairType):
                print 'ERROR: the distance bound file ', boundFile, ' does not have predicted distance information for ', atomPairType
                exit(1)

	return data['CbCb'][:,:,0], seq
	

def LoadNativeDistMatrix(distFile, atomPairType='CbCb'):
	with open(distFile, 'rb') as fh:
		data = cPickle.load(fh)
	atomDistMatrix = data['atomDistMatrix']

	if not atomDistMatrix.has_key(atomPairType):
		print 'ERROR: the native distance file ', distFile, ' does not has distance information for ', atomPairType
		exit(1)

	return atomDistMatrix['CbCb']

if predDistBoundFile.endswith('.pkl'):
	predDistMatrix, _ = LoadDistBoundFileInPKL(predDistBoundFile)
else:
	predDistMatrix, _ = LoadDistBoundFile(predDistBoundFile)
trueDistMatrix = LoadNativeDistMatrix(nativeDistMatrixFile)

offset = 0
distMatrix4Display = np.triu(predDistMatrix, offset) + np.tril(trueDistMatrix, -offset)
np.putmask(distMatrix4Display, distMatrix4Display>15, -1)

## do not display distance for residue pair (i, j) where |i-j| < 3
for offset in range(0, 3):
	np.fill_diagonal( distMatrix4Display[offset:, : -offset], -1)
	np.fill_diagonal( distMatrix4Display[: -offset, offset:], -1)

masked_array = np.ma.masked_where( distMatrix4Display<0, distMatrix4Display)

cmap = matplotlib.cm.spring  # Can be any colormap that you want after the cm
cmap.set_bad(color='white')

plt.imshow(masked_array, cmap=cmap)
plt.colorbar()
plt.title('Predicted and native distance matrix of ' + target)
#plt.show()

file4save = target + '.tiff'
plt.savefig(fname=file4save, format='tiff', dpi=300)
