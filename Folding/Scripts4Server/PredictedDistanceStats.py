import numpy as np
import sys
import os
import cPickle

distcbFile = sys.argv[1]

targetName = os.path.basename(distcbFile).split('.')[0]

distMatrix = np.loadtxt(distcbFile, dtype=np.float32)

thresholds = [ 8, 9, 10, 11, 12, 13, 14, 15 ]
count = []

M1s = np.ones_like(distMatrix, dtype = np.int8)
mask_LR = np.triu(M1s, 24)
mask_MLR = np.triu(M1s, 12)
mask_SMLR = np.triu(M1s, 6)
mask_MR = mask_MLR - mask_LR
mask_SR = mask_SMLR - mask_MLR

##for mask in [ mask_LR, mask_MR, mask_SR]:
for mask in [ mask_SR]:
	res = distMatrix[mask.nonzero()]

	for cutoff in thresholds:
		num = ( ( 0 < res ) & (res < cutoff) ).sum()
		count.append(num)

	## include the disordered regions here
	num = ( res < 0).sum() + count[-1]
	count.append(num)

ratio = count/np.float32(distMatrix.shape[0])

ratioStr = [ "{:.2f}".format(r) for r in ratio ]

line = '\t'.join( [ targetName, str(distMatrix.shape[0]) ] + ratioStr )

print line



