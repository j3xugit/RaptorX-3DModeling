import numpy as np
import sys
import os
import cPickle

distcbfile = sys.argv[1]
name = os.path.basename(distcbfile).split('.')[0]

distcbMatrix = np.loadtxt(distcbfile)

potentialFile = sys.argv[2]
fh = open(potentialFile, 'rb')
potMatrix = cPickle.load(fh)
fh.close()


sizePot = potMatrix.shape
sizeDist = distcbMatrix.shape


print 'number of residues in the distcb file: ', sizeDist[0]
print 'protein sequence length derived from the potential matrix: ', sizePot[0]

if sizePot[0] != sizeDist[0]:
	print 'Error: the number of residues in the distcb file is inconsistent with the size of the potential matrix'
	exit(-1)

energy = 0

for i in range(sizeDist[0]):
	for j in range(i+6, sizeDist[1]):

		distance = distcbMatrix[i, j]
		if distance < 0:
			distance = 16

		bin = 0
		if distance > 5:
			bin = max(11, np.int32(distance - 4) )

		if bin < 8:
			energy += potMatrix[i, j, bin]


print 'the energy of the protein structure for ' + name + ' is ', energy
