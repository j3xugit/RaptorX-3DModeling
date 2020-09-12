import numpy as np
import sys
import os

if  len(sys.argv) < 2:
	print 'please provide a true distcb file'
	exit(-1)

trueDistFile = sys.argv[1]

trueDist = np.loadtxt(trueDistFile, dtype=np.float32)

size = trueDist.shape
for offset in range(6):
	rng=np.arange(0, size[0]-offset)
	trueDist[rng, rng+offset] = -1.

for offset in range(1, 6):
	rng=np.arange(offset, size[0])
	trueDist[rng, rng-offset] = -1.

fakeDist = np.copy(trueDist)
fakeLower = np.ones_like(trueDist) * (-1.0)
fakeUpper = np.ones_like(trueDist) * (-1.0)

np.putmask(fakeDist, trueDist>15, -1)
np.putmask(fakeLower, fakeDist>=0, 0.5)
np.putmask(fakeUpper, fakeDist>=0, 0.5)


fakeDistFile = os.path.basename(trueDistFile).split('.')[0] + '.dist'
fakeLowerFile = os.path.basename(trueDistFile).split('.')[0] + '.lower'
fakeUpperFile = os.path.basename(trueDistFile).split('.')[0] + '.upper'

np.savetxt(fakeDistFile, fakeDist, fmt='%.4f', delimiter=' ')
np.savetxt(fakeLowerFile, fakeLower, fmt='%.4f', delimiter=' ')
np.savetxt(fakeUpperFile, fakeUpper, fmt='%.4f', delimiter=' ')
