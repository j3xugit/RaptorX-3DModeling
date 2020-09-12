import cPickle
import sys

import numpy as np

class prettyfloat(float):
    def __repr__(self):
        return "%0.2f" % self

inputFile = sys.argv[1]


fh = open(inputFile, 'rb')
data = cPickle.load(fh)
fh.close

for pair in data:
	seq, temp, features, distMatrix, xyresidues, coordinates = pair
	print '\n'
	print seq, temp

	for (f, res,point) in zip(features, xyresidues, coordinates):
		if res[3] == -1:
			print map(prettyfloat, f), res[0], chr(res[1]+ord('A')), res[2], '-', point
		else:
			print map(prettyfloat, f), res[0], chr(res[1]+ord('A')), res[2], chr(res[3]+ord('A')), point
	
