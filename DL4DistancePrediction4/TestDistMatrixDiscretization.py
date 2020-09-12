import os
import sys
import numpy as np

import config
import DistanceUtils
import cPickle

def Usage():
	print 'python TestDistMatrixDiscretization.py groundTruthFilePKL scheme'
	print '	scheme: identitication for discretization distance bins, e.g., 13C, 14C, 47CPlus'
	
if len(sys.argv)<3:
	Usage()
	exit(1)
gtfile = sys.argv[1]
scheme = sys.argv[2]

with open(gtfile, 'rb') as fh:
	data = cPickle.load(fh)

if not data.has_key('atomDistMatrix'):
	print 'ERROR: no distance matrix found in ', gtfile
	exit(1)

distMatrix = data['atomDistMatrix']

bins = config.distCutoffs[scheme]

for apt, matrix in distMatrix.iteritems():
	result, _, _ = DistanceUtils.DiscretizeDistMatrix(matrix, bins=bins, invalidDistanceSeparated = scheme.endswith('Plus') )
	print apt
	print np.diagonal(matrix)
	print np.diagonal(result)
	if apt == 'CaCa':
		print np.diagonal(matrix[:-1, 1:])
		print np.diagonal(result[:-1, 1:])
