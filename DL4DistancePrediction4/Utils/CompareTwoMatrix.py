import os
import sys
import numpy as np
import cPickle

if len(sys.argv) < 3:
	print 'python CompareTwoMatrix.py file1_PKL file2_PKL'
	exit(1)

file1=sys.argv[1]
file2=sys.argv[2]

with open(file1, 'rb') as fh:
	m1 = cPickle.load(fh)

with open(file2, 'rb') as fh:
	m2 = cPickle.load(fh)

	for apt in ['CaCa', 'CbCb', 'NO', 'CgCg', 'CaCg']:
		print 'atom pair: ', apt
		print 'shape of two matrices: ', m1[apt].shape, m2[apt].shape
		if not np.allclose(m1[apt], m2[apt], rtol=0.01, atol=0.01):
			print 'two matrices for ', apt, ' are not equal'
			diff = np.absolute(m1[apt] - m2[apt])
			print 'max diff: ', np.amax(diff), ' at position: ', np.unravel_index(np.argmax(diff), diff.shape)
			position = np.unravel_index(np.argmax(diff), diff.shape)
			print position
			print 'avg diff: ', np.mean(diff)
			print m1[apt][position]
			print m2[apt][position]

