#!/usr/bin/env python

"""Convert .pkl file created by deep learning to .epad_prob
Copyright
Author: zjw@ttic.edu (Jianwei Zhu)
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys

try:
	import cPickle as pickle
except:
  	import pickle

if __name__ == '__main__':
	if len(sys.argv) != 3:
    		sys.exit('Usage: %s <XXX.pkl> <XXX.epad_prob>' % sys.argv[0])
  	pkl_file, epad_file = sys.argv[1:6]

  	# read in pkl file
  	with open(pkl_file, 'r') as fin:
    		name, sequence, distMatrix = pickle.load(fin)[:3]

	if not distMatrix.has_key('CbCb_Discrete12C'):
		sys.exit('ERROR: the distance matrix PKL file does not have a key CbCb_Discrete12C')
  	dist = distMatrix['CbCb_Discrete12C'] * 100
  	m, n, l = dist.shape
  	assert m == n and l == 12

  	# output epad
  	with open(epad_file, 'w') as fout:
    		for i in range(m):
      			for j in range(i+1, m):
        			print(i, j, ' '.join(['%.6f'%_ for _ in dist[i,j]]), file=fout)
