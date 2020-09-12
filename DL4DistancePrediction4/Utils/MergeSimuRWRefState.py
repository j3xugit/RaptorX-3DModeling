import numpy as np
import sys
import os
import cPickle
import random

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import config

def str_display(ls):
        if not isinstance(ls, (list, tuple, np.ndarray)):
                str_ls = '{0:.7f}'.format(ls)
                return str_ls

        str_ls = ['{0:.7f}'.format(v) for v in ls ]
        str_ls2 = ' '.join(str_ls)
        return str_ls2

def Usage():
	print 'python MergeRWRefState.py refState1 refState2 refState3 ...'

def ReduceBins(refState, bins):
	assert bins[0] == 0

	newRefState = []
	start = 0
	for b in bins[1:]:
		end = int(b*2)
		newProb = np.sum(refState[:, start:end], axis=1, keepdims=True)
		newRefState.append(newProb)
		start = end
	newRefState.append(np.sum(refState[:, start:], axis=1, keepdims=True) )

	return np.concatenate(newRefState, axis=-1)

if len(sys.argv) < 2:
	Usage()
	exit(1)

refFiles = sys.argv[1:]

refStates = []

for f in refFiles:
	fh = open(f, 'rb')
	d = cPickle.load(fh)
	fh.close()
	refStates.append(d)

refState = np.average(refStates, axis=0)

"""
for row in refState:
	print str_display(row)
"""

## convert refState to different discretization schemes
schemes = ['25C', '34C', '64C', '63C', '62C']

result = dict()
for s in schemes:
	bins = config.distCutoffs[s]
	response = 'CbCb_Discrete' + s
	result[response] = ReduceBins(refState, bins)

savefile = 'RWRefState-final.pkl'
fh = open(savefile, 'wb')
cPickle.dump(result, fh,  protocol=cPickle.HIGHEST_PROTOCOL)
fh.close()

for row in result['CbCb_Discrete34C']:
	print str_display(row)
