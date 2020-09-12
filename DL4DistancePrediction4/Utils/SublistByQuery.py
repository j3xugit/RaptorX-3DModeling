import numpy as np
import sys
import os
import cPickle

if len(sys.argv) < 2:
	print 'Usage: python SublistByQuery.py PKLfile1 PKLfile2 PKLfile3 ...'
	exit(1)

PKLfiles = sys.argv[1:]

All = []

for file in PKLfiles:
	fh = open(file, 'rb')
	data = cPickle.load(fh)
	fh.close()
	All.extend(data)

print 'In total loaded ', len(All), ' entries'

SubsetByTarget=dict()

for d in All:
	target = d['queryName']
	if not SubsetByTarget.has_key(target):
		SubsetByTarget[target] = []
	SubsetByTarget[target].append(d)


## write
print 'In total find input features for ', len(SubsetByTarget.keys() ), ' targets'

for k, v in SubsetByTarget.iteritems():
	savefile = k + '.distanceFeatures.pkl'
	f= open(savefile, 'wb')
	cPickle.dump(v, f, cPickle.HIGHEST_PROTOCOL)
	f.close()

