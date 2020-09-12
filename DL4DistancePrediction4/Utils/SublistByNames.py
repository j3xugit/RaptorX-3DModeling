import numpy as np
import sys
import os
import cPickle

if len(sys.argv) < 2:
	print 'Usage: python SublistByNames.py proteinList PKLfile1 PKLfile2 PKLfile3 ... '
	exit(1)

listFile = sys.argv[1]
PKLfiles = sys.argv[2:]

fh = open(listFile, 'r')
proteins = set([ line.strip() for line in list(fh) ])
fh.close()


print '#need proteins: ', len(proteins)

All = []

for file in PKLfiles:
	fh = open(file, 'rb')
	data = cPickle.load(fh)
	fh.close()
	subset = [ d for d in data if d['name'] in proteins ]
	All.extend(subset)

allnames = [ d['name'] for d in All ]
allnames = list(dict.fromkeys(allnames))
print 'In total loaded ', len(All), ' records for ', len(allnames), ' proteins.'

## write
savefile = os.path.basename(listFile).split('.')[0] + '.distanceFeatures.pkl'
f= open(savefile, 'wb')
cPickle.dump(All, f, cPickle.HIGHEST_PROTOCOL)
f.close()

