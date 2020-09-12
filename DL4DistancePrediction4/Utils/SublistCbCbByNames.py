import numpy as np
import sys
import os
import cPickle


if len(sys.argv) < 2:
	print 'Usage: python SublistCbCbByNames.py proteinList PKLfile1 PKLfile2 PKLfile3 ... '
	print '     This script outputs all records with protein names in proteinList. Only Cb-Cb atomDistMatrix is saved'
        print '     This script also outputs a subset of proteinList not in the PKL files'
	exit(1)

listFile = sys.argv[1]
PKLfiles = sys.argv[2:]

fh = open(listFile, 'r')
proteins = set([ line.strip() for line in list(fh) ])
fh.close()

All = []

for file in PKLfiles:
	fh = open(file, 'rb')
	data = cPickle.load(fh)
	fh.close()
	All.extend(data)

print 'In total loaded ', len(All), ' entries from PKL files'

#subset = [ d for d in All if d['name'] in proteins ]

subset = []
validKeys = set(['CbCb', 'name', 'seq4matrix', 'pdbseq'])

for d in All:
	if d['name'] not in proteins:
		continue
	distMatrix =  d['atomDistMatrix']
	newDistMatrix = { k : v for k, v in distMatrix.iteritems() if k in validKeys }
	d['atomDistMatrix' ] = newDistMatrix

	subset.append(d)
	

remains = proteins - set([d['name'] for d in subset])
remainStr = '\n'.join(remains)
filename = os.path.basename(listFile).split('.')[0] + '-unfilled.list'

fh = open(filename, 'w')
fh.write(remainStr)
fh.close()


print 'Extracted records for ', len(subset), ' proteins'

## write

savefile = os.path.basename(listFile).split('.')[0] + '.' + str(os.getpid() ) + '.distanceFeatures.pkl'
f= open(savefile, 'wb')
cPickle.dump(subset, f, cPickle.HIGHEST_PROTOCOL)
f.close()

