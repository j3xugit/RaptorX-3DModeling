import os
import sys
import numpy as np
import cPickle

import glob

def Usage():
	print 'python CheckMisMatches4TPLNPDB.py folder4TPLCoordinates'
	print '	this script checks the number of mismatched residues between a template (*.hhm and *.tpl files) and its corresponding PDB file'
	print '	the input is a folder containing all the TPL coordinate files. The TPL coordinates are extracted from PDB files by the TPL sequence '

if len(sys.argv) < 2:
	Usage()
	exit(1)


tplCoordinateDir=sys.argv[1]
filepattern = tplCoordinateDir + '/*.pkl'
allfiles = glob.glob(filepattern)

print 'in total there are ', len(allfiles), ' files ending with pkl in ', tplCoordinateDir

numMisMatches = []
i = 0
for f in allfiles:
	with open(f, 'rb') as fh:
		data = cPickle.load(fh)
		numMisMatches.append((f, data['numMisMatches']))

	if i%1000 == 1:
		print '#files loaded: ', i

for fname, mismatch in numMisMatches:
	print fname, mismatch 

numMisMatches = [ a[1] for a in numMisMatches ]

bins = np.arange(10)
result = np.histogram(numMisMatches, bins)
print result
