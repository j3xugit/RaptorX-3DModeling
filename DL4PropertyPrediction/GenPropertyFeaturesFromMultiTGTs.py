import numpy as np
import sys
import os
import cPickle

import PropertyUtils
from Common.LoadTPLTGT import load_tgt as LoadTGT
from GenPropertyFeaturesFromTGT import LoadTestData4Properties

def Usage():
	print 'python GenPropertyFeaturesFromMultiTGTs.py proteinName tgtFile1 tgtFile2 tgtFile3 ...'
	print '\tThis script generates a property feature file (in PKL) from a set of tgt files of the same protein'
	print '\ttgtFiles: the tgt files generated for template-based modeling'

def main(argv):
	if len(argv) < 2:
		Usage()
		exit(1)

	proteinName = argv[0]
	tgtFiles = argv[1:]

	features = []
	for tgtFile in tgtFiles:
		if not os.path.isfile(tgtFile):
			print 'ERROR: invalid tgt file', tgtFile
			exit(1)
		
		feature = LoadTestData4Properties(tgtFile, proteinName)
		features.append(feature)

	savefile = proteinName + '.propertyFeatures.pkl'
	with open(savefile, 'wb') as fh:
		cPickle.dump(features, fh, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
        main(sys.argv[1:])

