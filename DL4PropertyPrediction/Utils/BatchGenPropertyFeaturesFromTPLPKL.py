import numpy as np
import sys
import os
import cPickle

from DL4PropertyPrediction import PropertyUtils
from GenPropertyFeaturesFromTPLPKL import LoadTrainData4Properties

def Usage():
	print 'python BatchGenPropertyFeaturesFromTPLPKL.py proteinListFile tplpklDir'
	print '	tplpklDir: the folder for tpl.pkl files generated for template-based modeling'

def main(argv):

	if len(argv) < 2:
		Usage()
		exit(1)

	proteinListFile = argv[0]
	tplDir = argv[1]

	with open(proteinListFile, 'r') as fh:
		names = [ line.strip() for line in list(fh) ]

	proteins = []
	for name in names:
		tplFile = os.path.join(tplDir, name + '.tpl.pkl')
		protein = LoadTrainData4Properties(tplFile)
		proteins.append(protein)

	savefile = os.path.basename(proteinListFile).split('.')[0] + '.' + os.getpid() + '.propertyFeatures.pkl'
	with open(savefile, 'wb') as fh:
		cPickle.dump( proteins, fh, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
        main(sys.argv[1:])

