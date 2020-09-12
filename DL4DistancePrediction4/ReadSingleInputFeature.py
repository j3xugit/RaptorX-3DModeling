import os
import sys

import numpy as np
import cPickle

from ReadProteinFeatures import ReadFeatures

def Usage():
	print 'python ReadSingleInputFeature.py proteinName featureFolder [savefolder]'
    	print '  This script reads one set of input feature for a single protein and save the result as dict() in a PKL file.'
    	print '  Only input features are loaded. To load the ground truth matrices, please use code in Util/.'
	print '	 The resultant file has name proteinName.inputFeatures.pkl and will be saved to savefolder'
	print '  savefolder is optional and set to current work directory if not specified'

def main(argv):
    	if len(argv) < 2:
        	Usage()
        	exit(1)

    	proteinName = argv[0]
    	featureDir = argv[1]
	savefolder ='./'
	if len(argv) >= 3:
		savefolder = argv[2]

    	d = featureDir
	if not os.path.isdir(d):
		print 'The feature directory does not exist: ', d
		exit(1)

	pFeature = ReadFeatures( p=proteinName, DataSourceDir=d)
	if pFeature is None:
		print 'Fatal ERROR: failed to read protein feature from ', d
		exit(1)

    	savefile = proteinName + '.inputFeatures.pkl'
	savefile = os.path.join(savefolder, savefile)
    	with open(savefile, 'wb') as savefh:
    		cPickle.dump( pFeature, savefh,  protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
	main(sys.argv[1:])
