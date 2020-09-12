import numpy as np
import cPickle

import sys
import os

sys.path.append(os.path.join(os.environ['ModelingHome'], 'Common') )

from ReadProteinFeatures import ReadFeatures

def Usage():
    print 'python ReadOneProteinFeatures.py proteinName featureFolder1 featureFolder2 featureFolder3 ... '
    print '  This script reads several sets of independent features for a single protein and save the results as a list of features in a PKL file.'
    print '  Each feature folder contains a set of features.'
    print '  Only input features for distance prediciton are loaded. To load the ground truth distance matrix, please use code in Util/.'


def main(argv):
    if len(argv) < 2:
        Usage()
        exit(1)

    proteinName = argv[0]
    featureDirs = argv[1:]

    #print proteinName
    #print featureDirs

    pFeatures = []
    for d in featureDirs:
	if not os.path.isdir(d):
		print 'The feature directory does not exist: ', d
		exit(1)

	pFeature = ReadFeatures( p=proteinName, DataSourceDir=d)
	if pFeature is None:
		print 'Fatal ERROR: failed to read protein features from ', d
		exit(1)

	#print pFeature.keys()
        pFeatures.append(pFeature)

    savefile = proteinName + '.distanceFeatures.pkl'
    print('In total %d feature sets loaded and they are saved to %s ' % (len(pFeatures), savefile ) )
    savefh = open(savefile, 'wb')
    cPickle.dump( pFeatures, savefh,  protocol=cPickle.HIGHEST_PROTOCOL)
    savefh.close()


if __name__ == "__main__":
    ## read several sets of protein features into a single PKL file
    main(sys.argv[1:])
