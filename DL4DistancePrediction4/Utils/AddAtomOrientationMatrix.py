import numpy as np
import sys
import os
import cPickle


"""
This script merges distance prediction input features with atom orientation matrices to form a new train/valid/test dataset for orientation prediction
Usage: python AddAtomOrientationMatrix.py distanceFeatures_PKL atomOrientationMatrixPKL_folder
"""
if len(sys.argv)<3:
	print 'python AddAtomOrientationMatrix.py distanceFeatures_PKL atomOrientationMatrixPKL_folder'
	print '       distanceFeautres_PKL is a PKL file containing a list of input features for distance/orientation prediction'
	print '       atomOrientationMatrixPKL_folder is a folder containing the true orientation matrices of all the proteins in distanceFeatures_PKL'
	exit(1)

originalFeatureFile = sys.argv[1]
atomOrientationMatrixDir = sys.argv[2]


##open the original input feature file
ofh = open(originalFeatureFile, 'rb')
orifeatures = cPickle.load(ofh)
ofh.close()

print 'loaded original features for ', len(orifeatures), ' proteins'


mergedFeaturesPool = []
missedProteins = []
##merge the features
for protein in orifeatures:
	name = protein['name']
	distfile = os.path.join(atomOrientationMatrixDir, name + '.atomOrientationMatrix.pkl')
	if not os.path.isfile(distfile):
		print 'WARNING: there is no atom-level orientation matrix file for ', name
		missedProteins.append(name)
		continue

	##load the orientation matrix
	distfh = open(distfile, 'rb')
	atomOrientationMatrices = cPickle.load(distfh)
	distfh.close()

	if atomOrientationMatrices['seq4matrix'] != protein['sequence']:
		print 'inconsistent sequence between the orientation matrix file and the input feature file for ', name
		print 'atomOrientationMatrix sequence: ', atomOrientationMatrices['seq4matrix']
		print 'original       sequence: ', protein['sequence']
		missedProteins.append(name)
		continue

	protein['atomOrientationMatrix'] = atomOrientationMatrices

	mergedFeaturesPool.append(protein)

print 'in total merged features for ', len(mergedFeaturesPool), ' proteins'
print '# unmerged proteins: ', len(missedProteins)

##save the merged features
filenameItems = os.path.basename(originalFeatureFile).split('.')[0: -2] + ['distOriFeatures', 'pkl']
filename = '.'.join( filenameItems)

print 'writing the merged features to file ', filename
fh = open(filename, 'wb')
cPickle.dump(mergedFeaturesPool, fh, protocol = cPickle.HIGHEST_PROTOCOL)
fh.close()
