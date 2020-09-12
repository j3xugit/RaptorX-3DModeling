import numpy as np
import sys
import os
import cPickle


"""
This script merges distance prediction input features with atom distance matrices to form a new train/valid/test dataset for distance prediction
Usage: python AddAtomDistMatrix.py distanceFeatures_PKL atomDistMatrixPKL_folder
"""
if len(sys.argv)<3:
	print 'python AddAtomDistMatrix.py distanceFeatures_PKL atomDistMatrixPKL_folder'
	print '       distanceFeautres_PKL is a PKL file containing a list of input features for distance prediction'
	print '       atomDistMatrixPKL_folder is a folder containing the true distance matrices of all the proteins in distanceFeatures_PKL'
	exit(1)

originalFeatureFile = sys.argv[1]
atomDistMatrixDir = sys.argv[2]


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
	distfile = os.path.join(atomDistMatrixDir, name + '.atomDistMatrix.pkl')
	if not os.path.isfile(distfile):
		print 'WARNING: there is no atom-level distane matrix file for ', name
		missedProteins.append(name)
		continue

	##load the dist matrix
	distfh = open(distfile, 'rb')
	atomDistMatrices = cPickle.load(distfh)
	distfh.close()

	##verify the legnth consistency
	if atomDistMatrices['CbCb'].shape[0] != len(protein['sequence']):
		print 'inconsistent protein lengths for ', name
		missedProteins.append(name)
		continue

	if atomDistMatrices['seq4matrix'] != protein['sequence']:
		print 'inconsistent sequence between the distance matrix file and the input feature file for ', name
		print 'atomDistMatrix sequence: ', atomDistMatrices['seq4matrix']
		print 'original       sequence: ', protein['sequence']
		missedProteins.append(name)
		continue

	if protein.has_key('DistMatrix'):
		del protein['DistMatrix']
	protein['atomDistMatrix'] = atomDistMatrices

	mergedFeaturesPool.append(protein)

print 'in total merged features for ', len(mergedFeaturesPool), ' proteins'
print '# unmerged proteins: ', len(missedProteins)

##save the merged features
filenameItems = os.path.basename(originalFeatureFile).split('.')[0: -2] + ['distanceFeatures', 'pkl']
filename = '.'.join( filenameItems)

print 'writing the merged features to file ', filename
fh = open(filename, 'wb')
cPickle.dump(mergedFeaturesPool, fh, protocol = cPickle.HIGHEST_PROTOCOL)
fh.close()
