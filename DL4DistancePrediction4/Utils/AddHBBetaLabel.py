import numpy as np
import sys
import os
import cPickle


"""
This script merges the distance prediction features with the HBBeta matrices to form a new training and validation set for atom-level distance, hydrogen-bonding, beta-pairing prediction
Usage: python AddHBBetaFeatures.py originalDistanceFeaturesPKL HBBetaMatrixPKL_folder
"""
if len(sys.argv)<3:
	print 'python AddHBBetaLabel.py originalDistanceFeaturesPKL HBBetaMatrixPKL_folder'
	exit(-1)

originalFeatureFile = sys.argv[1]
HBBetaMatrixDir = sys.argv[2]


##open the original distance feature file
ofh = open(originalFeatureFile, 'rb')
orifeatures = cPickle.load(ofh)
ofh.close()

print 'loaded original features for ', len(orifeatures), ' proteins'


mergedFeaturesPool = []
missedProteins = []
##merge the features
for protein in orifeatures:
	name = protein['name']
	hbfile = os.path.join(HBBetaMatrixDir, name + '.HBBeta.pkl')
	if not os.path.isfile(hbfile):
		print 'WARNING: there is no HBBeta matrix file for ', name
		missedProteins.append(name)
		continue

	##load the dist matrix
	hbfh = open(hbfile, 'rb')
	HBBetaMatrices = cPickle.load(hbfh)
	hbfh.close()

	##verify the length consistency
	if HBBetaMatrices['HB'].shape[0] != len(protein['sequence']):
		print 'inconsistent protein lengths for ', name
		exit(-1)
	if HBBetaMatrices['sequence'] != protein['sequence']:
		print 'inconsistent sequence between the HBBeta matrix file and the original distance feature file for ', name
		exit(-1)

	## in the HBBeta.pkl file, 1 for positive and 0 for negative
	## but in our distance prediction model, we use 0 for positive and 1 for negative
	if not protein.has_key('atomDistMatrix'):
		print 'protein ', name, ' has no atomDistMatrix in the keys. You may first run AddAtomDistance.py '
		exit(-1)

	protein['atomDistMatrix']['HB'] = 1 - HBBetaMatrices['HB'].toarray()
	protein['atomDistMatrix']['Beta'] = 1 - HBBetaMatrices['BetaPairing'].toarray()

	mergedFeaturesPool.append(protein)

print 'in total merged features for ', len(mergedFeaturesPool), ' proteins'
print '# unmerged proteins: ', len(missedProteins)

##save the merged features
filenameItems = os.path.basename(originalFeatureFile).split('.')[0: -2] + ['DistHBBetaFeatures', 'pkl']
filename = '.'.join( filenameItems)

print 'writing the merged features to file ', filename
fh = open(filename, 'wb')
cPickle.dump(mergedFeaturesPool, fh, protocol = cPickle.HIGHEST_PROTOCOL)
fh.close()
