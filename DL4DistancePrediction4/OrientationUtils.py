import os
import numpy as np
import cPickle

import config
import RangeNWeight

## please double check to make sure that these two constants are consistent with those in CalcOrientationMatrixFromSeqPDB.py
##NoCB = np.float16(540)
##NoValidCoordinates = np.float16(720)

from Common.PDBUtils import NoCBDegree, InvalidDegree

## load native distance matrix from a PKL file  
def LoadNativeOriMatrixFromFile(filename):
        if not os.path.isfile(filename):
                print 'ERROR: the native orientation matrix file does not exist: ', filename
                exit(1)

        fh = open(filename, 'rb')
        oriMatrix = cPickle.load(fh)
        fh.close()
        return oriMatrix

## load native distance matrix by protein name and location
def LoadNativeOriMatrix(name, location='AtomOrientationMatrix/'):
        filename = os.path.join(location, name+'.atomOrientationMatrix.pkl')
        return LoadNativeOriMatrixFromFile(filename)

## add the native orientation matrix to a protein data
def AddOrientationMatrix(protein, atomOriMatrix):
        if atomOriMatrix['seq4matrix'] != protein['sequence']:
                print 'ERROR: inconsistent sequence between the orientation matrix file and the input feature file for ', protein['name']
                print 'atomOriMatrix sequence: ', atomOriMatrix['seq4matrix']
                print 'original       sequence: ', protein['sequence']
                exit(1)

        protein['atomOrientationMatrix'] = atomOriMatrix


## treat all invalid orientation angles as an extra label
## treat NoCBDegree  differently depending on if this matrix is used as response or input features derived from templates?
def DiscretizeOrientationMatrix(orim, bins=None, distThreshold4Orientation=20., distMatrix=None, invalidEntrySeparated=False, asResponse=True):

        assert bins is not None
	assert distMatrix is not None
	assert orim.shape == distMatrix.shape

	## when the orientation angle is out of bounds defined by bins, it may be digitized into 0 or len(bins)
        result = np.digitize(orim, bins)

        ## merge the label 0 into label len(bins)
        np.putmask(result, result==0, len(bins) )
	
	## change the labels so that they start from 0
	result -= 1
        labels = range( len(bins)  )

	## do some postprocessing here

	if asResponse:
		#1. when one Cb is missing (i.e., the raw dihedral/angle value=+ or -540), set a random label other than the last one
		np.putmask(result, abs(orim)==abs(NoCBDegree), np.random.choice(labels[:-1]) )

	#2. when two residues are far away from each other, set the label to the largest
	np.putmask(result, distMatrix > distThreshold4Orientation, labels[-1])

	#3. When one residue has no valid 3D coordinates, set the label to the largest
	np.putmask(result, distMatrix < 0, labels[-1] )

	#3. When Plus or Minus is used, separate the invalid entry from the largest bin, but we do not recommend using Plus/Minus for orientation 
	if invalidEntrySeparated:
		#np.putmask(result, orim < -180, len(bins) )
		#np.putmask(result, orim > 180, len(bins) )
		np.putmask(result, abs(orim) == abs(InvalidDegree), len(bins) )
		labels.append(len(bins))
	
        return result.astype(np.int16), np.int16(labels), bins

def CalcLabelProb(data=None, numLabels=37, eps=np.finfo(np.float32).eps, numRanges=-1):

        freqs = [ ]
        #for separation in config.RangeBoundaries:
        for separation in RangeNWeight.GetRangeBoundaries(numRanges):
                #print 'separation=', separation
                freq = []
                for m in data:
                        index = np.triu_indices(m.shape[0], separation)
                        values = m[index]
                        res = np.bincount(values, minlength=numLabels )
                        freq.append(res)
                freqs.append(np.sum(freq, axis=0) )

        count = np.array(freqs)
        #print count.shape
        #print count

        ## the order of subtraction cannot be changed
        ## count[0], [1], [2], [3], [4] are for extra long-, long-, medium-, short- and near-range residue pairs, respectively
        for i in range(count.shape[0]-1, 0, -1):
                count[i] -= count[i-1]

        frequency = [ c/(eps + np.sum(c) ) for c in count ]
        return np.array(frequency)

def CalcOrientationProb(data=None, bins=None, distMatrices=None, invalidEntrySeparated=False, numRanges=-1):

	assert bins is not None
	assert distMatrices is not None

	labelMatrices = []
	for orim, dm in zip(data, distMatrices):
		m, labels, _ = DiscretizeOrientationMatrix(orim, bins=bins, distMatrix=dm, invalidEntrySeparated=invalidEntrySeparated)
		labelMatrices.append(m)

	probs = CalcLabelProb(labelMatrices, len(labels), numRanges=numRanges )

	return probs

def DeriveOriContactMatrix(predOriMatrix, response):
	labelName, labelType, subType = config.ParseResponse(response)
	symmetric = config.IsSymmetricLabel(labelName)	

	if labelName not in config.allOrientationNames:
		print 'ERROR: unsupported orientation label name in', response
		exit(1)

        if not config.IsDiscreteLabel(labelType):
                print 'ERROR: unsupported orientation label type in', response
                exit(1)

        numLabels = config.GetResponseProbDims(response)
        if subType.endswith('Plus') or subType.endswith('Minus'):
                largestValidLabel = numLabels-2
        else:
                largestValidLabel = numLabels-1

	contactMatrix = np.sum( predOriMatrix[:, :, :largestValidLabel], axis=2)

	return contactMatrix
		

## predOriMatrices is a python dict, each item is an predicted orientation matrix
def DeriveContactMatrices(predOriMatrices):

	contactMatrices = dict()
	for response, predOriMatrix in predOriMatrices.iteritems():
		contactMatrices[response] = DeriveOriContactMatrix(predOriMatrix, response)
	
	return contactMatrices
