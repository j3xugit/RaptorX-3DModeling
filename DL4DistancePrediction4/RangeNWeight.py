import numpy as np

def GetNumRanges(modelSpecs=None):
	if modelSpecs is None or (not modelSpecs.has_key('numRanges')):
		return 4
	return modelSpecs['numRanges']

def GetRangeBoundaries(numRanges=-1):
	## sequence separation of two residues for different ranges

	## five ranges: extreme-long, long, medium, short and near-ranges
	FiveRangeBoundaries = [48, 24, 12, 6, 0]

	## 4 ranges: long, medium, short and near-ranges
	FourRangeBoundaries = [24, 12, 6, 0]

	## only one range is used
	OneRangeBoundaries = [0]

	if numRanges == 1:
		return OneRangeBoundaries
	elif numRanges == 4:
		return FourRangeBoundaries
	elif numRanges == 5:
		return FiveRangeBoundaries
	else:
		print 'ERROR: unsupported numRanges: ', numRanges
		exit(1)

## this function determines the range into which one pair of residues fall 
## offset = abs(i-j) where i and j are the two residue positions
def GetRangeIndex(offset, numRanges=-1):
	RangeBoundaries = GetRangeBoundaries(numRanges)
        if offset < RangeBoundaries[-1]:
                return -1

        rangeIndex = 0
        for l in range(numRanges):
                if offset >= RangeBoundaries[l]:
                        break
                else:
                        rangeIndex += 1
        return rangeIndex

## set weight to different ranges
def SetWeight4Range(modelSpecs=None):
	if modelSpecs is None:
		print 'ERROR: modelSpecs is None for SetWeight4Range'
		exit(1)

	numRanges = GetNumRanges(modelSpecs)
	if modelSpecs.has_key('NoWeight4Range') and modelSpecs['NoWeight4Range']:
	#if config.NoWeight4Range(modelSpecs):
		weight4range = np.array( [1. ] * numRanges ).reshape( (-1, 1) ).astype(np.float32)
	elif numRanges == 5:
		weight4range = np.array([3.5, 3., 2.5, 1., 0.5]).reshape((-1,1)).astype(np.float32)
	elif numRanges == 4: 
		weight4range = np.array([3., 2.5, 1., 0.5]).reshape((-1,1)).astype(np.float32)
	elif numRanges == 1:
		weight4range = np.array([1.]).reshape((-1,1)).astype(np.float32)
	else:
		print 'ERROR: unsupported numRanges: ', numRanges
		exit(1)

	modelSpecs['weight4range'] = weight4range

	return weight4range


## set weight for 2C and/or 3C labels where 2C and 3C indicates 2 labels and 3 labels, respectively
## the weight matrix is 2D and has dimension numRanges * 2 or numRanges*3. Each row is for one specific range.
##for example, in [17, 4, 1], 17 is the weight for distance bin 0-8, 4 for distance bin 8-15 and 1 for distance >15 or -1
def SetWeight43C2C(modelSpecs=None):
	if modelSpecs is None:
		print 'ERROR: modelSpecs is None for SetWeight43C2C'
		exit(1)

	if not modelSpecs.has_key('weight4range'):
		print 'ERROR: please set weight4range in modelSpecs first'
		exit(1)

	numRanges = GetNumRanges(modelSpecs)
	weight43C = dict()

	if numRanges == 4:
		weight43C['low'] = np.array( [ [17, 4, 1], [5, 2, 1], [2.5, 0.6, 1], [0.2, 0.3, 1] ] ).astype(np.float32)
		weight43C['mid'] = np.array( [ [20.5, 5.4, 1], [5.4, 1.89, 1], [2.9, 0.7, 1], [0.2, 0.3, 1] ] ).astype(np.float32)
		weight43C['high']= np.array( [ [23, 6, 1], [6, 2.5 ,1], [3, 1, 1] ,[0.2, 0.3, 1] ] ).astype(np.float32)
		weight43C['veryhigh'] = np.array( [ [25, 6, 1], [7, 2.5 ,1], [3, 1, 1], [0.2, 0.3, 1] ] ).astype(np.float32)
		weight43C['exhigh']  =np.array( [ [28, 6, 1], [8, 2.5 ,1], [4, 1, 1], [0.2, 0.3, 1] ] ).astype(np.float32)

		weight4Beta2C = np.array( [ [360, 1], [70, 1], [50, 1], [120, 1] ] ).astype(np.float32) 
		weight4HB2C = np.array( [ [600., 1], [120., 1], [90., 1], [5., 1] ] ).astype(np.float32) 

	elif numRanges == 5:
		weight43C['low'] = np.array( [[25, 6, 1],  [17, 4, 1], [5, 2, 1], [2.5, 0.6, 1], [0.2, 0.3, 1] ] ).astype(np.float32)
		weight43C['mid'] = np.array( [ [30, 8.1, 1], [20.5, 5.4, 1], [5.4, 1.89, 1], [2.9, 0.7, 1], [0.2, 0.3, 1] ] ).astype(np.float32)
		weight43C['high']= np.array( [ [34, 9, 1], [23, 6, 1], [6, 2.5 ,1], [3, 1, 1] ,[0.2, 0.3, 1] ] ).astype(np.float32)
		weight43C['veryhigh'] = np.array( [ [37.5, 9, 1], [25, 6, 1], [7, 2.5 ,1], [3, 1, 1], [0.2, 0.3, 1] ] ).astype(np.float32)
		weight43C['exhigh']  =np.array( [ [42, 9, 1],  [28, 6, 1], [8, 2.5 ,1], [4, 1, 1], [0.2, 0.3, 1] ] ).astype(np.float32)

		# weight for Beta-pairing, only two labels, 0 for positive and 1 for negative
		weight4Beta2C = np.array( [ [420, 1], [360, 1], [70, 1], [50, 1], [120, 1] ] ).astype(np.float32)
		# weight for hydrogen-bonding, only two labels, 0 for positive and 1 for negative
		weight4HB2C = np.array( [ [700, 1], [600., 1], [120., 1], [90., 1], [5., 1] ] ).astype(np.float32)

	elif numRanges == 1:
		weight43C['low'] = np.array( [ [17, 4, 1] ] ).astype(np.float32)
		weight43C['mid'] = np.array( [ [20.5, 5.4, 1] ] ).astype(np.float32)
		weight43C['high']= np.array( [ [23, 6, 1] ] ).astype(np.float32)
		weight43C['veryhigh'] = np.array( [ [25, 6, 1] ] ).astype(np.float32)
		weight43C['exhigh']  =np.array( [ [28, 6, 1] ] ).astype(np.float32)

		weight4Beta2C = np.array( [ [360, 1] ] ).astype(np.float32) 
		weight4HB2C = np.array( [ [600., 1] ] ).astype(np.float32) 

	else:
		print 'ERROR: unsupported numRanges: ', numRanges

	if modelSpecs.has_key('LRbias'):
                modelSpecs['weight4Discrete3C']= np.multiply(weight43C[modelSpecs['LRbias'] ], modelSpecs['weight4range'])
        else:
                modelSpecs['weight4Discrete3C']= np.multiply(weight43C['mid'], modelSpecs['weight4range'])

	modelSpecs['weight4HB_Discrete2C'] = np.multiply(weight4HB2C, modelSpecs['weight4range'])
        modelSpecs['weight4Beta_Discrete2C'] = np.multiply(weight4Beta2C, modelSpecs['weight4range'])

        ## weight for real value
        modelSpecs['weight4continuous'] = np.multiply(np.array([1.] * numRanges).reshape((-1, 1)).astype(np.float32), modelSpecs['weight4range'])

	return modelSpecs

