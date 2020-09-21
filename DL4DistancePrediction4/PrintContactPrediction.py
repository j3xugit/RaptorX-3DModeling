import os
import sys
import cPickle
import numpy as np
import getopt

import config
import ContactUtils
import DistanceUtils

## find all strings in strList that start with prefix
def FindStringsStartWith(strList, prefix):
	results = []
	for str in strList:
		if str.startswith(prefix):
			results.append(str)
	return results

def Usage():
	print 'python PrintContactPrediction.py [-a | -c | -d savefolder ] predictedDistMatrixPKL'
	print '	This script prints predicted contacts from a predicted distance/orientation file and saves the result in two text formats: matrix and CASP'
	print '	-a: if is specified, also prints contact information for other atom pairs, e.g., CaCa, defaut no'
	print '	  	for other atom pairs, the result files are named after targetName.XX.CASP.rr and targetName.XX.CM.txt where XX represents the atom pair'
	print '	  	for Cb-Cb pairs, the result files are named after targetName.CASP.rr and targetName.CM.txt'
	print '	-c: print contact probability only, but not distance. By default, both contact and distance probability will be printed'
	print '	-d: folder for result saving, default current work directory'

if len(sys.argv) < 2:
	Usage()
	exit(1)

savefolder = os.getcwd()
bPrintOtherAtomPairs = False
contactOnly = False

try:
	opts, args = getopt.getopt(sys.argv[1:], "d:ca", ["savefolder=", "contactOnly=", "allatompairs="])
except getopt.GetoptError as err:
	print err
	Usage()
	exit(1)

for opt, arg in opts:
	if opt in ("-d", "--savefolder"):
		savefolder = arg

	elif opt in ("-c", "--contactOnly"):
		contactOnly = True
		
	elif opt in ("-a", "--allatompairs"):
		bPrintOtherAtomPairs = true
	else:
		Usage()
		exit(1)

if len(args) < 1:
	Usage()
	exit(1)

if not os.path.isdir(savefolder):
	os.mkdir(savefolder)

predFile = args[0]

with open(predFile, 'rb') as fh:
	pred = cPickle.load(fh)
targetName, sequence, distProbMatrix, contactMatrices = pred[:4]

if not contactMatrices.has_key('CbCb'):
	print 'ERROR: the predicted matrix file does not have information for Cb-Cb atom pairs'
	exit(1)

filename = os.path.basename(predFile).split('.')[0]

for apt, m in contactMatrices.iteritems():
	if apt not in config.allDistLabelNames:
		continue

	if apt == 'CbCb':
        	contactFileName = filename + '.CM.txt'
                contactCASPFileName = filename + '.CASP.rr'

        elif bPrintOtherAtomPairs:
                contactFileName = filename + '.' + apt + '.CM.txt'
                contactCASPFileName = filename + '.' + apt + '.CASP.rr'
	else:
		continue

	contactFile = os.path.join(savefolder, contactFileName)
        np.savetxt(contactFile, m, fmt='%1.6f', delimiter=' ')

	contactCASPFile = os.path.join(savefolder, contactCASPFileName)
	if contactOnly:
        	ContactUtils.SaveContactMatrixInCASPFormat(targetName, sequence, m, contactCASPFile, distMatrix=None, probScaleFactor=1)
		continue

	responses = FindStringsStartWith(distProbMatrix.keys(), apt)
	if len(responses) != 1:
		## right now for one apt, only one response is allowed
		print 'ERROR: incorrect distance information for', apt, 'in', predFile
		exit(1)

	response = responses[0]
	labelName, labelType, subType = config.ParseResponse(response)

	if not config.IsDiscreteLabel(labelType):
		print 'ERROR: right now only discrete distance matrix is supported'
		exit(1)

	## convert distance matrix to what's needed by CASP
	distMatrix = DistanceUtils.MergeDistanceBinsBySum(distProbMatrix[response], config.distCutoffs[subType], config.distCutoffs['10C'])
        ContactUtils.SaveContactMatrixInCASPFormat(targetName, sequence, m, contactCASPFile, distMatrix=distMatrix, probScaleFactor=1)
