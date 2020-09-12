import os
import sys
import numpy as np
import cPickle

from Common.LoadHHM import load_hhm as LoadHHM
from DL4PropertyPrediction import PropertyUtils

ccmpredZSuffix=".ccmpred_zscore"
ccmpredSuffix=".ccmpred"
psicovSuffix=".psicov_zscore"
otherPairFeatureSuffix=".pot"

def Usage():
	print 'python ReadSingleInputFeature2.py proteinName featureFolder [savefolder]'
    	print '	This script generates input feature for contact/distance/orientation prediction for one protein and save the result as dict() in a PKL file.'
	print '	featurefolder: it shall contain the following files: .hhm, .predictedProperties.pkl, .ccmpred, .ccmpred_zscore, and .pot '
	print '	savefolder: optional, if not specified set to current work directory'
	print '	The resultant file has name proteinName.inputFeatures.pkl and will be saved to savefolder'

def LoadECMatrix(file, seqName=None, seq=None):
	fh=open(file, 'r')
    	allECs = []
    	for line in list(fh):
        	ECs = [ np.float16(x) for x in line.split() ]
        	if seq is not None:
            		assert len(seq) == len(ECs)
        	allECs.append(ECs)
    	fh.close()

    	if seq is not None:
        	assert len(seq) == len(allECs)

    	ECMatrix = np.array(allECs)
    	if np.isnan( np.sum(ECMatrix.astype(np.float32) ) ):
        	print 'ERROR: there is at least one NaN in ', file
        	exit(1)

    	return ECMatrix

def LoadOtherPairFeatures(file, seqName=None, seq=None):
    	if seq is None:
        	print 'ERROR: Please provide the sequence content for: ', seqName
        	exit(1)

    	indexList = []
    	valueList = []
    	fh = open(file, 'r')
    	for line in list(fh):
        	fields = line.split()
        	indices = [ np.int16(x)-1 for x in fields[0:2] ]
        	values = [ np.float16(x) for x in fields[2:] ]

        	indexList.append(indices)
        	valueList.append(values)

    	fh.close()

    	indexArr = np.transpose(np.array(indexList) )
    	if np.amin(indexArr) < 0 or np.amax(indexArr) >= len(seq):
        	print 'ERROR: In LoadOtherPairFeatures: index out of seq length, protein Name: ', seqName, file
        	exit(1)

    	allPairs = np.zeros((len(seq), len(seq), len(valueList[0]) ), dtype=np.float16 )
    	allPairs[ indexArr[0], indexArr[1] ] = valueList

    	##add the below statement to make the matrix symmetric
    	allPairs[ indexArr[1], indexArr[0] ] = valueList

    	if np.isnan( np.sum(allPairs.astype(np.float32) ) ):
        	print 'ERROR: there are NaNs in file ', file
        	exit(1)

    	return allPairs

def ReadFeatures(p=None, DataSourceDir=None):
	if p is None:
        	print 'ERROR: please specify a valid target name!'
        	exit(1)
    	if DataSourceDir is None:
        	print 'ERROR: please specify a folder containing all the features for the target!'
        	exit(1)
    	if not os.path.isdir(DataSourceDir):
        	print 'ERROR: the feature directory does not exist: ', DataSourceDir
        	return None

    	OneProtein=dict()
    	OneProtein['name'] = p

    	seqFile=os.path.join(DataSourceDir, p + ".seq")
    	if not os.path.isfile(seqFile):
        	print 'ERROR: the sequence file does not exist: ', seqFile
        	return None
    	OneProtein['sequence'] = SequenceUtils.LoadFASTAFile(seqFile)

	propertyFile = os.path.join(DataSourceDir, p + ".predictedProperties.pkl")
	name, sequence, predProb, predString = PropertyUtils.LoadPredictedProperties(propertyfile)

	if np.isnan(np.sum(predProb['SS8_Discrete8C']) ) or np.isnan(np.sum(predProb['SS3_Discrete3C']) ) or np.isnan(np.sum(predProb['ACC_Discrete3C']) ):
		print 'ERROR: there are NaNs in the property file: ', propertyFile
		return None

    	OneProtein['SS3'] = predProb['SS3_Discrete3C']
	OneProtein['SS8'] = predProb['SS8_Discrete8C']
    	OneProtein['ACC'] = predProb['ACC_Discrete3C']

	hhmfile = os.path.join(DataSourceDir, p + ".hhm")
    	hhm = LoadHHM(hhmfile)
    	if np.isnan( np.sum(hhm['PSFM']) ) or np.isnan(np.sum(hhm['PSSM']) ):
        	print 'ERROR: There are NaNs in the hhm file: ', hhmfile
		return None

    	OneProtein['PSFM'] = hhm['PSFM']
    	OneProtein['PSSM'] = hhm['PSSM']

	ccmpredZf = os.path.join(DataSourceDir, p + ccmpredZSuffix)
    	if os.path.isfile(ccmpredZf):
        	OneProtein['ccmpredZ'] = LoadECMatrix(ccmpredZf, seqName=p, seq=OneProtein['sequence'])

    	ccmpredf = os.path.join(DataSourceDir, p + ccmpredSuffix)
    	OneProtein['ccmpred'] = LoadECMatrix(ccmpredf, seqName=p, seq=OneProtein['sequence'])

    	psicovf = os.path.join(DataSourceDir, p + psicovSuffix)
    	if os.path.isfile(psicovf):
        	OneProtein['psicovZ'] = LoadECMatrix(psicovf, seqName=p, seq=OneProtein['sequence'])

   	 ## OtherPairs is a L*L*3 matrix consisting of contact potential, MI and APC-corrected MI
    	otherPairFeaturesf = os.path.join(DataSourceDir, p + otherPairFeatureSuffix)
    	OneProtein['OtherPairs'] = LoadOtherPairFeatures(otherPairFeaturesf, seqName=p, seq=OneProtein['sequence'])

	return OneProtein

def main(argv):
    	if len(argv) < 2:
        	Usage()
        	exit(1)

    	proteinName = argv[0]
    	featureDir = argv[1]
	savefolder =os.getcwd()
	if len(argv) >= 3:
		savefolder = argv[2]
	if not os.path.isdir(savefolder):
		os.mkdir(savefolder)

    	d = featureDir
	if not os.path.isdir(d):
		print 'ERROR: the feature directory does not exist: ', d
		exit(1)

	pFeature = ReadFeatures( p=proteinName, DataSourceDir=d)
	if pFeature is None:
		print 'Fatal ERROR: failed to read protein feature for ', proteinName,' from ', d
		exit(1)

    	savefile = proteinName + '.inputFeatures.pkl'
	savefile = os.path.join(savefolder, savefile)
    	with open(savefile, 'wb') as savefh:
    		cPickle.dump( pFeature, savefh,  protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
	main(sys.argv[1:])
