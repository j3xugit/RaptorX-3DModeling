import sys
import os
import numpy as np
import cPickle

from Common import SequenceUtils
from Common import LoadTPLTGT

ccmpredZSuffix=".ccmpred_zscore"
ccmpredSuffix=".ccmpred"
psicovSuffix=".psicov_zscore"
otherPairFeatureSuffix=".pot"

##in the future we may add some code to check correctness
def LoadSS3(file, seqName=None, seq=None):
    fh=open(file, 'r')
    allprobs = []
    AAs = ""

    for line in list(fh)[3:]:
        probs = [ np.float32(x) for x in line.split()[3:] ]
        allprobs.append(probs)

        AAs += line.split()[1]

    fh.close()

    if seq is not None:
        assert len(seq) == len(allprobs)
        assert seq == AAs

    SS3probs = np.array(allprobs)

    if np.isnan( np.sum(SS3probs) ):
	print 'ERROR: there are NaNs in file ', file
	exit(1)

    return SS3probs

def LoadACC(file, seqName=None, seq=None):
    fh=open(file, 'r')
    allprobs = []
    AAs = ""
    for line in list(fh)[5:]:
        probs = [ np.float32(x) for x in line.split()[3:6] ]
        allprobs.append(probs)

        AAs += line.split()[1]

    fh.close()

    if seq is not None:
        assert len(seq) == len(allprobs)
        assert seq == AAs

    ACCprobs = np.array(allprobs)

    if np.isnan( np.sum(ACCprobs) ):
	print 'ERROR: there are NaNs in file ', file
	exit(1)

    return ACCprobs

def LoadDISO(file, seqName=None, seq=None):
    fh=open(file, 'r')
    allprobs = []
    AAs = ""
    for line in list(fh)[4:]:
        probs = [ np.float32(x) for x in line.split()[3:] ]
        allprobs.append(probs)

        AAs += line.split()[1]

    fh.close()

    if seq is not None:
        assert len(seq) == len(allprobs)
        assert seq == AAs

    DISOprobs = np.array(allprobs)

    if np.isnan( np.sum(DISOprobs) ):
	print 'ERROR: there is NaN in file ', file
	exit(1)

    return DISOprobs

##this function will return both profile and predicted ss8
def LoadProfile(file, seqName=None, seq=None):
    fh=open(file, 'r')
    content = list(fh)
    fh.close()

    seqLen = int(content[1].strip())

    assert len(content) == (3 + 3*seqLen)
  
    if seqName is not None:
        assert content[0].strip() == seqName

    if seq is not None:
        assert content[2].strip() == seq

    allprobs = []
    for line in content[3 : 3+seqLen]:
	probs = [ np.float32(x) for x in line.split(',') ]
        allprobs.append(probs)

    allscores = []   
    for line in content[3+seqLen : 3+2*seqLen ]:
	scores = [ np.float32(x) for x in line.split(',') ]
	allscores.append(scores)

    allss8s = []
    for line in content[3+2*seqLen : 3+3*seqLen ]:
	ss8s = [ np.float32(x) for x in line.split(',')[0:8] ]
	allss8s.append(ss8s)

    return np.array(allprobs), np.array(allscores), np.array(allss8s)

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
        print 'Please provide the sequence content of seq with name: ', seqName
	sys.exit(1)

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
        print 'In LoadOtherPairFeatures: index out of seq length'
	sys.exit(1)
 
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
    ###
    fastafh=open(seqFile, "r")
    fastacontent = [ line.strip() for line in list(fastafh) ]
    if not fastacontent[0].startswith('>'):
        OneProtein['sequence'] = ''.join(fastacontent)
    else:
        OneProtein['sequence'] = ''.join(fastacontent[1: ])
    fastafh.close()
    ###

    ssf = os.path.join(DataSourceDir, p + ".ss3")
    OneProtein['SS3'] = LoadSS3(ssf, seqName=p, seq=OneProtein['sequence'])

    accf = os.path.join(DataSourceDir, p + ".acc")
    OneProtein['ACC'] = LoadACC(accf, seqName=p, seq=OneProtein['sequence'])

    disof = os.path.join(DataSourceDir, p + ".diso")
    OneProtein['DISO'] = LoadDISO(disof, seqName=p, seq=OneProtein['sequence'])

    tgtf = os.path.join(DataSourceDir, p + ".tgt")

    tgt = LoadTPLTGT.load_tgt(tgtf)
    OneProtein['PSFM'] = tgt['PSFM']
    OneProtein['PSSM'] = tgt['PSSM']
    OneProtein['SS8'] = tgt['SS8'] 

    if np.isnan( np.sum(tgt['PSFM']) ) or np.isnan(np.sum(tgt['PSSM']) ) or np.isnan(np.sum(tgt['SS8']) ):
	print 'ERROR: There are NaNs in the tgt file: ', tgtf
	exit(1)

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

def Usage():
    print 'python ReadProteinFeatures.py proteinListFile featureMetaFolder '
    print '  proteinListFile: the file containing a list of proteins, each protein name in a line'
    print '  featureMetaFolder: specify a folder containing all the features, under which each protein has an independent feature folder named after feat_proteinName_contact'
    print '  This script reads only protein features for distance prediciton. To load the true distance matrix, please use AddAtomDistMatrix.py in Utils/ '


def main(argv):
    if len(argv) < 2:
        Usage()
        sys.exit(1)

    listFile = argv[0]
    featureMetaDir = argv[1]

    print 'listFile=', listFile
    print 'feature folder=', featureMetaDir

    if not os.path.isdir(featureMetaDir):
        print 'the provided feature folder is invalid: ', featureMetaDir
        sys.exit(1)
    if not os.path.isfile(listFile):
	print 'the provided protein list file is invalid: ', listFile
	sys.exit(1)

    lfh = open(listFile, 'r')
    proteins = [ p.strip() for p in list(lfh) ]
    lfh.close()

    pFeatures = []
    for p in proteins:
	thisFeatureDir = os.path.join(featureMetaDir, 'feat_' + p + '_contact/')
	pFeature = ReadFeatures( p=p, DataSourceDir=thisFeatureDir)
	if pFeature is None:
		print 'WARNING: skip protein ', p, ' due to failure in reading its features from ', thisFeatureDir
		continue
        pFeatures.append(pFeature)

	if len(pFeatures)%500 == 0:
		print 'finished loading features for ', len(pFeatures), ' proteins'

    savefile = os.path.basename(listFile) + '.' + featureMetaDir.rstrip('/').split('/')[-1] + '.' + str(os.getpid()) + '.distanceFeatures.pkl'
    print 'Writing the result to ', savefile
    savefh = open(savefile, 'wb')
    cPickle.dump( pFeatures, savefh,  protocol=cPickle.HIGHEST_PROTOCOL)
    savefh.close()


if __name__ == "__main__":
    ## read a list of protein features into a single PKL file
    main(sys.argv[1:])
