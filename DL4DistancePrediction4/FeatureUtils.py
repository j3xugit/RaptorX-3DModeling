import os
import sys
import cPickle
import numpy as np

import config
import DataProcessor
import CCMpredUtils
import MSAUtils

def LoadFeaturePKL(name, location='Feature4Train_2017_E001_PKL/', modelSpecs=None):

	## load up the basic input features
	filename = os.path.join(location, name + '.inputFeatures.pkl')
	if not os.path.isfile(filename):
		print 'ERROR: the input feature file does not exist: ', filename
		exit(1)
	with open(filename) as fh:
		feature = cPickle.load(fh)

	## check to see if loading up extra features	
	bUseCCMFnorm, bUseCCMsum, bUseCCMraw, bUseFullMI, bUseFullCov = config.ParseExtraCCMmode(modelSpecs)

	bUseExtraCCM = bUseCCMFnorm or bUseCCMsum or bUseCCMraw
	if bUseExtraCCM:
		extrafile = os.path.join(location, name + '.extraCCM.pkl')
		if not os.path.isfile(extrafile):
			print 'ERROR: the file for extra CCM information does not exist: ', extrafile
			exit(1)
		with open(extrafile) as fh:
			extra=cPickle.load(fh)

	if bUseCCMFnorm:
		feature['CCMFnorm'] = extra['Fnorm']
		feature['CCMFnormZ'] = extra['FnormZ']

	seqLen = len(feature['sequence'])
	if bUseCCMsum:
		if not extra.has_key('sumCCM'):
			print 'ERROR: CCM summary is requested, but the file does not have it: ', extrafile
			exit(1)
		feature['sumCCM'] = CCMpredUtils.ExpandMatrix(extra['sumCCM'], seqLen)

	if bUseCCMraw:
		if not extra.has_key('rawCCM'):
			print 'ERROR: CCM raw matrix is requested, but the file does not have it: ', extrafile
			exit(1)
		feature['rawCCM'] = CCMpredUtils.ExpandMatrix(extra['rawCCM'], seqLen)

	if bUseFullMI:
		alnfile = os.path.join(location, name + '.a2m')
		if not os.path.isfile(alnfile):
			print 'ERROR: the a2m file does not exist: ', alnfile
			exit(1)
		feature['fullMI'] = MSAUtils.CalcPairMatrixFromFile(alnfile, matrixType='mi')

	if bUseFullCov:
		covfile = os.path.join(location, name + '.cov.pkl')
		if not os.path.isfile(covfile):
			alnfile = os.path.join(location, name + '.a2m')
			if not os.path.isfile(alnfile):
				print 'ERROR: the a2m file does not exist:', alnfile
				exit(1)
			feature['fullCov'] = MSAUtils.CalcPairMatrixFromFile(alnfile)
		else:
			with open(covfile, 'rb') as fh:
				feature['fullCov'] = cPickle.load(fh)

	## check to see if we shall load up ESM information
	layers = config.ParseESMmode(modelSpecs)
	if layers is not None:
		esmfile = os.path.join(location, name + '.esm2.pkl')
		if not os.path.isfile(esmfile):
			print 'ERROR: the file for ESM information does not exist: ', esmfile
			exit(1)

		with open(esmfile, 'rb') as fh:
			esm = cPickle.load(fh)

		esmfeature = []
		for layer in layers:
			layer4key = layer % (esm['numModelLayers'] + 1 )
			if not esm.has_key(layer4key):
				print 'ERROR: attention weight for layer ', layer, ' requested but not available in ', esmfile
				exit(1)
			esmfeature.append(esm[layer4key])

		feature['ESM'] = np.concatenate(esmfeature, axis=2)

		#print 'ESM feature has shape', feature['ESM'].shape
	
	return feature


##d is a dictionary for a protein
def LocationFeature(d):
        ##add one specific location feature here, i.e., posFeature[i, j]=min(1, abs(i-j)/50.0 )
        seqLen = len(d['sequence'])
        posFeature = np.ones((seqLen, seqLen), dtype=config.MyFloat)
        separation_cutoff = 50
        end = min(separation_cutoff - 1, posFeature.shape[0])

        for offset in xrange(0, end):
                i = np.arange(0, posFeature.shape[0]-offset)
                j = i + offset
                posFeature[i, j] = offset/(1. * separation_cutoff)

        for offset in xrange(1, end):
                i = np.arange(offset, posFeature.shape[0])
                j = i - offset
                posFeature[i, j] = offset/(1. * separation_cutoff)

        return posFeature

def CubeRootFeature(d):
        ##the cube root of the difference between the two positions, the radius of protein is related to this feature
        seqLen = len(d['sequence'])
        feature = []
        for i in range(seqLen):
                dVector = [ abs(j-i) for j in range(seqLen) ]
                feature.append(dVector)
        posFeature = np.cbrt( np.array( feature ).astype(config.MyFloat) )
        #np.putmask(posFeature, posFeature > np.cbrt(300.), np.cbrt(300.) )
        #posFeature /= np.cbrt(100)
        return posFeature

def LogFeature(d):
        ##the logarithm of the difference between the two positions, the radius of protein is related to this feature
        seqLen = len(d['sequence'])
        feature = []
        for i in range(seqLen):
                dVector = [ abs(j-i)+1 for j in range(seqLen) ]
                feature.append(dVector)
        posFeature = np.log10( np.array( feature ).astype(config.MyFloat) )
        np.putmask(posFeature, posFeature > np.log10(300.), np.log10(300.) )

        return posFeature

def NewLocationFeature(d):
        seqLen = len(d['sequence'])
        mg = np.mgrid[0:seqLen, 0:seqLen]
        diff = abs(mg[0] - mg[1])
        f1 = (1./(1+diff)).astype(np.float16)
        f2 = (1./np.sqrt(1+diff)).astype(np.float16)

        return [f1, f2]

def Rg(d, normLen=100, AsSequentialFeature=True):
        seqLen = len(d['sequence'])

        ## estimate Rg and normalize it by the Rg of a protein with 100 AAs
        rg = 0.395*(seqLen/normLen)**(3./5)+7.257
        if AsSequentialFeature:
                ret = np.array([rg]*seqLen).astype(np.float16)
        else:
                ret = np.float16(rg)

        return ret

def CalcFeatureExpect4OneProtein(d):
	seqfeatures = np.mean(d['seqFeatures'].astype(np.float32), axis=0)
	seqfeatures_weight = d['seqFeatures'].shape[0]

	matrixfeatures = np.mean(d['matrixFeatures'].astype(np.float32), axis=(0, 1) ) 
	matrixfeatures_weight = np.prod(d['matrixFeatures'].shape[:2])

	if d.has_key('embedFeatures'):
		embedfeatures = np.mean(d['embedFeatures'].astype(np.float32), axis=0) 
		embedfeatures_weight = d['embedFeatures'].shape[0] 

		return (seqfeatures, seqfeatures_weight, matrixfeatures, matrixfeatures_weight, embedfeatures, embedfeatures_weight)

	return (seqfeatures, seqfeatures_weight, matrixfeatures, matrixfeatures_weight)


## trainData is a list of real data
def CalcExpectedValueOfFeatures(trainData, modelSpecs):
	seqfeatures = []
	seqweights = []
	
	matrixfeatures = []
	matrixweights = []

	embedfeatures = []
	embedweights = []

	for d in trainData:
		res = CalcFeatureExpect4OneProtein(d)
		seqfeatures.append(res[0])
		seqweights.append(res[1])
		matrixfeatures.append(res[2])
		matrixweights.append(res[3])
		
		if len(res) == 6:
			embedfeatures.append(res[4])
			embedweights.append(res[5])

	modelSpecs['seqFeatures_expected'] = np.average(seqfeatures, axis=0, weights=seqweights)
	modelSpecs['matrixFeatures_expected'] = np.average(matrixfeatures, axis=0, weights=matrixweights)
	modelSpecs['embedFeatures_expected'] = np.average(embedfeatures, axis=0, weights=embedweights)

## if seqFeatures is an ndarray, then return shape[1]
## if it is a list, then we need to sum up
def DetermineNumSeqFeatures(seqFeatures):
	if isinstance(seqFeatures, (list, tuple) ):
		numSeqFeatures = 0
		for f in seqFeatures:
			if len(f.shape) == 1:
				numSeqFeatures += 1
			elif len(f.shape) == 2:
				numSeqFeatures += f.shape[1]
			else:
				print 'wrong shape for input sequential feature: ', f.shape
				exit(1)
	else:
		numSeqFeatures = seqFeatures.shape[1]
	return numSeqFeatures

def DetermineNumMatrixFeatures(matrixFeatures):
	if isinstance(matrixFeatures, (list, tuple) ):
		numMatrixFeatures = 0
		for f in matrixFeatures:
			if len(f.shape) == 2:
				numMatrixFeatures += 1
			elif len(f.shape) == 3:
				numMatrixFeatures += f.shape[2]
			else:
				print 'wrong shape for input matrix feature: ', f.shape
				exit(1)
	else:
		numMatrixFeatures = matrixFeatures.shape[2]
	return numMatrixFeatures

## for training purpose, this function will sample input feature for one training protein and determine the dimension of its input features
def DetermineFeatureDimensionBySampling(metaData, modelSpecs):

	protein = DataProcessor.SampleProteinInfo(metaData, numSamples=1)[0]
        d = DataProcessor.LoadRealData(protein, modelSpecs, loadLabel=False, returnMode='list')

        ## obtain the dimension of each type of input feature
	modelSpecs['n_in_seq'] = DetermineNumSeqFeatures(d['seqFeatures'])
	modelSpecs['n_in_matrix'] = DetermineNumMatrixFeatures(d['matrixFeatures']) + DetermineNumMatrixFeatures(d['matrixFeatures_nomean'])

        if d.has_key('embedFeatures'):
                modelSpecs['n_in_embed'] = d['embedFeatures'].shape[1]

def CalcFeatureExpectBySampling(metaData, modelSpecs):
	seqfeatures = []
	seqweights = []
	
	matrixfeatures = []
	matrixweights = []

	embedfeatures = []
	embedweights = []

	dataLocation = DataProcessor.SampleProteinInfo(metaData)
	for loc in dataLocation:
		d = DataProcessor.LoadRealData(loc, modelSpecs, loadLabel=False)
		res = CalcFeatureExpect4OneProtein(d)
		seqfeature, seqweight, matrixfeature, matrixweight = res[:4]
		seqfeatures.append(seqfeature)
		matrixfeatures.append(matrixfeature)
		seqweights.append(seqweight)
		matrixweights.append(matrixweight)

		if len(res) == 6:
			embedfeature, embedweight = res[5:]
			embedfeatures.append(embedfeature)
			embedweights.append(embedweight)

	modelSpecs['seqFeatures_expected'] = np.average(seqfeatures, axis=0, weights=seqweights)
	modelSpecs['matrixFeatures_expected'] = np.average(matrixfeatures, axis=0, weights=matrixweights)
	modelSpecs['embedFeatures_expected'] = np.average(embedfeatures, axis=0, weights=embedweights)

## check consistency between model specification and input features of training/validation/prediction data
## numChecks: the number of randomly sampled proteins to be checked
def CheckModelNDataConsistency(model, data, numChecks=1):
        for i in range(numChecks):
                rindex = np.random.randint(0, high=len(data) )
                dataNumSeqFeatures = DetermineNumSeqFeatures(data[rindex]['seqFeatures'])
                if (model['n_in_seq'] != dataNumSeqFeatures):
                        print 'ERROR: inconsistency in n_in_seq between model and data: ', model['n_in_seq'], dataNumSeqFeatures
                        exit(1)

                rindex = np.random.randint(0, high=len(data) )
                dataNumMatrixFeatures = DetermineNumMatrixFeatures(data[rindex]['matrixFeatures']) + DetermineNumMatrixFeatures(data[rindex]['matrixFeatures_nomean'])
                if (model['n_in_matrix'] != dataNumMatrixFeatures):
                        print 'ERROR: inconsistency in n_in_matrix between model and data: ', model['n_in_matrix'], dataNumMatrixFeatures
                        exit(1)

                if data[0].has_key('embedFeatures'):
                        rindex = np.random.randint(0, high=len(data) )
                        res = (model['n_in_embed'] == data[rindex]['embedFeatures'].shape[1])
                        if not res:
                                print 'ERROR: inconsistency in n_in_embed between model and data: ', model['n_in_embed'], data[rindex]['embedFeatures'].shape[1]
                                exit(1)
        return True
