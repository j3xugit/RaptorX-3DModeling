import numpy as np

MyFloat=np.float16

## this constant is used to scale up predicted contact probability to maximize MCC and F1 values when p=0.5 is used as cutoff for binary contact classification.
## A probability value p is scaled to p^ProbScaleFactor, i.e., 0.23 is scaled to 0.5. 
## this scale factor is only used in saving the predicted contact matrix into a CASP submission file (i.e., in generating CASP.rr file)
## in CASP12, we did not do this and our MCC and F1 score are not the best (when p=0.5 is used as cutoff), although I do think such a scale-up is meaningless.
## this scale constant depends on the weight factor used in calculating loss function, so whenever the weight factor changes, this constant shall be adjusted.
ProbScaleFactor = np.log(0.5)/np.log(0.23)

## ResNet2DV21 and ResNet2DV22 are added on March 8, 2018
## ResNet2DV21 is the same as ResNet2D. ResNet2DV23 is recommended.
## ResNet2DV23 is almost same as ResNet2D except that the former has removed unused batch normalization layers (and parameters)
## ResNet2DV22 differs from ResNet2DV21 in that the former has two batch norm layers in each residual block while the latter has only one
## ResNet2DV22 seems to be better than ResNet2DV21, but maybe this depends on training algorithm and learning rate
allNetworks = ['ResNet2D', 'ResNet2DV21', 'ResNet2DV22', 'ResNet2DV23', 'DilatedResNet2D']

allActivations = ['RELU', 'TANH', 'ELU'] 

allAlgorithms = ['SGDM', 'SGDM2', 'Adam', 'SGNA', 'AdamW', 'AdamWAMS', 'AMSGrad']

## embedding modes are used to convert 1d sequence to 2d matrix
#allEmbeddingModes = ['SeqOnly', 'Seq+SS', 'Profile', 'OuterCat']
allEmbeddingModes = ['SeqOnly', 'Seq+SS', 'OuterCat']

def UseSampleWeight(modelSpecs):
	return modelSpecs.has_key('UseSampleWeight') and modelSpecs['UseSampleWeight']

def UseOneHotEncoding(modelSpecs):
	return modelSpecs.has_key('OneHotEncoding') and modelSpecs['OneHotEncoding']

def UseRgAsSequentialFeature(modelSpecs):
	return modelSpecs.has_key('UseRgAsSeqFeature') and modelSpecs['UseRgAsSeqFeature']

def NoOldLocationFeatures(modelSpecs):
	return modelSpecs.has_key('NoOldPosFeatures') and modelSpecs['NoOldPosFeatures']

def UseNewLocationFeatures(modelSpecs):
	return modelSpecs.has_key('UseNewPosFeatures') and modelSpecs['UseNewPosFeatures']

def UseSequentialFeatures(modelSpecs):
	return modelSpecs.has_key('UseSequentialFeatures') and modelSpecs['UseSequentialFeatures']

def UseCCMZ(modelSpecs):
	return modelSpecs.has_key('UseCCMZ') and modelSpecs['UseCCMZ']

def UseRawCCM(modelSpecs):
	return modelSpecs.has_key('UseRawCCM') and modelSpecs['UseRawCCM']

def UsePSICOV(modelSpecs):
	return modelSpecs.has_key('UsePSICOV') and modelSpecs['UsePSICOV']

def UseMI(modelSpecs):
	if modelSpecs.has_key('UseOtherPairs') and modelSpecs['UseOtherPairs']:
		return True
	if not modelSpecs.has_key('UseMI'):
		return False
	return modelSpecs['UseMI']

def UseContactPotential(modelSpecs):
	if modelSpecs.has_key('UseOtherPairs') and modelSpecs['UseOtherPairs']:
		return True
	if not modelSpecs.has_key('UseContactPotential'):
		return False
	return modelSpecs['UseContactPotential']

def NoWeight4Range(modelSpecs):
	if not modelSpecs.has_key('NoWeight4Range'):
		return False
	return modelSpecs['NoWeight4Range']

def NoWeight4Label(modelSpecs):
	if not modelSpecs.has_key('NoWeight4Label'):
		return False
	return modelSpecs['NoWeight4Label']

def TrainByRefLoss(modelSpecs):
	if not modelSpecs.has_key('TrainByRefLoss'):
		return False
	return modelSpecs['TrainByRefLoss']

def UseRefState(modelSpecs):
	if not modelSpecs.has_key('UseRefState'):
		return False
	return modelSpecs['UseRefState']

def EmbeddingUsed(modelSpecs):
	if not modelSpecs.has_key('seq2matrixMode'):
		return False
	return any( k in modelSpecs['seq2matrixMode'] for k in ('SeqOnly', 'Seq+SS') )

def UseTemplate(modelSpecs):
	return modelSpecs.has_key('UseTemplate') and modelSpecs['UseTemplate']

def InTPLMemorySaveMode(modelSpecs):
	if not modelSpecs.has_key('TPLMemorySave'):
		return False
	return modelSpecs['TPLMemorySave']

## get the usage mode of Evoltionary Coupling information. extraCCMmode is the old key name
## ECInfo is the new key name
def ParseExtraCCMmode(modelSpecs):
	if (not modelSpecs.has_key('extraCCMmode')) and (not modelSpecs.has_key('ECInfo')):
		return False, False, False, False, False

	if modelSpecs.has_key('ECInfo'):
		exCCMmode = np.int32(modelSpecs['ECInfo'])
	else:
		exCCMmode = np.int32(modelSpecs['extraCCMmode'])

	bUseCCMFnorm = (exCCMmode &1 )>0
	bUseCCMsum = (exCCMmode &2 )>0
	bUseCCMraw = (exCCMmode &4 )>0
	bUseFullMI = (exCCMmode &8 )>0
	bUseFullCov = (exCCMmode &16 )>0

	return bUseCCMFnorm, bUseCCMsum, bUseCCMraw, bUseFullMI, bUseFullCov

def ParseESMmode(modelSpecs):
	if not modelSpecs.has_key('ESM'):
		return None

	if modelSpecs['ESM'] == '':
		return None

	if not isinstance(modelSpecs['ESM'], str):
		return [ np.int32(modelSpecs['ESM']) ]

	import re
	#layers = [ np.int32(field) for field in re.split(',| ', modelSpecs['ESM']) ]
	layers = []
	segments = re.split(',| ', modelSpecs['ESM'])
	for seg in segments:
		fields = seg.split('t')
		if len(fields) == 1:
			layers.append(np.int32(fields[0]))
		elif len(fields) == 2:
			start = np.int32(fields[0])
			end = np.int32(fields[1])
			layers.extend( np.arange(start, end+1) )
		else:
			print 'ERROR: inocrrect format for ESM option:', seg
			exit(1)
			
	if not bool(layers):
		return None
	return layers

def ParseAttentionMode(modelSpecs):
	if not modelSpecs.has_key('Attention'):
		return None
	if modelSpecs['Attention'] is None:
		return None

	fields = modelSpecs['Attention'].split('+')
	UseAvg = True

	UseMax = False
	if 'UseMax' in fields:
		UseMax= True
	assert UseAvg or UseMax

	UseFC = True
	if 'UseConv' in fields:
		UseFC = False
	return (UseAvg, UseMax, UseFC)

def GetBoundingBoxOffset(modelSpecs):
	if not modelSpecs.has_key('boundingBoxOffset'):
		return None
	return modelSpecs['boundingBoxOffset']

## when the pairwise feature has hundreds of channels, maybe we shall compress it before feeding it into 2D CNN
def CompressMatrixInput(modelSpecs):
	return modelSpecs.has_key('CompressMatrixInput') and modelSpecs['CompressMatrixInput']

## a response is what we want to predict using our deep network, e.g., distance and orientation

## definition for responses that can be predicted by this deep learning package
allAtomPairNames = ['CbCb', 'CaCa', 'CgCg', 'CaCg', 'NO']
allAtomPairTypes = allAtomPairNames
allDistLabelNames = allAtomPairNames + ['HB', 'Beta']

## symmetric label names
symAtomPairNames = ['CbCb', 'CaCa', 'CgCg', 'Beta']

## inter-residue orientation
TwoRAngleNames = ['Ca1Cb1Cb2']
TwoRDihedralNames = ['Ca1Cb1Cb2Ca2','N1Ca1Cb1Cb2']
TwoROriNames = ['Ca1Cb1Cb2Ca2','N1Ca1Cb1Cb2','Ca1Cb1Cb2']

## orientation formed by four Ca atoms. The first two Ca atoms are sequentially adjacent, so are the second two.
FourCaDihedralNames = ['Ca1Ca2Ca3Ca4', 'Ca1Ca2Ca4Ca3']
FourCaAngleNames = ['Ca1Ca2Ca3', 'Ca1Ca2Ca4']
FourCaOriNames = ['Ca1Ca2Ca3Ca4', 'Ca1Ca2Ca3', 'Ca1Ca2Ca4Ca3', 'Ca1Ca2Ca4']

allDihedralNames = ['Ca1Cb1Cb2Ca2','N1Ca1Cb1Cb2','Ca1Ca2Ca3Ca4', 'Ca1Ca2Ca4Ca3']
allAngleNames = ['Ca1Cb1Cb2', 'Ca1Ca2Ca3', 'Ca1Ca2Ca4']

allOrientationNames = allDihedralNames + allAngleNames
symOrientationNames = ['Ca1Cb1Cb2Ca2', 'Ca1Ca2Ca4Ca3']

def IsSymmetricLabel( labelName ):
	return  (labelName in symAtomPairNames) or (labelName in symOrientationNames)

allPairwiseLabelNames = allDistLabelNames + allOrientationNames
allLabelNames = allPairwiseLabelNames

## a label name string can be CbCb or CbCb+CaCa or CbCb+CaCa+CgCg or AllAP or Ca1Cb1Cb2Ca2 or Ca1Cb1Cb2Ca2+N1Ca1Cb1Cb2 or AllOri where AllAP represents allAtomPairNames and AllOri represents allOrientationNames
def IsValidLabelName(labelName):
	allnames = set( allLabelNames )
	return labelName in allnames

## some abbreviations are defined in this function for convenience
def ParseLabelNames(nameStr):
        fields = nameStr.split('+')
	names = []
	for f in fields:
		if f.upper() == 'AllAP'.upper():
        		names.extend(allAtomPairNames)
		elif f.upper() == 'AllOri'.upper():
			names.extend(allOrientationNames)

		elif f.upper() == 'TwoROri'.upper():
			names.extend(TwoROriNames)
		elif f.upper() == 'TwoRDihedral'.upper():
			names.extend(TwoRDihedralNames)
		elif f.upper() == 'TwoRAngle'.upper():
			names.extend(TwoRAngleNames)

		elif f.upper() == 'FourCaOri'.upper():
			names.extend(FourCaOriNames)
		elif f.upper() == 'FourCaDihedral'.upper():
			names.extend(FourCaDihedralNames)
		elif f.upper() == 'FourCaAngle'.upper():
			names.extend(FourCaAngleNames)

		elif f.upper() == 'AllDihedral'.upper():
			names.extend(allDihedralNames)
		elif f.upper() == 'AllAngle'.upper():
			names.extend(allAngleNames)
		elif IsValidLabelName(f):
			names.append(f)
		else:
			print 'ERROR: unsupported label name ', f, ' in ', nameStr 
			exit(1)
	return names

## check if a list or set of label names has any orientation labels
def HasOrientationNames(labelNames):
	lnames = set(labelNames).intersection( set(allOrientationNames) )
	return bool(lnames)

def HasTwoROriNames(labelNames):
	lnames = set(labelNames).intersection( set(TwoROriNames) )
	return bool(lnames)

def HasFourCaOriNames(labelNames):
	lnames = set(labelNames).intersection( set(FourCaOriNames) )
	return bool(lnames)
	

## 1) In a distance matrix, -1 represents an invalid distance (i.e, at least one residue has no valid 3D coodinates in PDB file) and a positive value represents a valid distance
## 2) In the beta-pairing (Beta) or hydrogen-bonding (HB) matrix, -1 is used to indicate no valid distance between two Cbeta atoms
## we also use a value 100 + real_distance to indicate that one entry does not form a beta pairing or hydrogen bond but has valid distance
## when one entry in a Beta or HB matrix forms a beta pairing or hydrogen bond, this entry contains the real distance of the Cbeta atoms.
## 3) In the dihedral/angle matrices, a value > 180 is used to indicate an invalid entry


## XXC represents one discretization scheme for distance/orientation.
discreteLabelsBasic = ['56C', '52C', '49C', '47C', '42C', '37C', '31C', '25C', '24C', '19C', '10C', '16C', '14C', '13C', '12C', '3C', '2C']

## treat all invalid entry in the distance/orientation matrix as a separate label and counted in calculating loss function
discreteLabelsPlus = [ label + 'Plus' for label in discreteLabelsBasic ]

## treat all invalid entry in the distance/orientation matrix as a separate label, but not counted in calculating loss function
discreteLabelsMinus = [ label + 'Minus' for label in discreteLabelsBasic ]

discreteLabels = discreteLabelsBasic + discreteLabelsPlus + discreteLabelsMinus

allDiscreteLabelTypes = [ 'Discrete' + label for label in discreteLabels ]
allContinuousLabelTypes = ['Normal', 'LogNormal'] 

allLabelTypes = allDiscreteLabelTypes + allContinuousLabelTypes

def IsContinuousLabel(labelType):
	ret = [ labelType.startswith(cl) for cl in allContinuousLabelTypes ]
	return any(ret)

def IsDiscreteLabel(labelType):
	return labelType.startswith('Discrete')

#a repsonse is defined by one label name and one label type, e.g., CbCb_Discrete34C
def Response2LabelName(response):
        return response.split('_')[0]

def Response2LabelType(response):
        return response.split('_')[1]

def Response2subType(response):
        return response.split('_')[1][len('Discrete'):]

def ParseResponse(response):
	fields = response.split('_')
	labelName = fields[0]
	labelType = fields[1]
	if labelType.startswith('Discrete'):
		subType = labelType[len('Discrete'):]
	elif labelType.startswith('Normal'):
		subType = labelType[len('Noraml'):]
	elif labelType.startswith('LogNormal'):
		subType = labelType[len('LogNoraml'):]
	else:
		print 'unsupported labelType: ', labelType
		exit(1)

	return labelName, labelType, subType

def GetAllLabelNames(responses):
	labelNames = [ Response2LabelName(response) for response in responses ]
	return labelNames
	
def HasOrientationResponses(responses):
	labelNames = [ Response2LabelName(response) for response in responses ]
	return HasOrientationNames(labelNames)

## the number of dimensions for a predicted value of a response.
## currently only 1d Normal or LogNormal is implemented. 
## to support the 2d Normal or LogNormal, need to check out other places to make sure it is correct
## Normal indicates 1d normal variable and Normal2d indicates 2d normal variable
responseValueDims = dict()
responseValueDims['Normal'] = 1
responseValueDims['Normal2d'] = 2
responseValueDims['Normal2d2'] = 2
responseValueDims['Normal2d4'] = 2

responseValueDims['LogNormal'] = 1
responseValueDims['LogNormal2d'] = 2
responseValueDims['LogNormal2d2'] = 2
responseValueDims['LogNormal2d4'] = 2

## for a discrete response, only one label is needed for one predicted value 
for labelType in allDiscreteLabelTypes:
	responseValueDims[labelType]=1

def GetResponseValueDims(response):
        labelName, labelType, subType = ParseResponse(response)
        if IsContinuousLabel(labelType):
                return responseValueDims[labelType]

	return 1


## the number of parameters needed to define a probability distribution function of a response
## A Normal/LogNormal distribution may have different number of parameters depending on if we want to predict its variance and correlation
## When the response is discrete, it is the number of discrete labels
responseProbDims = dict()
responseProbDims['Normal']=2
responseProbDims['Normal2d']=5
responseProbDims['Normal2d2']=2
responseProbDims['Normal2d4']=4

responseProbDims['LogNormal']=2
responseProbDims['LogNormal2d']=5
responseProbDims['LogNormal2d2']=2
responseProbDims['LogNormal2d4']=4

## discretization schemes for distance and orientation

## When Plus or Minus is used, the invalid distance/dihedral/angle is separated from the last distance/dihedral/angle bin.
## Otherwise, the invalid entry is merged into the last bin of distance/dihedral/angle.
## When Plus is used, the invalid entry is counted towards loss function in training, 
## but when Minus is used, the invalid entry is not taken into consideration by the loss function

## usually, we use neither Plus or Minus

distCutoffs = {}

distCutoffs['56C'] = np.array( [0] + np.linspace(2.0, 20.0, num=55).tolist() ).astype(np.float32)
distCutoffs['47C'] = np.array( [0] + np.linspace(2.0, 20.0, num=46).tolist() ).astype(np.float32)
distCutoffs['46C'] = np.array( [0] + np.linspace(2.0, 20.0, num=46).tolist() ).astype(np.float32)

distCutoffs['42C'] = np.array( [0] + np.linspace(2.0, 22.0, num=41).tolist() ).astype(np.float32)
distCutoffs['38C'] = np.array( [0] + np.linspace(2.0, 20.0, num=37).tolist() ).astype(np.float32)
distCutoffs['36C'] = np.array( [0] + np.linspace(3.0, 20.0, num=35).tolist() ).astype(np.float32)
distCutoffs['34C'] = np.array( [0] + np.linspace(4.0, 20.0, num=33).tolist() ).astype(np.float32)

distCutoffs['10C'] = np.array( [0] + np.linspace(4.0, 20.0, num=9).tolist() ).astype(np.float32)

distCutoffs['52C'] = np.array( [0] + np.linspace(4.0, 16.5, num=51).tolist() ).astype(np.float32)
distCutoffs['25C'] = np.array( [0] + np.linspace(4.5, 16.0, num=24).tolist() ).astype(np.float32)
distCutoffs['14C'] = np.array( [0] + np.linspace(4.0, 16.0, num=13).tolist() ).astype(np.float32)
distCutoffs['13C'] = np.array( [0] + np.linspace(5.0, 16.0, num=12).tolist() ).astype(np.float32)
distCutoffs['12C'] = np.array( [0] + np.linspace(5.0, 15.0, num=11).tolist() ).astype(np.float32)

distCutoffs['3C']  = np.array( [0, 8, 15] ).astype(np.float32)
distCutoffs['2C']  = np.array( [0, 8 ] ).astype(np.float32)

distCutoffKeys = distCutoffs.keys()
for k in distCutoffKeys:
	distCutoffs[ k + 'Plus' ] = distCutoffs[k]
	distCutoffs[ k + 'Minus' ] = distCutoffs[k]

## for hydron-bonding: the maximum Cbeta distance of two residues forming a hydrogen bond is slightly more than 9 Angstrom
MaxHBDistance = 9.5
distCutoffs_HB = {}
distCutoffs_HB['2C'] = np.array( [0, MaxHBDistance] ).astype(np.float32)
distCutoffs_HB['2CPlus'] = distCutoffs_HB['2C']
distCutoffs_HB['2CMinus'] = distCutoffs_HB['2C']

## for 2-body orientation, when the distance between two residues (i.e., Cb atoms) > a threshold, all the dihedral/angles are grouped into the largest bin of valid entry
## for 4-body orientation, when the distance between the 1st and 3rd Ca atoms > a threshold, all the dihedral/angles are grouped into the largest bin of valid entry
## When Plus or Minus is used, an invalid entry is separated from the largest bin of valid entry; otherwise, merged into the largest bin of valid entry
dihedralCutoffs = {}
dihedralCutoffs['49C'] = np.array( np.linspace(-180., 180., num=49).tolist() ).astype(np.float32)
dihedralCutoffs['41C'] = np.array( np.linspace(-180., 180., num=41).tolist() ).astype(np.float32)
dihedralCutoffs['37C'] = np.array( np.linspace(-180., 180., num=37).tolist() ).astype(np.float32)
dihedralCutoffs['31C'] = np.array( np.linspace(-180., 180., num=31).tolist() ).astype(np.float32)
dihedralCutoffs['25C'] = np.array( np.linspace(-180., 180., num=25).tolist() ).astype(np.float32)

dihedralKeys = dihedralCutoffs.keys()
for k in dihedralKeys:
	dihedralCutoffs[ k + 'Plus' ] = dihedralCutoffs[k]
	dihedralCutoffs[ k + 'Minus' ] = dihedralCutoffs[k]

angleCutoffs = {}
angleCutoffs['25C'] = np.array( np.linspace(0., 180., num=25).tolist() ).astype(np.float32)
angleCutoffs['19C'] = np.array( np.linspace(0., 180., num=19).tolist() ).astype(np.float32)
angleCutoffs['16C'] = np.array( np.linspace(0., 180., num=16).tolist() ).astype(np.float32)
angleCutoffs['13C'] = np.array( np.linspace(0., 180., num=13).tolist() ).astype(np.float32)

angleKeys = angleCutoffs.keys()
for k in angleKeys:
	angleCutoffs[ k + 'Plus' ] = angleCutoffs[k]
	angleCutoffs[ k + 'Minus' ] = angleCutoffs[k]

def GetCutoffs(response):
	labelName, labelType, subType = ParseResponse(response)
	if not labelType.startswith('Discrete'):
		return None

	if labelName in allAtomPairNames:
		return distCutoffs[subType]
	elif labelName in ['HB']:
		return distCutoffs_HB[subType] 
	elif labelName in allDihedralNames:
		return dihedralCutoffs[subType]
	elif labelName in allAngleNames:
		return angleCutoffs[subType]
	else:
		print 'ERROR: unsupported response in GetCutoffs: ', response
		exit(1)

def GetResponseProbDims(response):
	labelName, labelType, subType = ParseResponse(response)
	if IsContinuousLabel(labelType):
		return responseValueDims[labelType]

	cutoff = GetCutoffs(response)
	numLabels = len(cutoff)
	if subType.endswith('Plus') or subType.endswith('Minus'):
		numLabels += 1
	return numLabels

def GetLabelMinMaxValues(labelName):
	if labelName in allAtomPairNames:
		return 0, np.finfo(np.float16).max
	elif labelName in allDihedralNames:
		return -180, 180
	elif labelName in allAngleNames:
		return 0, 180
	else:
		print 'ERROR: unsupported Label in GetLabelMinMaxValues: ', labelName
		exit(1)


## the maximum Cbeta distance of two residues forming a beta pair is approximately 8 Angstrom
MaxBetaDistance = 8.0

## the distance cutoff for Cbeta-Cbeta contact definition
ContactDefinition = 8.00001

## when the distance between two atoms is beyond this cutoff, we assume they have no interaction at all
InteractionLimit = 15.00001

## the top ratio*L predicted values are evaluated where L is the sequence length
topRatios = dict()
for apt in allAtomPairNames:
	topRatios[apt] = 0.5
topRatios['HB'] = 0.1
topRatios['Beta'] = 0.1
for apt in allOrientationNames:
	topRatios[apt] = 5.

## the size of queue for multiprocessing
def QSize(modelSpecs):
	if modelSpecs.has_key('QSize'):
		return modelSpecs['QSize']
	return 200

def NumTrainDataLoaders(modelSpecs):
	if modelSpecs.has_key('NumDataLoaders'):
		return np.int32(modelSpecs['NumDataLoaders'])
	return 1
