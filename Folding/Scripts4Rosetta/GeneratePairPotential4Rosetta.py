import os
import sys
import cPickle
import numpy as np
import getopt

from DL4DistancePrediction4 import DistanceUtils
import DL4DistancePrediction4.config as config
from DL4DistancePrediction4.config import ParseResponse
from Common.SelectAtoms import SelectAtomPair, SelectAtoms4Orientation
from Common.SequenceUtils import LoadFASTAFile

def Usage():
	print 'python GeneratePairPotential4Rosetta.py [-a labelNames] [ -f funcType] [-s minSeqSep] [ -c potentialCutoff ] [-t topRatio] [-q querySeqFile] [-d savefolder] potential_PKL'
	print '     potential_PKL: the raw potential matrix in PKL format, a tuple of 5 items: target, sequence, potential (dict), cutoffs (dict), validProbs(dict)'
	print '     -a: labelNames, e.g., CbCb, CbCb+TwoROri(default), CbCb+AllOri where TwoROri represents inter-residue orientation and AllOri denotes all orientation types'
	print '     -f: potential function type: SPLINE(default), only Rosetta SPLINE is supported now'
	print ''
	print '     -s: use potential for atom pairs with at least this sequence separation, i.e., only atoms of two residues i and j with abs(i-j)>=this value is considered'
	print '		e.g., 1+2, where the first is for distance (default 1) and the 2nd for orientation (default 2)'
	print '		When only one number is specified, it is used for distance only'
	print ''
	print '     -c: threshold for distance and orientation potential, e.g., 0.04+20, where the first value is for orientation and the 2nd for distance'
	print'		for distance potential, when all potential of one atom/residue pair is worse (larger) than this value (default max of float32), ignore this pair'
	print '		for orientation potential, when the absolute of all potential of one pair is smaller than this value (default 0.04), ignore this pair'
	print '	    -t: the average number of constraints per residue/atom employed for each type of potential, e.g., 20, 20+30, where the first number (default 25) is for orientation and the 2nd (default seq length) for distance'
	print '	     the -c and -t options can be combined to filter atom/residue pairs. However, currently they are not implemented for distance potential'
	print ''
	print '     -q: the protein sequence file. If provided, check the consistency between this file and the sequence in potential_PKL'		
	print '	    -d: the folder for saving the result file, default current work directory'
	print '		one file XXX.pairPotential4Rosetta.SPLINE.txt and one subfolder SplinePotential4XXX will be generated under this folder where XXX is the base name of potential_PKL'
	print '	    Example: python GeneratePairPotential4Rosetta.py -a CbCb+AllOri -s 2+2 T0950.pairPotential.DFIRE.18.1.61.pkl'

## this function writes the constraints into Rosetta format
## target is the protein name
## constraints is a list of python dict and each dict corresponds to one constraint
def WriteSplineConstraints(constraints, savefile=None, savefolder4histfile=None):
	if savefile is None:
		print 'ERROR: please specify the save file for constaints!'
		exit(1)
		
	if savefolder4histfile is None:
		print 'ERROR: please specify the save file for constaints!'
		exit(1)
	histfileDir = savefolder4histfile
	if not os.path.isdir(histfileDir):
		os.mkdir(histfileDir)

	expVal = 0.
	weight = 1.

	numIgnored = 0

	potStrs = []
	for constraint in constraints:
		## write histogram to histfile
		response = constraint['response']
		labelName, _, _ = config.ParseResponse(response)
		x = constraint['x']
		y = constraint['y']
		if not np.isfinite(y).all():
			print 'WARNING: ignore one constraint since it may have an NaN or infinite value:', constraint
			numIgnored += 1
			continue
			
		atomNums = [ str(i) for i in constraint['atomNums'] ]
		atomNumStr = '-'.join(atomNums)

		histfile = os.path.join(histfileDir, response + '-' + atomNumStr + '.potential.txt')
		xStr = '\t'.join(['x_axis'] + [ "{:.4f}".format(e) for e in x ] )
		yStr = '\t'.join(['y_axis'] + [ "{:.4f}".format(e) for e in y ] )
		with open(histfile, 'w') as fh:
			fh.write('\n'.join([xStr, yStr]) + '\n')

                #potStr = ' '.join(['Angle', atom1.upper(), str(i+1), atom2.upper(), str(i+2), atom3.upper(), str(j+1), 'SPLINE', description, histfile] + [ "{:.4f}".format(e) for e in [expVal, weight, binWidth] ] )
		potStrList = [ constraint['type'] ]
		for name, number in zip(constraint['atoms'], atomNums):
			potStrList.extend([name.upper(), number])
		potStrList.append('SPLINE')
		potStrList.append(labelName)
		potStrList.append(histfile)

		potStrList.extend( ['0', '1', "{:.6f}".format(constraint['binWidth']) ])
		potStr = ' '.join(potStrList)

		potStrs.append(potStr)

	if numIgnored > 100:
		print 'ERROR: too many constraints are ignored:', numIgnored
		exit(1)

	if len(potStrs)>0:
        	with open(savefile, 'w') as fh:
        		fh.write('\n'.join(potStrs) + '\n')

	return potStrs
		
## generate orientation constraints. Only top topRatio *L constraints are generated where L is the sequence length
def GenerateSplinePotential4Orientation(potData, labelNames=None, topRatio=20, minSeqSep=2, potThreshold=0.04):

	if len(potData) < 5:
		print 'ERROR: validProb is needed for orientation potential data'
		exit(1)
	target, sequence, potential, cutoffs, validProbs = potData[:5]

	allConstraints = []

        for response, pot in potential.iteritems():
		description = response
		#print 'response=', response
		labelName, labelType, subType = ParseResponse(response)
		if labelName not in config.allOrientationNames:
			continue
		if labelName not in labelNames:
			continue

		## cutoffs has the discretization scheme of orientation angles. Here we require all bins have the same width
                x = cutoffs[response] 
		binWidths = [ b-a for a, b in zip(x[:-1], x[1:]) ]
                binWidth = binWidths[0]
		assert all([binWidth == b for b in binWidths])	

		binWidth = binWidth * np.float32(np.pi)/180
		x = [ e * np.float32(np.pi) / 180 for e in x ]

		validProb = validProbs[response]
		assert pot.shape[:2] == validProb.shape

		ind = np.unravel_index(np.argsort(-validProb, axis=None), validProb.shape)
		ratio = min(topRatio, validProb.shape[0]-3)
		topK = np.int32(ratio * validProb.shape[0])

		for i, j in zip(ind[0][:topK], ind[1][:topK]):
			offset = abs(i-j)
			if offset < minSeqSep:
				continue
			if config.IsSymmetricLabel(labelName) and i>j:
				continue

        		atoms = SelectAtoms4Orientation(sequence, i, j, labelName)
			if atoms is None:
				continue

                        y = pot[i, j]
			## if the absolute value of the best potential is too small, then skip
			ymax = np.amax( np.absolute(y) )
			if (ymax < potThreshold):
				continue

			## each value in original y is the potential for one interval, corresponding to the middle of an interval. 
			## the cutoff points in x is the boundary of intervals instead of middle points, so we need to shift y by half of the binwidth
			## here we calculate the approximate value for each cutoff point in x. y2 and y3 are cycle shifts of y
			y2 = [ y[-1] ] + y.tolist()
			y3 = y.tolist() + [ y[0] ]

			## when one angle is symmetric, it appears only once in the constraint set, 
			## when one angle is not symmetric, it appears twice in the constraint set, so we divide its potential by 2 more
			if config.IsSymmetricLabel(labelName):
				y = [ (a+b)/2 for a, b in zip(y2, y3) ]
			else:
				y = [ (a+b)/4 for a, b in zip(y2, y3) ]

			constraint = dict()
			constraint['x'] = x
			constraint['y'] = y
			constraint['response'] = response
			constraint['binWidth'] = binWidth 

			if labelName in config.allAngleNames:
				constraint['type'] = 'Angle'
				assert len(atoms) == 3
				atom1, atom2, atom3 = atoms
				constraint['atoms'] = [ atom1, atom2, atom3]
				if labelName == 'Ca1Cb1Cb2':
					constraint['atomNums'] = [ i+1, i+1, j+1]
				elif labelName == 'Ca1Ca2Ca3':
					constraint['atomNums'] = [ i+1, i+2, j+1]
				elif labelName == 'Ca1Ca2Ca4':
					constraint['atomNums'] = [ i+1, i+2, j+2]
				else:
					print 'ERROR: unsupported angle names: ', labelName
					exit(1)

			elif labelName in config.allDihedralNames:
				constraint['type'] = 'Dihedral'
				assert len(atoms) == 4
				atom1, atom2, atom3, atom4 = atoms
				constraint['atoms'] = [ atom1, atom2, atom3, atom4]
				if labelName == 'Ca1Cb1Cb2Ca2':
					constraint['atomNums'] = [ i+1, i+1, j+1, j+1]
				elif labelName == 'N1Ca1Cb1Cb2':
					constraint['atomNums'] = [ i+1, i+1, i+1, j+1]
				elif labelName == 'Ca1Ca2Ca3Ca4':
					constraint['atomNums'] = [ i+1, i+2, j+1, j+2]
				elif labelName == 'Ca1Ca2Ca4Ca3':
					constraint['atomNums'] = [ i+1, i+2, j+2, j+1]
				else:
					print 'ERROR: unsupported dihedral names: ', labelName
					exit(1)

			else:
				print 'ERROR: unsupported labelName in GeneratePotential for orientation: ', labelName
				exit(1)

                        allConstraints.append(constraint)

	return allConstraints

## generate distance potential constraints for PyRosetta
## currently topRatio and potThreshold are not used for distance potential
def GenerateSplinePotential4Distance(potData, labelNames=['CbCb'], topRatio=None, minSeqSep=1, potThreshold=np.finfo(np.float32).max, barrier=1.0):

	target, sequence, potential, distCutoffs = potData[:4]
	if len(potData)>4:
		validProb = potData[4]
	else:
		print 'WARNING: it is better to provide validProb for distance potential'
		validProb = None

	allConstraints = []
        for response, pot in potential.iteritems():
		#print 'response=', response

		description = response
		labelName, labelType, subType = ParseResponse(response)

		if labelName not in config.allAtomPairNames:
			continue
		if labelName not in labelNames:
			continue

                #x = distCutoffs[response][1:]
                x = distCutoffs[response]
		binWidths = [ b-a for a, b in zip(x[1:-1], x[2:]) ]
                binWidth = np.average(binWidths)
		#assert all([ (binWidth-b)<0.0001 for b in binWidths])	

		## here we add repulsion to reduce steric clashes
		## for CaCa and CbCb, the minimum distance is 3.6A. The penalty is 3 in [2, 3.6] and 10 in [0, 2]

		## for NO, the minimum distance is 2, and the penalty is 1 for [2, 2.4] and 5 for [0, 2]
		firstMinDist = 2
		secondMinDist = 3.6
		yPenalty = [10, 4, 0.5]

		if labelName == 'NO':
			secondMinDist = 2.4
			yPenalty = [8, 3, 0.2]
		xPrefix = [ 0, firstMinDist, secondMinDist ]

		## find the index of the 2nd min distance in x, i.e., x[secondLabel] <=secondMinDist < x[secondLabel+1]
		secondLabel = DistanceUtils.LabelsOfOneDistance(secondMinDist + 0.0001, x)
		#print 'secondLabel=', secondLabel
		assert secondLabel >= 1
		assert secondLabel < len(distCutoffs[response])

		xk = [ (a+b)/2 for a, b in zip(x[secondLabel:-1], x[secondLabel+1:]) ]
		xk.append(x[-1] + binWidth/2.)
		xk = xPrefix + xk
		#print 'xk=', xk

		#LabelOfAdjacentCAs = DistanceUtils.LabelsOfOneDistance(4.50001, x)
		"""
		## the first interval of distCutoffs, i.e., [0, x[0]), usually has width > binWidth
		## Here we split the first interval into several bins with the same binwidth
		## Assume that x[0] = binWidth times an integer
		bstep = barrier * binWidth
		xPrefix = np.linspace(-binWidth, x[0]-binWidth, np.rint(x[0]/binWidth).astype(np.int32)+1 )
		xPrefix = xPrefix.tolist()

                xk = xPrefix + x.tolist()
		xk2 = [ (a+b)/2 for a, b in zip(xk[:-1], xk[1:]) ]
		xk2.append(x[-1]+binWidth/2.)
		xk = xk2

		yPrefix = np.arange( len(xPrefix)-1, -1, -1) * bstep
		yPrefix = yPrefix.tolist()
		"""

                size = pot.shape
		residuePairs = []

                for i in xrange(size[0]):
			jstart = i+minSeqSep
			if not config.IsSymmetricLabel(labelName):
				jstart=0

                        for j in xrange(jstart, size[1]):
				offset = abs(i-j)
				if offset < minSeqSep:
					continue
				residuePairs.append( (i, j) )

		## always use repulsion potential for two sequentially adjacent Ca atoms	
		if minSeqSep > 1 and labelName == 'CaCa':
			residuePairs.extend( [(i, i+1) for i in xrange(size[0]-1) ] )

		for i, j in residuePairs:
                	y = pot[i, j]

			"""
			## y[0] is the potential for the first interval [0, x[0]). We increase potential for distance < x[0] for every binWidth Angstrom
			yPrefix2 = [ y[0] + ye for ye in yPrefix ]
			yk = yPrefix2 + y[1:].tolist()
			"""
			yPrefix = [ max(y[secondLabel], 0) + ye for ye in yPenalty ]
			y2 = y.tolist()
			yk = yPrefix + y2[secondLabel:]

                        assert len(xk) == len(yk), 'xk and yk length does not match for ' + labelName + ' and residues ' + str(i) + ' ' + str(j)

			## when one atom pair is not symmetric (e.g., NO), it appears twice in the constraint set, so we divide its potential by 2
			if not config.IsSymmetricLabel(labelName):
				yk = [ ye/2. for ye in yk]

                        atom1, atom2 = SelectAtomPair(sequence, i, j, labelName)

			constraint = dict()
			constraint['x'] = xk
			constraint['y'] = yk
			constraint['response'] = response
			constraint['binWidth'] = binWidth

			constraint['type'] = 'AtomPair'
			constraint['atoms'] = [ atom1, atom2]
			constraint['atomNums'] = [ i+1, j+1]

                        allConstraints.append(constraint)

	return allConstraints

## convert potential matrix to Rosetta constraints; potData is a tuple of 5 items: name, sequence, pairPotential, discretization scheme and validProbs
def GenerateSplinePotential(potData, labelNames=['CbCb', 'CaCa', 'NO'] + config.TwoROriNames, topRatio4Dist=np.iinfo(np.int32).max, topRatio4Ori=25, minSeqSep4Dist=1, minSeqSep4Ori=2, distPotThreshold=np.finfo(np.float32).max, oriPotThreshold=0.04):
	allConstraints = []

        distLabelNames = [ lname for lname in labelNames if lname in config.allAtomPairNames ]
	if len(distLabelNames) < 1:
		print 'WARNING: there is no distance potential in the potential data'
	else:
		## generate distance potential for PyRosetta
        	allDistConstraints = GenerateSplinePotential4Distance(potData, labelNames=distLabelNames, topRatio=topRatio4Dist, minSeqSep=minSeqSep4Dist, potThreshold=distPotThreshold)
		allConstraints.extend(allDistConstraints)

        oriLabelNames = [ lname for lname in labelNames if lname in config.allOrientationNames ]
        if len(oriLabelNames)>0:
        	## generate orientation potential for PyRosetta
                allOriConstraints = GenerateSplinePotential4Orientation(potData, labelNames=oriLabelNames, topRatio=topRatio4Ori, minSeqSep=minSeqSep4Ori, potThreshold=oriPotThreshold)
                allConstraints.extend( allOriConstraints)

	if len(allConstraints) < 1:
		print 'ERROR: there is no Rosetta constraints generated from the potential data'
		exit(1)

	return allConstraints

def main(argv):

	inputFile = None

        labelNames = ['CbCb'] + config.TwoROriNames
	seqSep4Dist = 1
	seqSep4Ori = 2
	distPotThreshold = np.finfo(np.float32).max
	oriPotThreshold = 0.04
	barrier = 1.0

	funcType = 'SPLINE'
	allFuncTypes = set(['SPLINE'])

	topRatio4Ori = 25
	topRatio4Dist = np.iinfo(np.int32).max

	savefolder = os.getcwd()

	if len(argv) < 1:
		Usage()
		exit(1)
        try:
                opts, args = getopt.getopt(argv,"a:f:s:c:b:t:q:d:",["labelNames=", "funcType=", "minSeqSep=", "potentialCutoff=", "barrier=", "topRatio=", "querySeqFile=", "savefolder="])
                #print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

        if len(args) != 1:
                Usage()
                exit(1)

	inputFile = args[0]
	querySeqFile = None
	querySeq = None

        for opt, arg in opts:
                if opt in ("-a", "--labelNames"):
			labelNames = config.ParseLabelNames(arg)

	  	elif opt in ("-s", "--minSeqSep"):
			fields = arg.split('+')
			seqSep4Dist= np.int32(fields[0])
			assert seqSep4Dist > 0, "The sequence separation for a valid distance potential shall be at least 1"

			if len(fields) > 1:
				seqSep4Ori= np.int32(fields[1])
				assert seqSep4Ori > 1, "The sequence separation for a valid orientation shall be at least 2"

                elif opt in ("-c", "--potentialCutoff"):
			fields = arg.split('+')
                        oriPotThreshold = np.float32(fields[0])
			if len(fields) > 1:
				distPotThreshold= np.float32(fields[1])

                elif opt in ("-f", "--funcType"):
                        funcType = arg.upper()
			if funcType not in allFuncTypes:
				print 'ERROR: unsupported potential func type:', funcType
				exit(1)

		elif opt in ("-b", "--barrier"):
			barrier = np.float32(arg)
			assert barrier>=0

		elif opt in ("-t", "--topRatio"):
			fields = arg.split('+')
			topRatio4Ori = np.float32(fields[0])
			if len(fields) > 1:
				topRatio4Dist = np.float32(fields[1])

		elif opt in ("-q", "--querySeqFile"):
			querySeqFile = arg

		elif opt in ("-d", "--savefolder"):
			savefolder = arg
                else:
                        Usage()
                        exit(1)

	if inputFile is None:
                print 'ERROR: Please provide a generic distance/orientation potential file for input'
                exit(1)

        if not os.path.isfile(inputFile):
                print 'ERROR: the input potential file does not exist: ', inputFile
                exit(1)

	if querySeqFile is not None and os.path.isfile(querySeqFile):
		querySeq = LoadFASTAFile(querySeqFile)

	if not os.path.isdir(savefolder):
		os.mkdir(savefolder)

	## load up the potential file
	with open(inputFile, 'r') as fh:
		potData = cPickle.load(fh)

	if querySeq is not None and querySeq != potData[1]:
		print 'ERROR: inconsistent sequences in', querySeqFile, inputFile 
		exit(1)

	allConstraints = GenerateSplinePotential(potData, labelNames=labelNames, topRatio4Dist=topRatio4Dist, topRatio4Ori=topRatio4Ori, minSeqSep4Dist=seqSep4Dist, minSeqSep4Ori=seqSep4Ori, distPotThreshold=distPotThreshold, oriPotThreshold=oriPotThreshold)

	## save the Rosetta constraints into files
	target = os.path.basename(inputFile).split('.')[0]
	rosettaPotentialFileName = os.path.join(savefolder, target + '.pairPotential4Rosetta.SPLINE.txt')
	savefolder4histfile = os.path.join(savefolder, 'SplinePotential4' + target + '/')
	WriteSplineConstraints(allConstraints, savefile=rosettaPotentialFileName, savefolder4histfile=savefolder4histfile)

if __name__ == "__main__":
        main(sys.argv[1:])

