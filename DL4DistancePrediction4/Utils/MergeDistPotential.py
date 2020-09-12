import numpy as np
import cPickle
import sys
import os
import getopt

import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import config
from SelectAtoms import SelectAtomPair

def Usage():
    	print 'python MergeDistPotential.py [-t] distPotential_PKL1 distPotential_PKL2...'
	print '	 This script calculates the average potential of a set of potential files'
    	print '  -t: if specified, output the distance potential matrix in text format; otherwise in PKL file (default)'

def str_display(ls):
        if not isinstance(ls, (list, tuple, np.ndarray)):
                str_ls = '{0:.3f}'.format(ls)
                return str_ls

        str_ls = ['{0:.3f}'.format(v) for v in ls ]
        str_ls2 = ' '.join(str_ls)
        return str_ls2


def LoadDistPotential(file):
	fh = open(file, 'rb')
	pot = cPickle.load(fh)
	fh.close()
	return pot

def MergePotentials(inputFiles):

	targetName = None
	targetSeq = None
	targetDistCutoff = None
	potentials = dict()

	for inputFile in inputFiles:
    		if not os.path.isfile(inputFile):
			print 'The input file does not exist: ', inputFile
			exit(1)
    		content = LoadDistPotential(inputFile)
		assert len(content) >=4

    		name, sequence, potential, distCutoffs = content[:4]

		if targetName is None:
			targetName = name
			targetSeq = sequence
			targetDistCutoff = distCutoffs

			for k, v in potential.iteritems():
				potentials[k] = [v]

		else:
			## check consistency
			if sequence != targetSeq:
				print 'ERROR: inconsistent sequence in merging distance-based potential from file: ', inputFile
				print 'target seq: ', targetSeq
				print 'subjct seq: ', sequence
				exit(1)

			for k, v in potential.iteritems():
				if not potentials.has_key(k):
					potentials[k] = [v]
					targetDistCutoff[k] = distCutoffs[k]
				else:
					if np.any( targetDistCutoff[k] != distCutoffs[k] ):		
						print 'ERROR: inconsistent distance discretization in merging distance-based potential from file: ', inputFile
						print 'target distance cutoffs: ', targetDistCutoff
						print 'subjct distance cutoffs: ', distCutoffs
						exit(1)
					potentials[k].append(v)

	## now calculate average
	avgPotential = dict()
	for k, vlist in potentials.iteritems():
		avgPotential[k] = np.average(vlist, axis=0)
		
	return (targetName, targetSeq, avgPotential, targetDistCutoff)


def main(argv):

    	inputFiles = None
	saveInTxt = False

	if len(argv) < 2:
		Usage()
		exit(1)

    	try:
        	opts, args = getopt.getopt(argv,"t",["textFormat="])
        	print opts, args

    	except getopt.GetoptError:
        	Usage()
        	exit(1)


    	if len(args) < 2:
        	Usage()
        	exit(1)

	inputFiles = args

    	for opt, arg in opts:

		if opt in ("-t", "--textFormat"):
	    		saveInTxt = True

		else:
	    		Usage()
	    		exit(1)

	avgPotential = MergePotentials(inputFiles)


	savefile = avgPotential[0] + '.mergedDistPotential.pkl'
	fh = open(savefile, 'wb')
	cPickle.dump(avgPotential, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	fh.close()

	if not saveInTxt:
		return

	targetName, sequence, potential, distCutoff = avgPotential[:4]

	## save to text file
	potentialFileName = targetName + '.mergedDistPotential.txt'
	fh = open(potentialFileName, 'w')
	fh.write('#TARGET\t' + targetName + '\n')
	fh.write('#SEQ\t' + sequence + '\n')
	fh.write('#Distance cutoffs\t' + str_display(distCutoff) + '\n')

	for response, pot in potential.iteritems():
		apt = config.Response2LabelName(response)
                distLabelType = config.Response2LabelType(response)

		size = pot.shape
		for i in xrange(size[0]):
			rawPotStrs = []

			for j in xrange(i+ minSeqSep, size[1]):
				atom1, atom2 = SelectAtomPair(sequence, i, j, apt)
				y = pot[i, j]

				rawPotStr = ' '.join(['AtomPair', atom1.upper(), str(i+1), atom2.upper(), str(j+1), distLabelType] + [ "{:.4f}".format(e) for e in y ] )
				rawPotStrs.append(rawPotStr)


			if len(rawPotStrs) >0:
				fh.write('\n'.join(rawPotStrs) + '\n')

				
	fh.close()


if __name__ == "__main__":
    	main(sys.argv[1:])
