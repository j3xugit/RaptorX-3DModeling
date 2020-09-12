import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np

import config
import DistanceUtils

from utilsNoT import str_display

import getopt

def Usage():

    	print 'python BatchEvaluateProbAccuracy.py poteinList predProbMatrix_folder ground_truth_folder [ -s minSeqSep]  [-c probCutoff for prediction] [ -d distCutoff for native]'
	print '  This script evaluate predicted dist prob accuracy for a list of proteins in their predicted dist matrix files '
    	print '  predProbMatrix_folder: a folder containing predicted distance prob with name like XXX.predictedDistMatrix.pkl'
    	print '     A predicted distance prob matrix file contains a tuple of 6 or 7 items: name, primary sequence, predDistProbMatrix, ...'
	print '     It is a dict() and each item is a matrix with dimension L*L*numDistBins'
	print '	 -s: optional. The minimum sequence separation between two residues for which its distance is evaluated. default 12'
	print '  -c: specify the prob cutoff for prediction (default 0.5)'
        print '  -d: specify the dist cutoff for native (default 15.0)'
	print '  This script will output precision, recall, and F1' 

"""
def str_display(ls):
        if not isinstance(ls, (list, tuple, np.ndarray)):
                str_ls = '{0:.4f}'.format(ls)
                return str_ls

        str_ls = ['{0:.4f}'.format(v) for v in ls ]
        str_ls2 = ' '.join(str_ls)
        return str_ls2
"""

def main(argv):


	if len(argv)<3:
		Usage()
		exit(1)

	proteinListFile = argv[0]
	predFolder = argv[1]
	nativefolder = argv[2]
	fileSuffix = '.predictedDistMatrix.pkl'

	minSeqSep = 12
	probCutoff4Pred = 0.5
        distCutoff4Native = 15.0

	try:
                opts, args = getopt.getopt(argv[3:], "s:c:d:", ["minSeqSep=", "probCutoff4Pred=", "distCutoff4Native="])
                #print opts, args
        except getopt.GetoptError as err:
                print err
                Usage()
                exit(1)


        for opt, arg in opts:
                if opt in ("-s", "--minSeqSep"):
                        minSeqSep = np.int32(arg)
                        if minSeqSep < 2:
                                print 'ERROR: minSeqSep too small'
                                exit(1)

                elif opt in ("-c", "--probCutoff4Pred"):
                        probCutoff4Pred = np.float32(arg)
                        if probCutoff4Pred < 0.1 :
                                print 'ERROR: prob cutoff for predicted is too small '
                                exit(1)

                elif opt in ("-d", "--distCutoff4Native"):
                        distCutoff4Native = np.float32(arg)
                        if distCutoff4Native < 1.0 :
                                print 'ERROR: dist cutoff for native is too small '
                                exit(1)
                else:
                        Usage()
                        exit(1)

	if not os.path.isfile(proteinListFile):
		print 'the protein list file does not exist: ', proteinListFile
		exit(1)

	if not os.path.isdir(predFolder):
		print 'the folder for predicted bound matrix files does not exist: ', predFolder
		exit(1)

	if not os.path.isdir(nativefolder):
		print 'the folder for native distance matrix files does not exist: ', nativefolder
		exit(1)

	fh = open(proteinListFile, 'r')
	proteins = [ line.strip() for line in list(fh) ]
	fh.close()

	AccPerProtein = dict()
	accs = dict()
	for protein in proteins:
		predFile = os.path.join( predFolder, protein + fileSuffix ) 
		if not os.path.isfile(predFile):
			print 'the distance bound file does not exist: ', predFile
			exit(1)

		fh = open(predFile, 'rb')
		pred = cPickle.load(fh)[2]
		fh.close()

		nativeFile = os.path.join(nativefolder, protein + '.atomDistMatrix.pkl' )
		if not os.path.isfile(nativeFile):
			print 'the native atomDistMatrix file does not exist: ', nativeFile
			exit(1)

		fh = open(nativeFile, 'rb')
		native = cPickle.load(fh)
		fh.close()		

                acc = DistanceUtils.EvaluateProbAccuracy(pred, native, minSeqSep=minSeqSep, probCutoff=probCutoff4Pred, distCutoff4Native=distCutoff4Native)
		AccPerProtein[protein] = acc

		for k, v in acc.iteritems():
			if not accs.has_key(k):
				accs[k] = [ v ]
			else:
				accs[k].append(v)

	for k, v in accs.iteritems():
		accs[k] = np.average(v, axis=0)
		print 'average', k, str_display(accs[k])


	for k, v in AccPerProtein.iteritems():
		for apt, value in v.iteritems():
			print k, apt, str_display(value)


if __name__ == "__main__":
    	main(sys.argv[1:])
