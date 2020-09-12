import sys
import os
import cPickle
import numpy as np

import ContactUtils

if __name__ == "__main__":

	if len(sys.argv) != 3:
    		print 'python CalcCASPContactPredAccuracy.py pred_CASP_file nativeDistMatrixFilePKL'
		print '	the native file may end with .native.pkl or .atomDistMatrix.pkl'
    		exit(1)

	predFile = sys.argv[1]
	distFile = sys.argv[2]

	pred, target, sequence = ContactUtils.LoadContactMatrixInCASPFormat(predFile)

	with open(distFile, 'rb') as fh:
		nativeMatrix = cPickle.load(fh)

	if distFile.endswith('.native.pkl'):
		native = nativeMatrix['atomDistMatrix']
	else:
		native = nativeMatrix
		
	accs = ContactUtils.EvaluateCbCbContactPrediction(pred, native)
 
	accsStr = [ str(a) for a in accs ]
	resultStr = target + ' ' + str(pred.shape[0]) + ' TopAcc '
	resultStr += (' '.join(accsStr) )
	print resultStr


