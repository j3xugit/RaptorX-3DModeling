import numpy as np
import sys
import os

from ContactUtils import TopAccuracy
from ContactUtils import LoadContactMatrix

if __name__ == "__main__":

	if len(sys.argv) < 3:
    		print 'python EvaluateContactAccuracyTXT.py predContactMatrixFile_txt nativeContactMatrixFile_txt [targetName]'
		print '	This script evaluates the accuracy of a predicted contact matrix by comparing it with native contact matrix'
		print '	Both matrix files are in text format with L lines and each line has L columns where L is the protein sequence length'
    		exit(1)

	predFile = sys.argv[1]
	nativeFile = sys.argv[2]

	if len(sys.argv)>3:
		target = sys.argv[3]
	else:
		target = os.path.basename(nativeFile).split('.')[0]

	pred = LoadContactMatrix(predFile)
	truth = LoadContactMatrix(nativeFile)

	accs = TopAccuracy(pred, truth)
	accsStr = [ str(a) for a in accs ]
	resultStr = target + ' ' + str(pred.shape[0]) + ' TopAcc '
	resultStr += (' '.join(accsStr) )
	print resultStr

