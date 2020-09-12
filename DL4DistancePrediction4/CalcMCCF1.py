import numpy as np
import sys
import os

from ContactUtils import CalcMCCF1
from ContactUtils import LoadContactMatrix

from utilsNoT import str_display

"""
def str_display(ls):
        if not isinstance(ls, (list, tuple, np.ndarray)):
                str_ls = '{0:.4f}'.format(ls)
                return str_ls

        str_ls = ['{0:.4f}'.format(v) for v in ls ]
        str_ls2 = ' '.join(str_ls)
        return str_ls2
"""

if __name__ == "__main__":

	if len(sys.argv) != 4:
    		print 'python CalcMCCF1.py pred_matrix_file distcb_matrix_file target'
		print '      Both matrix files are text format with L lines and each line has L columns where L is the protein sequence length'
    		exit(1)

	predFile = sys.argv[1]
	distcbFile = sys.argv[2]
	target = sys.argv[3]

	pred = LoadContactMatrix(predFile)
	truth = LoadContactMatrix(distcbFile)

	for prob in np.arange(20, 60, 2):
		#print "prob=", prob
		accs = CalcMCCF1(pred=pred, truth=truth, probCutoff=prob/100.)
		resultStr = target + ' ' + str(pred.shape[0]) + ' cutoff=' + str(prob) + ' ' + str_display(accs)
		print resultStr

