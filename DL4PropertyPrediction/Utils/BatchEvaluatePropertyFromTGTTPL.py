import numpy as np
import sys
import os
import cPickle

from DL4PropertyPrediction import PropertyUtils
from Common.LoadTPLTGT import load_tpl as LoadTPL
from Common.LoadTPLTGT import load_tgt as LoadTGT

## this script evaluates the predicted secondary structure in TGT file using information in TPL file as ground truth

def Usage():
	print 'python BatchEvaluatePropertyFromTGTTPL.py proteinList tgt_dir tpl_dir'
	print '	This script evaluates the accuracy of the predicted secondary structure in tgt/tgt.pkl files using ground truth in .tpl.pkl or .native.pkl files'
	print '	tgt_dir: the folder for tgt or tgt.pkl files'
	print '	tpl_dir: the folder for tpl.pkl files or .native.pkl files'

def main(argv):

	if len(argv) < 3:
		Usage()
		exit(1)

	proteinListFile = argv[0]
	tgtDir = argv[1]
	tplDir = argv[2]

	with open(proteinListFile, 'r') as fh:
		names = [ line.strip() for line in list(fh) ]

	errors = []
	for name in names:
		tgtFile = os.path.join(tgtDir, name + '.tgt')
		if not os.path.isfile(tgtFile):
			tgtFile = os.path.join(tgtDir, name + '.tgt.pkl')

		tplFile = os.path.join(tplDir, name + '.tpl.pkl')
		if not os.path.isfile(tplFile):
			tplFile = os.path.join(tplDir, name + '.native.pkl')

		tgt = LoadTGT(tgtFile)
		with open(tplFile, 'r') as fh:
			tpl = cPickle.load(fh)

		tgtSS = list(tgt['SSEseq'])
		tplSS = [PropertyUtils.SS8Letter2SS3Letter[c] for c in tpl['SS_str'] ] 
		tplMissing = tpl['missing']

		numResidues = sum( [ m==0 for m in tplMissing ] )
                totalError = sum( [ p!=t for p, t, m in zip(tgtSS, tplSS, tplMissing) if m==0 ] )
                tmpError = np.array([numResidues, totalError])

		errors.append(tmpError)

	err = np.array(errors)
	err_avg = np.average(err, axis=0)
        err2 = err_avg[1:]*1./err_avg[0]

        ind_err = np.divide(err[:, 1:]*1.0, err[:, 0:1])
        err1 = np.average( ind_err, axis=0)

        print 'avg by target: ', err1, ' avg by residue: ', err2

if __name__ == "__main__":
        main(sys.argv[1:])
