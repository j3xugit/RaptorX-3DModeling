import os
import sys
import cPickle

from LoadTPLTGT import load_tgt as LoadTGT
from LoadTPLTGT import load_tpl as LoadTPL
from LoadHHM import load_hhm as LoadHHM

## this script checks if the following types of files have Xs in sequence: .hhm, .hhm.pkl, .tgt, .tgt.pkl, .tpl, .tpl.pkl

def Usage():
	print 'Usage: python CheckIfSequenceHasXs.py inputFile'
	print '	this script checks if the following types of files have Xs in sequence: .hhm, .hhm.pkl, .tgt, .tgt.pkl, .tpl, .tpl.pkl'
	print '	if Xs are contained in the sequence, the file name will be printed out as well as the sequence'

if len(sys.argv) < 2:
	Usage()
	exit(1)

infile = sys.argv[1]

if infile.endswith('.pkl'):
	with open(infile, 'rb') as fh:
		seq = cPickle.load(fh)['sequence']
elif infile.endswith('.hhm'):
	hhm = LoadHHM(infile)
	seq = hhm['sequence']
elif infile.endswith('.tgt'):
	tgt = LoadTGT(infile)
	seq = tgt['sequence']
elif infile.endswith('.tpl'):
	tpl = LoadTPL(infile)
	seq = tpl['sequence']
else:
	Usage()
	exit(1)

if 'X' in seq:
	print infile, seq
