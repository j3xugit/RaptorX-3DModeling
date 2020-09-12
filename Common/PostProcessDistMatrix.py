import os
import sys
import cPickle

from PDBUtils import PostProcessDistMatrix

def Usage():
	print "python PostProcessDistMatrix.py groundTruth or TPLPKL_file [savefolder]"
	print "	the input file shall end with .native.pkl or .tpl.pkl"
	print "	savefolder: the folder for result save, default current work directory"
	print "		the result file has the same name as the input file, so if you do not want to overwrite the input file, "
	print "		please make sure that savefolder is different from the folder of the input file"

if len(sys.argv)< 2:
	Usage()
	exit(1)

infile = sys.argv[1]
savefolder = os.getcwd()

if len(sys.argv) > 2:
	savefolder = sys.argv[2]
	
if not os.path.isdir(savefolder):
	os.mkdir(savefolder)

with open(infile, 'rb') as fh:
	protein = cPickle.load(fh)

if protein.has_key('atomDistMatrix'):
	protein['atomDistMatrix'] = PostProcessDistMatrix(protein['sequence'], protein['atomDistMatrix'])

savefile = os.path.join(savefolder, os.path.basename(infile) )
with open(savefile, 'wb') as fh:
	cPickle.dump(protein, fh, protocol=cPickle.HIGHEST_PROTOCOL)
