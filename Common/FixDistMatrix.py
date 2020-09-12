import os
import sys
import cPickle

import PDBUtils
from PDBUtils import CalcDistMatrix, PostProcessDistMatrix

def Usage():
	print "python FixDistMatrix.py groundTruth or TPLPKL_file pdbfile [savefolder]"
	print "	the input file shall end with .native.pkl or .tpl.pkl"
	print "	pdbfile: the structure file"
	print "	savefolder: the folder for result save, default current work directory"
	print "		the result file has the same name as the input file, so if you do not want to overwrite the input file, "
	print "		please make sure that savefolder is different from the folder of the input file"

if len(sys.argv)< 3:
	Usage()
	exit(1)

infile = sys.argv[1]
pdbfile = sys.argv[2]

savefolder = os.getcwd()
if len(sys.argv) > 3:
	savefolder = sys.argv[3]
if not os.path.isdir(savefolder):
	os.mkdir(savefolder)

with open(infile, 'rb') as fh:
	protein = cPickle.load(fh)

sequence = protein['sequence']
result, pdbseq, numMisMatches, numMatches = PDBUtils.ExtractCoordinatesNDSSPBySeq(sequence, pdbfile)
if numMisMatches > 5:
	print 'ERROR: too many mismatches between query sequence and ATOM record in ', pdbfile
        exit(1)

if numMatches < min(30, 0.5*len(sequence)):
        print 'ERROR: more than half of query sequence not covered by ATOM record in ', pdbfile
        exit(1)

coordInfo, dssp = result
coordinates, numInvalidAtoms = coordInfo
if numInvalidAtoms.has_key('CA') and numInvalidAtoms['CA']>10:
	print 'ERROR: too many Ca atoms do not have valid 3D coordinates in ', pdbfile
        exit(1)
if numInvalidAtoms.has_key('CB') and numInvalidAtoms['CB']>10:
        print 'ERROR: too many Cb atoms do not have valid 3D coordinates in ', pdbfile
        exit(1)

distMatrix = CalcDistMatrix(coordinates)
protein['atomDistMatrix'] = PostProcessDistMatrix(sequence, distMatrix)

savefile = os.path.join(savefolder, os.path.basename(infile) )
with open(savefile, 'wb') as fh:
	cPickle.dump(protein, fh, protocol=cPickle.HIGHEST_PROTOCOL)
