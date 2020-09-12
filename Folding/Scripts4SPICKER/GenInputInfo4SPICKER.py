import os
import sys
import glob
import getopt
import numpy as np
import json

from Common.SequenceUtils import LoadFASTAFile, AA1LetterCode23LetterCode
from Common.PDBUtils import ExtractCoordinatesBySeq

def Usage():
	print 'python GenInputInfo4SPICKER.py [-l ] [-s savefolder] [-c cutoff-method] seqFile [modelListFile / modelFolder]'
	print '	seqFile: the primary sequence file in FASTA format'
	print '	modelListFile: a file for a list of models to be clustered'
	print '	modelFolder: the folder of decoys to be clustered by spicker'
	print '		only one of modelListFile and modelFolder shall be provided'
	print '	-l: if specified, the input modelListFile; otherwise modelFolder'
	print '	-c: cutoff method for SPICKER(see readme in SPICKER for explanation), default -1'
	print '	-s: the folder for results, default current work directory, but better provide a folder'

## extract Ca coordinates and save them into a file, which be returned
def ExtractCaCoordinates(sequence, decoySet, savefile):

	fh = open(savefile, 'w')

	for decoyFile in decoySet:
		coords, pdbseq, numMisMatches, numMatches= ExtractCoordinatesBySeq(sequence, decoyFile, atoms=['CA'])
		if coords is None:
			print 'ERROR: failed to extract Ca coordinates from ', decoyFile
			exit(1)

		if numMisMatches > 0:
			print 'ERROR: there is at least one mismatch between sequence and decoy file: ', decoyFile
			exit(1)

		atoms, numInvalidAtoms = coords
		if numInvalidAtoms['CA'] >0:
			print 'ERROR: there is at least one invalid Ca atom in decoy file: ', decoyFile
			exit(1)

		fh.write('     %d%9.3f      1      1\n' % (len(atoms), 100))
		## print out the Ca coordinates
		for atom in atoms:
			x, y, z = atom['CA']
			fh.write('%9.3f%9.3f%9.3f\n' % (x, y, z))

	fh.close()
		
	
if len(sys.argv)< 3:
	Usage()
	exit(1)

try:
	opts, args = getopt.getopt(sys.argv[1:],"lc:s:",["modelList=", "cutoff=", "savefolder="])
except getopt.GetoptError:
        Usage()
        exit(1)

if len(args) != 2:
        Usage()
        exit(1)

seqFile = args[0]
target = os.path.basename(seqFile).split('.')[0]
sequence = LoadFASTAFile(seqFile)

modelListOrFolder = args[1]
if not (os.path.isdir(modelListOrFolder) or os.path.isfile(modelListOrFolder) ):
	print 'ERROR: invalid model list file or folder:', modelListOrFolder
	exit(1)

InputIsModelList = False

savefolder = os.getcwd()
cutoff=-1
cutoffs=[1, -1, -2]

for opt, arg in opts:
	if opt in ("-l", "--modelList"):
		InputIsModelList = True

        elif opt in ("-s", "--savefolder"):
		savefolder = arg
		if not os.path.isdir(savefolder):
			os.mkdir(savefolder)

        elif opt in ("-c", "--cutoff"):
		cutoff = np.int32(arg)
		if cutoff not in set(cutoffs):
			print 'ERROR: cutoff shall have value in ', cutoffs
			exit(1)
        else:
                Usage()
                exit(1)

if InputIsModelList:
	with open(modelListOrFolder, 'r') as fh:
		decoys = [ line.strip() for line in list(fh) ]

else:
	namePattern=modelListOrFolder + '/' + target + '*.pdb'
	decoys = glob.glob(namePattern)

if len(decoys) < 1:
	print 'ERROR: there are no decoys for clustering in', modelListOrFolder
	exit(1)

if len(decoys) < 10:
	print 'WARNING: there are very few decoys for clustering in', modelListOrFolder

## generate rmsinp
rmsinp = os.path.join(savefolder, 'rmsinp')
with open(rmsinp, 'w') as fh:
	fh.write('%d %d\n' % (1, len(sequence)) )
	fh.write('%d\n' % len(sequence) )

## generate seq.dat
seqdat = os.path.join(savefolder, 'seq.dat')
with open(seqdat, 'w') as fh:
	for i, aa in zip(range(1, 1+len(sequence)), sequence):
		fh.write('%5d   %s\n' % (i, AA1LetterCode23LetterCode[aa]) )

## extract Ca coodinates from each decoy file
mapping_CA2decoy = dict()
#mapping_decoy2CA = dict()

savefiles = []
numDecoysPerFile=100
for i in range(0, len(decoys), numDecoysPerFile):
	if i+numDecoysPerFile <= len(decoys):
		decoySet = decoys[i:i+numDecoysPerFile]
	else:
		decoySet = decoys[i:]
	savefile = os.path.join(savefolder, 'MSet' + str(i) + '.CA')
	ExtractCaCoordinates(sequence, decoySet, savefile)
	savefiles.append( os.path.basename(savefile) )
	
	for j, decoy in zip(range(len(decoySet)), decoySet):
		mapping_CA2decoy[ str(i) + '-' + str(j) ] = decoy
		#mapping_decoy2CA[ decoy ] = (i, j)

## generate tra.in
train= os.path.join(savefolder, 'tra.in')
with open(train, 'w') as fh:
	fh.write('%d %d %d\n' % (len(savefiles), cutoff, 1) )
	fh.writelines('\n'.join(savefiles) +'\n' )


mappingfile = os.path.join(savefolder, 'mapping.json')
with open(mappingfile, 'w') as fh:
	json.dump( (numDecoysPerFile, mapping_CA2decoy), fh )
