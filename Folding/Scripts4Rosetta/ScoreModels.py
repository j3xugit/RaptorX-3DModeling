import os
import sys
import numpy as np
import getopt
import glob
import tempfile
import shutil

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta import *

from ScoreOneModel import Score as ScoreOneModel
from FoldNRelax import DeriveRosettaCSTFile, CheckCSTFile

def Usage():
	print 'python ScoreModels.py inputFolder cstFile [ -I extraInput | -s savefile]'
	print '	inputFolder: a folder for all PDB files to be scored'
	print '	cstFile: a Rosetta constraint file or a file for predicted distance/orientation probability in PKL format'
	print ' -I: when this option specified, cstFile shall be interpreted as the predicted dist/orientation matrix file instead of Rosetta constraint file'
        print '           extraInput is used to inform how to derive Rosetta constrains from predicted dist/ori matrix, e.g., potential type and some hyperparamters'
        print '           It also provides predicted Phi/Psi angle file, e.g., -I T1111.predictedProperties.pkl,alpha=1.61 or -I T111.predictedProperties.pkl'
	print '	savefile: write reslts to this file. If not specified, then write to screen'

def main(argv):

	if len(argv) < 2:
		Usage()
		exit(1)

	inputFolder = argv[0]
	pairFile = argv[1]

	pairFileIsCST = True
        PhiPsiFile = None
        param4Potential = 1.61

	savefile = None

	if not os.path.isdir(inputFolder):
		print 'ERROR: invalid folder for models: ', inputFolder
		exit(1)

	if not os.path.isfile(pairFile):
                print 'ERROR: invalid file for predicted distance/orienation information:', pairFile
                exit(1)

	try:
                opts, args = getopt.getopt(argv[2:],"I:s:",["extraInput=", "savefile=" ])
                #print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

	for opt, arg in opts:
                if opt in ("-I", "--extraInput" ):
                        pairFileIsCST = False

                        fields = arg.split(',')
                        for field in fields:
                                columns = field.split('=')
                                assert len(columns) > 0
                                if len(columns) == 1:
                                        if PhiPsiFile is None:
                                                PhiPsiFile = columns[0]
                                                continue
                                        else:
                                                print 'ERROR: wrong format in the -I option'
                                                exit(1)

                                if columns[0].upper() == 'phipsi'.upper():
                                        PhiPsiFile = columns[1]

                                elif columns[0].upper() == 'alpha'.upper():
                                        param4Potential = np.float32(columns[1])
                                else:
                                        print 'ERROR: unsupported arguments for -I option:', arg
                                        exit(1)

                elif opt in ("-s", "--savefile"):
                        savefile = arg
                else:
                        Usage()
                        exit(1)

	init('-out:level 100')
	#init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 100')
	rosetta.basic.options.set_boolean_option( 'run:nblist_autoupdate', True )

	if not pairFileIsCST:
                assert PhiPsiFile is not None
                if not os.path.isfile(PhiPsiFile):
                        print 'ERROR: invalid file for predicted Phi/Psi angles: ', PhiPsiFile
                        exit(1)

                if param4Potential > 10:
                        param4Potential = random.uniform(1.57, 1.63)
                assert param4Potential <= 1.63 and (param4Potential >= 1.57), 'param for distance potential is out of range'
                print 'alpha for DFIRE is', param4Potential

                ## create a CST folder and generate the CST file
                if os.path.isdir('/dev/shm'):
                        cstfolder = tempfile.mkdtemp(prefix='cstDir4'+os.path.basename(pairFile).split('.')[0], dir='/dev/shm')
                else:
                        cstfolder = tempfile.mkdtemp(prefix='cstDir4'+os.path.basename(pairFile).split('.')[0])

                cstFile = DeriveRosettaCSTFile(pairFile, PhiPsiFile, saveFolder=cstfolder, param4Potential=param4Potential)
        else:
                cstFile = pairFile

        if not CheckCSTFile(cstFile):
                print 'ERROR: incorrect cstFile: ', cstFile
                if not pairFileIsCST:
                        shutil.rmtree(cstfolder)
                exit(1)

        constraints = protocols.constraint_movers.ConstraintSetMover()
        constraints.add_constraints(True)
        constraints.constraint_file(cstFile)

        scorefxn = create_score_function('ref2015')
        scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, 1)
        scorefxn.set_weight(rosetta.core.scoring.dihedral_constraint, 1)
        scorefxn.set_weight(rosetta.core.scoring.angle_constraint, 1)

	namePattern = inputFolder + '/*.pdb' 
	files = glob.glob(namePattern)

	scores = []
	for f in files:
		pose = pose_from_pdb(f)
		score = ScoreOneModel(pose, scorefxn, constraints)
		score['modelFile'] = f
		scores.append(score)

	## remove cstfolder if needed
	if not pairFileIsCST:
                shutil.rmtree(cstfolder)

	scores = sorted(scores, key=lambda x: x['totalPot'] )

	outStrs = []
	for s in scores:
		overall = [ s['totalPot'], s['totalDistPot'], s['totalDihedralPot'], s['totalAnglePot'] ]
		average = [ e/len(s['distPots']) for e in overall ]
		average = [ '{:.2f}'.format(e) for e in average ]
		rawOut = [ os.path.basename(s['modelFile']) ] + average 
		rawOutStr = ' '.join([ str(e) for e in rawOut ])
		outStrs.append(rawOutStr)
	finalOutStr = '\n'.join(outStrs)
	
	if savefile is None:
		print finalOutStr
		return

	with open(savefile, 'w') as fh:
		fh.writelines(finalOutStr)

if __name__ == "__main__":
        main(sys.argv[1:])

