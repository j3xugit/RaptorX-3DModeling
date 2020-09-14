import os
import sys
import numpy as np
import getopt
import tempfile
import shutil
import random
import socket

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta import *

from FoldNRelax import Relax, RelaxWithoutConstraints, CheckCSTFile, DeriveRosettaCSTFile, ReadFASTAFile

def Usage():
	print 'python Relax.py inputFile cstFile [-r ncycles] [-I extraInput ] [ -d weight4dihedral ] [ -w weight4distance ] [ -a weight4angle] [-e weight4beta] [ -n nrepeats] [-s savefolder]  '
	print '       inputFile: a file with name ending with .pdb'
	print '       cstFile: the constraint file including all kinds of constraints acceptable by Rosetta'
  	print '       -I: when this option specified, cstFile shall be interpreted as the predicted dist/orientation matrix file instead of Rosetta constraint file'
        print '		extraInput is used to inform how to derive Rosetta constrains from predicted dist/ori matrix, e.g., potential type and some hyperparamters'
        print '		It also provides predicted Phi/Psi angle file, e.g., -I T1111.predictedProperties.pkl,alpha=1.61 or -I T111.predictedProperties.pkl'
	print '	      -r: the number of cycles (<=3) for relaxing without any predicted constraints, default 0'
	print '		When it is 0, will not do relax tionwithout any predicted constraints'
	print '       -w: weight for distance potential, default 0.2'
	print '       -d: weight for dihedral of orientation, default 0.2'
	print '	      -a: weight for angle of orientation, default 0.2'
	print '	      -e: weight for beta formation, default None'
	print '       -n: the number of cycles for energy minimization with predicted constraints, default 5'
	print '	      -s: the folder for result save, default current work directory'

def InitializePose(inputFile):
	pose = pose_from_pdb(inputFile)
	switch = SwitchResidueTypeSetMover("fa_standard")
	switch.apply(pose)
	return pose

def main(argv):

	if len(argv) < 2:
		Usage()
		exit(1)

	inputFile = argv[0]
	pairFile = argv[1]

	pairFileIsCST = True
	PhiPsiFile = None
	param4Potential = 1.61

	if not os.path.isfile(inputFile):
                print 'ERROR: invalid input file for relax:', inputFile
                exit(1)
        if not os.path.isfile(pairFile):
                print 'ERROR: invalid file for predicted distance/orienation information:', pairFile
                exit(1)

	try:
                opts, args = getopt.getopt(argv[2:],"rI:d:a:w:e:n:s:",["warmup=", "extraInput=", "w4dihedral=", "w4angle=", "w4distance=", "w4beta=", "ncycles=", "savefolder=" ])
                #print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

	warmup=0
	w4distance = 0.2
	w4dihedral = 0.2
	w4angle = 0.2
	w4beta = None

	ncycles = 5

	savefolder = os.getcwd()

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

		elif opt in ("-r", "--warmup"):
			warmup=np.int32(arg)
			assert warmup>=0 and warmup<=3 

		elif opt in ("-d", "--w4dihedral"):
			w4dihedral = np.float32(arg)
		elif opt in ("-a", "--w4angle"):
			w4angle = np.float32(arg)
                elif opt in ("-w", "--w4distance"):
                        w4distance = np.float32(arg)
		elif opt in ("-e", "--w4beta"):
			w4beta = np.float32(arg)

		elif opt in ("-n", "--ncycles"):
			ncycles = np.int32(arg)
			assert ncycles>=0

		elif opt in ("-s", "--savefolder"):
			savefolder = arg
                else:
                        Usage()
                        exit(1)

	if not inputFile.endswith('.pdb'):
		print 'ERROR: the input file for Relax shall have a name ending with .pdb'
		exit(1)

	init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 100')
	rosetta.basic.options.set_boolean_option( 'run:nblist_autoupdate', True )

	pose = InitializePose(inputFile)
	if warmup > 0:
		pose = RelaxWithoutConstraints(pose, ncycles=warmup)

	if not pairFileIsCST:
                assert PhiPsiFile is not None
                if not os.path.isfile(PhiPsiFile):
                        print 'ERROR: invalid file for predicted Phi/Psi angles: ', PhiPsiFile
                        exit(1)

                if param4Potential > 10:
                        param4Potential = random.uniform(1.57, 1.63)

                print 'alpha for DFIRE:', param4Potential
                assert param4Potential <= 1.63 and (param4Potential >= 1.57), 'In Relax.py, param for DFIRE distance potential is out of range'

		seq = None
                seqLen = pose.total_residue()

                ## create a CST folder and generate the CST file
		machine = socket.gethostname().split('.')[0]
                print 'Running Relax.py on', machine

		## machines with large CPU memory
		LargeRAMs = ['raptorx6', 'raptorx7', 'raptorx8', 'raptorx9', 'raptorx10']

                if machine in LargeRAMs:
                        cstfolder = tempfile.mkdtemp(prefix='cstDir4'+os.path.basename(pairFile).split('.')[0], dir='/dev/shm')

                elif seqLen<400 and os.path.isdir('/dev/shm'):
                        cstfolder = tempfile.mkdtemp(prefix='cstDir4'+os.path.basename(pairFile).split('.')[0], dir='/dev/shm')
                else:
                        cstfolder = tempfile.mkdtemp(prefix='cstDir4'+os.path.basename(pairFile).split('.')[0])

                cstFile = DeriveRosettaCSTFile(seq, pairFile, PhiPsiFile, saveFolder=cstfolder, param4Potential=param4Potential)
        else:
                cstFile = pairFile

  	if not CheckCSTFile(cstFile):
                print 'ERROR: incorrect cstFile: ', cstFile
                if not pairFileIsCST:
                        shutil.rmtree(cstfolder)
                exit(1)

        constraints = protocols.constraint_movers.ConstraintSetMover()
        constraints.constraint_file(cstFile)
        constraints.add_constraints(True)
	constraints.apply(pose)

	## remove cstfolder after adding the constraints
	if not pairFileIsCST:
		shutil.rmtree(cstfolder)

	pose = Relax(pose, ncycles=ncycles, w4distance=w4distance, w4dihedral=w4dihedral, w4angle=w4angle, w4beta=w4beta)

	targetNames = os.path.basename(inputFile).split('.')[:-1]
	finalModelFile = '.'.join( targetNames + [ 'relaxed', str(os.getpid()), 'pdb' ] )
	finalModelFile = os.path.join(savefolder, finalModelFile)
	pose.dump_pdb(finalModelFile)

if __name__ == "__main__":
        main(sys.argv[1:])

