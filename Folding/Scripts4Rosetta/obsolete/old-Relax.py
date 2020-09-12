import os
import sys
import numpy as np
import getopt

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta import *
from FoldNRelax import Relax, RelaxWithoutConstraints, CheckCSTFile

def Usage():
	print 'python Relax.py inputFile cstFile [-r] [ -d weight4dihedral ] [ -w weight4distance ] [ -a weight4angle] [-e weight4beta] [ -n nrepeats] [-s savefolder]  '
	print '       inputFile: a file with name ending with .pdb'
	print '       cstFile: the constraint file including all kinds of constraints acceptable by Rosetta'
	print '	      -r: warm up by relaxing without constraints, default True'
	print '       -w: weight for distance potential, default 0.2'
	print '       -d: weight for dihedral of orientation, default 0.2'
	print '	      -a: weight for angle of orientation, default 0.2'
	print '	      -e: weight for beta formation, default None'
	print '       -n: the number of repeats for energy minimization, default 5'
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
	cstFile = argv[1]

	if not os.path.isfile(inputFile):
                print 'ERROR: invalid input file for relax:', inputFile
                exit(1)
        if not os.path.isfile(cstFile):
                print 'ERROR: invalid constraint file for relax:', cstFile
                exit(1)

	try:
                opts, args = getopt.getopt(argv[2:],"rd:a:w:e:n:s:",["warmup", "w4dihedral=", "w4angle=", "w4distance=", "w4beta=", "ncycles=", "savefolder=" ])
                #print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

	warmup=True
	w4distance = 0.2
	w4dihedral = 0.2
	w4angle = 0.2
	w4beta = None

	ncycles = 5

	savefolder = os.getcwd()

        for opt, arg in opts:
		if opt in ("-r", "--warmup"):
			warmup=True
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

		elif opt in ("-s", "--savefolder"):
			savefolder = arg
                else:
                        Usage()
                        exit(1)

	if not inputFile.endswith('.pdb'):
		print 'ERROR: the input file for Relax shall have a name ending with .pdb'
		exit(1)

	if not os.path.isfile(cstFile):
                print 'ERROR: invalid constraint file: ', cstFile
                exit(1)

        CheckCSTFile(cstFile)

	init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 200')
	rosetta.basic.options.set_boolean_option( 'run:nblist_autoupdate', True )
	pose = InitializePose(inputFile)
	if warmup:
		pose = RelaxWithoutConstraints(pose)

        constraints = protocols.constraint_movers.ConstraintSetMover()
        constraints.constraint_file(cstFile)
        constraints.add_constraints(True)
	pose = Relax(pose, constraints, ncycles=ncycles, w4distance=w4distance, w4dihedral=w4dihedral, w4angle=w4angle, w4beta=w4beta)

	targetNames = os.path.basename(inputFile).split('.')[:-1]
	finalModelFile = '.'.join( targetNames + [ 'relaxed', str(os.getpid()), 'pdb' ] )
	finalModelFile = os.path.join(savefolder, finalModelFile)
	pose.dump_pdb(finalModelFile)

if __name__ == "__main__":
        main(sys.argv[1:])

