import os
import sys
import numpy as np
import getopt

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta import *
from FoldNRelax import Relax

def Usage():
	print 'python Relax.py inputFile cstFile [ -d weight4dihedral ] [ -w weight4distance ] [ -a weight4angle] [-p weight4phipsi] [ -s w4beta] [ -n ncycles] [-t tolerance][-v savefolder]  '
	print '       inputFile can be a file ending with .pdb'
	print '       cstFile is the constraint file including all kinds of constraints acceptable by Rosetta'
	print '       -w: weight for pairwise distance-based constraints, default 1.0'
	print '       -d: weight for dihedral of orientation, default 1.0'
	print '	      -a: weight for angle of orientation, default 1.0'
	print '       -p: weight for backbone torsion angle-based constraints, default 5.0'
	print '       -s: weight for beta terms in the rosetta score, default 2.0'
	print '       -n: the number of iterations for energy minimization (default 5), not used now'
	print '	      -t: the tolerance for termination of the optimization procedure (default 0.0001, i.e., 0.01%), not used now'
	print '	      -v: the folder for result save, default current work directory'

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

	try:
                opts, args = getopt.getopt(argv[2:],"d:a:p:w:n:t:s:v:",["w4dihedral=", "w4angle=", "w4phipsi=", "w4distance=", "ncycles=", "tolerance=", "w4beta=", "savefolder=" ])
                print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

	w4phipsi = 5.
	w4distance = 1.0
	w4beta = 2.
	w4dihedral = 1.
	w4angle = 1.

	ncycles = 5
	tolerance = 0.0001

	savefolder = './'

        for opt, arg in opts:
                if opt in ("-p", "--w4phipsi"):
                        w4phipsi = np.float32(arg)
		elif opt in ("-d", "--w4dihedral"):
			w4dihedral = np.float32(arg)
		elif opt in ("-a", "--w4angle"):
			w4angle = np.float32(arg)
                elif opt in ("-w", "--w4distance"):
                        w4distance = np.float32(arg)
		elif opt in ("-s", "--w4beta"):
			w4beta = np.float32(arg)

		elif opt in ("-n", "--ncycles"):
			ncycles = np.int32(arg)
		elif opt in ("-t", "--tolerance"):
			tolerance = np.float32(arg)
			assert tolerance>0
		elif opt in ("-v", "--savefolder"):
			savefolder = arg
                else:
                        Usage()
                        exit(1)

	if w4dihedral <= 0:
		w4dihedral = w4phipsi

	targetNames = os.path.basename(inputFile).split('.')[:-1]

	init()
	rosetta.basic.options.set_boolean_option( 'run:nblist_autoupdate', True )
	pose = InitializePose(inputFile)
	pose = Relax(pose, cstFile, ncycles=ncycles, w4phipsi=w4phipsi, w4distance=w4distance, w4dihedral=w4dihedral, w4angle=w4angle, w4beta=w4beta)

	finalModelFile = '.'.join( targetNames + [ 'relaxed', str(os.getpid()), 'pdb' ] )
	finalModelFile = os.path.join(savefolder, finalModelFile)
	pose.dump_pdb(finalModelFile)

if __name__ == "__main__":
        main(sys.argv[1:])

