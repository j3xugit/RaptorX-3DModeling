import os
import sys
import numpy as np
import getopt

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta import *

def Usage():
	print 'python FoldNRelax.py inputFile cstFile [ -w weight4distance ] [ -d weight4dihedral ] [-a weight4angle] [-p weight4phipsi] [-s weight4beta] [ -n ncycles4initModel ] [-t tolerance] [-v savefolder ] [-r] [-b] [-q]'
	print '       inputFile can be a file ending with .fasta, .seq, and .pdb'
	print '       cstFile is the constraint file including all kinds of constraints acceptable by Rosetta'
	print '       -w: weight for pairwise distance-based constraints, default 1.0'
	print '       -d: weight for orientation dihedral, default 1.0'
	print '	      -a: weight for orientation angle, default 1.0'
	print '	      -p: weight for backbone phi/psi, default 5.0. when -d is used, the weight for phi/psi is same as other dihedral'
	print '	      -s: weight for beta strand and long-range hydrogen-bond, default 2.0. When <0, not used'
	print '       -n: the number of iterations for energy minimization (default 2000) at the initial folding stage. '
	print '	      -t: the tolerance for termination of the optimization procedure (default 0.0001, i.e., 0.01%) '
	print '	      -v: the folder for result save, default current work directory'
	print '       -r: if specified, do not sample phi/psi angles by the predicted distribution, otherwise sample'
	print '       -b: if specified, do not use neighbor list in minimization, otherwise use neighbor list'
	print '	      -q: if specified, do not run FastRelax (i.e. fast mode), otherwise run FastRelax'

def ReadFASTAFile(inputFile):
	f = open(inputFile, 'r')    # open the file
    	sequence = f.readlines()    # read the text
    	f.close()    # close it
    	# removing the trailing "\n" and any header lines
    	sequence = [line.strip() for line in sequence if not '>' in line]
    	sequence = ''.join( sequence )    # combine into a single sequence
	return sequence
	
## extract phi/psi distribution from cstFile, which has the Rosetta format with AMBER function
def ExtractPhiPsiDistribution(seq=None, cstFile=None):

	assert seq is not None
	seqLen = len(seq)	
	PhiDistribution = [None] * seqLen 
	PsiDistribution = [None] * seqLen 
	rows = None

	with open(cstFile, 'r') as f:
		content = f.readlines()
		rows = [ c for c in content if 'AMBER' in c and 'Dihedral' in c ]

	for r in rows:
		fields = r.split()
		assert len(fields)==13
		assert fields[0] == 'Dihedral'
		x0 = np.float32(fields[-3])
		n = np.float32(fields[-2])
		k = np.float32(fields[-1])
		idx = np.int32(fields[4])-1

		if fields[1]=='N' and fields[3]=='CA' and fields[5]=='C' and fields[7]=='N':
			## psi
			PsiDistribution[idx] = (x0, n, k)
		elif fields[1]=='C' and fields[3]=='N' and fields[5]=='CA' and fields[7]=='C':
			##phi
			PhiDistribution[idx] = (x0, n, k)
		else:
			print 'WARNING: unknown duhedral type in line: ', r
			continue

	distribution = dict()
	distribution['phi'] = PhiDistribution
	distribution['psi'] = PsiDistribution

	return distribution
			

## distribution is a list of tuples. Each tuple=(x0, n_period, k) defines one AMBER function for phi and/or psi
## Here an AMBER function is defined as f(x) = k * (1 + cos( ( n_period * x ) - x0 ) ) where k is the amplitude, n_period is the period and x0 is the maximum point
## x and x0 are radians and range from 0 to 2*pi
## the prob of x is proportional to exp( -f(x) )
def SampleDihedralsByAMBER(distribution):
	## generate a list of discrete angles from 0 to 360
	step = 8.
	anchors_degree = np.arange(step/2.-180, 180.-0.1, step)
	anchors = anchors_degree / 180. * np.pi
	#print anchors_degree
	#print anchors
	
	## random sample pertubations between -step/2 and step/2
	samples = np.random.uniform(-step/2, step/2, size=len(distribution)).astype(np.float32) 

	for idx, d in zip(range(len(distribution)), distribution):
		if d is None:
			samples[idx] = np.random.uniform(-180., 180.)
			continue

		x0, n, k = d
		#print idx, d
		## ll represents log-likelihood
		ll = -k * ( 1 + np.cos( (n*anchors) - x0 ) )
	
		## calculate the probability of the anchors by d
		prob = np.exp(ll)
		prob = prob / np.sum(prob)

		#print idx, prob

		## random sample one anchor by prob
		sample = np.random.choice(anchors_degree, p=prob)
		samples[idx] += sample

	return samples

def InitializePose(inputFile, PhiPsiDistribution=None):
	if inputFile.endswith('.fasta') or inputFile.endswith('.seq'):
		sequence = ReadFASTAFile(inputFile)
		pose = pose_from_sequence(sequence)

		if PhiPsiDistribution is None:
			for i in range(1, pose.total_residue() + 1):
				pose.set_phi(i, np.random.uniform(-180, 180))
				pose.set_psi(i, np.random.uniform(-180, 180))
		else:
			PhiDistribution = PhiPsiDistribution['phi']
			phis = SampleDihedralsByAMBER(PhiDistribution)
			PsiDistribution = PhiPsiDistribution['psi']
			psis = SampleDihedralsByAMBER(PsiDistribution)
		
			for i, phi, psi in zip(range(1, pose.total_residue() + 1), phis, psis):
				pose.set_phi(i, phi)
				pose.set_psi(i, psi)

	else:
		pose = pose_from_pdb(inputFile)

	return pose

def Fold(pose=None, cstFile=None, ncycles=2000, tolerance=0.0001, w4distance=1., w4angle=1., w4dihedral=1., w4phipsi=5., w4beta=2., UseNBList=True):
	assert pose is not None
	if not os.path.isfile(cstFile):
		print 'ERROR: invalid constraint file: ', cstFile
		exit(1)

	constraints = protocols.constraint_movers.ConstraintSetMover()
	constraints.constraint_file(cstFile)
	constraints.add_constraints(True)
	constraints.apply(pose)

	mm = MoveMap()
	mm.set_bb(True)
	mm.show()

	if w4beta >0:
		scorefxn =  create_score_function("score4_smooth")
		scorefxn.set_weight(rosetta.core.scoring.ss_pair, w4beta)
		scorefxn.set_weight(rosetta.core.scoring.sheet, w4beta)
		scorefxn.set_weight(rosetta.core.scoring.hbond_lr_bb, w4beta)
	else:
		scorefxn = ScoreFunction()

	scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, w4distance)
	scorefxn.set_weight(rosetta.core.scoring.dihedral_constraint, w4dihedral)
	scorefxn.set_weight(rosetta.core.scoring.angle_constraint, w4angle)

	scorefxn.show(pose)

	#minmover = protocols.minimization_packing.MinMover(mm, scorefxn, 'lbfgs_armijo_nonmonotone', 0.01, True)
	#minmover.cartesian(True)
	minmover = protocols.minimization_packing.MinMover()
	minmover.movemap(mm)
	minmover.nb_list(UseNBList)

	minmover.score_function(scorefxn)
	minmover.min_type('lbfgs_armijo_nonmonotone')

	minmover.max_iter(ncycles)
	minmover.tolerance(tolerance)	
	minmover.show()
	minmover.apply(pose)

	scorefxn.show(pose)

	switch = SwitchResidueTypeSetMover("fa_standard")
	switch.apply(pose)

	return pose

def Relax(pose=None, cstFile=None, w4distance=1., w4angle=1., w4dihedral=1., w4phipsi=10., w4beta=1.5, ncycles=5):

        assert pose is not None
        if not os.path.isfile(cstFile):
                print 'ERROR: invalid constraint file: ', cstFile
                exit(1)

        constraints = protocols.constraint_movers.ConstraintSetMover()
        constraints.add_constraints(True)
        constraints.constraint_file(cstFile)
        constraints.apply(pose)

        scorefxn = get_fa_scorefxn()
        scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, w4distance)
        scorefxn.set_weight(rosetta.core.scoring.dihedral_constraint, w4dihedral)
        scorefxn.set_weight(rosetta.core.scoring.angle_constraint, w4angle)

        scorefxn.set_weight(rosetta.core.scoring.hbond_lr_bb, w4beta)
        scorefxn.set_weight(rosetta.core.scoring.fa_elec, w4beta)
        scorefxn.set_weight(rosetta.core.scoring.ss_pair, w4beta)
        scorefxn.set_weight(rosetta.core.scoring.sheet, w4beta)

        scorefxn.show(pose)

        fastrelax = rosetta.protocols.relax.FastRelax(scorefxn, ncycles)
        fastrelax.apply(pose)

        scorefxn.show(pose)
        return pose


def main(argv):

	if len(argv) < 2:
		Usage()
		exit(1)

	inputFile = argv[0]
	cstFile = argv[1]

	try:
                opts, args = getopt.getopt(argv[2:],"w:d:a:p:n:t:s:v:rbq",[ "w4distance=", "w4dihedral=", "w4angle=", "w4phipsi", "ncycles=", "tolerance=", "w4beta=", "savefolder=", "randomPhiPsi=", "nonblist=", "quick="])
                print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

	w4phipsi = 5.
	w4dihedral = 1.
	w4angle = 1.
	w4distance = 1.0
	w4beta = 2
	ncycles = 2000
	tolerance = 0.0001

	UseNBList = True
	Quick = False
	SampleByPhiPsiDistribution = True

	savefolder = './'

        for opt, arg in opts:
                if opt in ("-w", "--w4distance"):
                        w4distance = np.float32(arg)
                elif opt in ("-d", "--w4dihedral"):
                        w4dihedral = np.float32(arg)
		elif opt in ("-a", "--w4angle"):
			w4angle = np.float32(arg)
		elif opt in ("-p", "--w4phipsi"):
			w4phipsi = np.float32(arg)
		elif opt in ("-s", "--w4beta"):
			w4beta = np.float32(arg)
	
		elif opt in ("-n", "--ncycles"):
			ncycles = np.int32(arg)
		elif opt in ("-t", "--tolerance"):
			tolerance = np.float32(arg)
			assert tolerance>0

		elif opt in ("-v", "--savefolder"):
			savefolder = arg

		elif opt in ("-r", "--randomPhiPsi"):
			SampleByPhiPsiDistribution = False

		elif opt in ("-b", "--nonblist"):
			UseNBList = False

		elif opt in ("-q", "--quick"):
			Quick = True

                else:
                        Usage()
                        exit(1)

	if w4dihedral <= 0:
		w4dihedral = w4phipsi

	target = os.path.basename(inputFile).split('.')[0]

	PhiPsiDistribution = None

	if SampleByPhiPsiDistribution:
		assert inputFile.endswith('.fasta') or inputFile.endswith('.seq')
		seq = ReadFASTAFile(inputFile)
		PhiPsiDistribution = ExtractPhiPsiDistribution(seq, cstFile)

	init()
	rosetta.basic.options.set_boolean_option( 'run:nblist_autoupdate', True )
	pose = InitializePose(inputFile, PhiPsiDistribution)

	switch = SwitchResidueTypeSetMover("centroid")
	switch.apply(pose)

	pose = Fold(pose, cstFile, w4phipsi=w4phipsi, w4angle=w4angle, w4dihedral=w4dihedral, w4distance=w4distance, w4beta=w4beta, tolerance=tolerance, ncycles=ncycles, UseNBList=UseNBList)

	modelFile = os.path.join(savefolder, target + '.fold.' + str(os.getpid()) + '.pdb')
	pose.dump_pdb(modelFile)

	if not Quick:
		pose = Relax(pose, cstFile, w4phipsi=w4phipsi, w4angle=w4angle, w4dihedral=w4dihedral, w4distance=w4distance, w4beta=w4beta)
		modelFile = os.path.join(savefolder, target + '.relaxed.' + str(os.getpid()) + '.pdb')
		pose.dump_pdb(modelFile)

if __name__ == "__main__":
        main(sys.argv[1:])

