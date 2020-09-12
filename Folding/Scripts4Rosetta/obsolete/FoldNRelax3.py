import os
import sys
import numpy as np
import getopt

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta import *

from pyrosetta.rosetta.protocols.minimization_packing import MinMover

def Usage():
	print 'python FoldNRelax.py inputFile cstFile [ -w weight4distance ] [ -d weight4dihedral ] [-a weight4angle] [-e w4beta] [ -n ncycles4initFolding ] [-t tolerance] [-s savefolder ] [-r] [-b] [-q] [-p]'
	print '       inputFile: a file ending with .fasta, .seq, and .pdb'
	print '       cstFile: a constraint file including all kinds of constraints acceptable by Rosetta'
	print '       -n: the number of iterations for each MinMover (default 1000) at initial folding stage'
	print '	      -t: the tolerance for termination of initial folding, default 0.0001, i.e., 0.01%'
	print '	      -s: the folder for result save, default current work directory'
	print '	'
	print '       -w: weight for distance constraints used by FastRelax, default 1.0'
	print '       -d: weight for orientation dihedral used by FastRelax, default 1.0'
	print '	      -a: weight for orientation angle used by FastRelax, default 1.0'
	print '	      -e: weight for beta sheet used by FastRelax, default None'
	print ' '
	print '       -r: if specified, do not sample phi/psi angles by the predicted distribution, otherwise sample (default)'
	print '       -b: if specified, do not use neighbor list in minimization, otherwise use neighbor list (default)'
	print '	      -q: if specified, do not run FastRelax, otherwise run FastRelax (default)'
	print '	      -p: if specified, use perturbation in the first stage of folding, default No'

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
			print 'WARNING: unknown dihedral type in line: ', r
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

def GetAngles(pose):
	phis = np.array([ pose.phi(i) for i in range(1, pose.total_residue() + 1 ) ])
	psis = np.array([ pose.psi(i) for i in range(1, pose.total_residue() + 1 ) ])

	return phis, psis

def SetAngles(pose, phis, psis):	
	for i, phi, psi in zip(range(1, pose.total_residue() + 1), phis, psis):
		pose.set_phi(i, phi)
		pose.set_psi(i, psi)

def InitializePose(inputFile, PhiPsiDistribution=None):
	if inputFile.endswith('.fasta') or inputFile.endswith('.seq') or inputFile.endswith('.fa'):
		sequence = ReadFASTAFile(inputFile)
		pose = pose_from_sequence(sequence)

		if PhiPsiDistribution is None:
			phis = np.random.uniform(-180, 180, pose.total_residue() )
			psis = np.random.uniform(-180, 180, pose.total_residue() )
			SetAngles(pose, phis, psis)
			"""
			for i in range(1, pose.total_residue() + 1):
				pose.set_phi(i, np.random.uniform(-180, 180))
				pose.set_psi(i, np.random.uniform(-180, 180))
			"""
		else:
			PhiDistribution = PhiPsiDistribution['phi']
			phis = SampleDihedralsByAMBER(PhiDistribution)
			PsiDistribution = PhiPsiDistribution['psi']
			psis = SampleDihedralsByAMBER(PsiDistribution)

			SetAngles(pose, phis, psis)
			"""
			for i, phi, psi in zip(range(1, pose.total_residue() + 1), phis, psis):
				pose.set_phi(i, phi)
				pose.set_psi(i, psi)
			"""
	else:
		pose = pose_from_pdb(inputFile)

	return pose

## angle is an array and each element is degree of an angle instead of radian
def CorrectAngles(angles):
	temp = angles % 360
	np.putmask(temp, temp>180, temp-360)
	return temp

## angles is an array and each element is degree but not radian
## perturb angles by adding a noise sampled from normal distribution N(0, sigma)
def PerturbAngles(angles, sigma):
	noise = np.random.normal(0, sigma, angles.shape).astype(np.float32)
	return CorrectAngles( angles + noise )

## minimize the potential by repeating perturbation and optimization
def MinimizeEnergyByPerturbation(pose, mover, scorefunc, sigmas=[10, 7.5, 3, 2], show=True):

	initf = scorefunc(pose)
	initphis, initpsis = GetAngles(pose)

	bestphis, bestpsis = initphis, initpsis
	bestf = initf

	numRounds = 0
	i=0
	converged = 0
	while i < len(sigmas):
		s = sigmas[i]
		newphis = PerturbAngles(bestphis, s)
		newpsis = PerturbAngles(bestpsis, s)
		SetAngles(pose, newphis, newpsis)
		mover.apply(pose)
		numRounds += 1

		f = scorefunc(pose)
		if f >= bestf:
			converged += 1
		else:
			converged = 0
			bestf = f
			bestphis, bestpsis = GetAngles(pose)

		if converged > 1:
			## skip some iterations
			i = max(i, len(sigmas)-2 )
		i += 1

	## set Phi/Psi to the best
	SetAngles(pose, bestphis, bestpsis)

	if show:
		print 'init score: ', initf
		print 'final score: ', bestf
		print 'improvement: ', np.round( (initf - bestf)/abs(initf) * 100, 3)
		print '#rounds of perturbation: ', numRounds

	return pose

def Fold(pose=None, cstFile=None, ncycles=1000, tolerance=0.0001, UseNBList=True, UsePerturbation=False):
	assert pose is not None
	if not os.path.isfile(cstFile):
		print 'ERROR: invalid constraint file: ', cstFile
		exit(1)

	scriptdir = os.path.dirname(os.path.realpath(__file__))

	constraints = protocols.constraint_movers.ConstraintSetMover()
	constraints.constraint_file(cstFile)
	constraints.add_constraints(True)
	constraints.apply(pose)

	mmap = MoveMap()
	mmap.set_bb(True)
	mmap.set_chi(False)
	mmap.set_jump(True)
	mmap.show()

	sf = ScoreFunction()
    	sf.add_weights_from_file(scriptdir + '/params/scorefxn.wts')

    	sf1 = ScoreFunction()
    	sf1.add_weights_from_file(scriptdir + '/params/scorefxn1.wts')

    	sf_vdw = ScoreFunction()
    	sf_vdw.add_weights_from_file(scriptdir + '/params/scorefxn_vdw.wts')

    	sf_cart = ScoreFunction()
    	sf_cart.add_weights_from_file(scriptdir + '/params/scorefxn_cart.wts')

	min_mover = MinMover(mmap, sf, 'lbfgs_armijo_nonmonotone', tolerance, True)
    	min_mover.max_iter(ncycles)

    	min_mover1 = MinMover(mmap, sf1, 'lbfgs_armijo_nonmonotone', tolerance, True)
    	min_mover1.max_iter(ncycles)

    	min_mover_vdw = MinMover(mmap, sf_vdw, 'lbfgs_armijo_nonmonotone', tolerance, True)
    	min_mover_vdw.max_iter(500)

    	min_mover_cart = MinMover(mmap, sf_cart, 'lbfgs_armijo_nonmonotone', tolerance, True)
    	min_mover_cart.max_iter(ncycles)
    	min_mover_cart.cartesian(True)

    	repeat_mover = RepeatMover(min_mover, 3)
	repeat_mover.apply(pose)

	if UsePerturbation:
		pose = MinimizeEnergyByPerturbation(pose, min_mover, sf, sigmas=[10, 7.5, 3, 2])

	min_mover_cart.apply(pose)

	switch = SwitchResidueTypeSetMover("fa_standard")
	switch.apply(pose)

	return pose

def Relax(pose=None, cstFile=None, w4distance=1.0, w4angle=1., w4dihedral=1., w4beta=None, ncycles=5):
        assert pose is not None
        if not os.path.isfile(cstFile):
                print 'ERROR: invalid constraint file: ', cstFile
                exit(1)

        constraints = protocols.constraint_movers.ConstraintSetMover()
        constraints.add_constraints(True)
        constraints.constraint_file(cstFile)
        constraints.apply(pose)

        scorefxn = create_score_function('ref2015')
        scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, w4distance)
        scorefxn.set_weight(rosetta.core.scoring.dihedral_constraint, w4dihedral)
        scorefxn.set_weight(rosetta.core.scoring.angle_constraint, w4angle)

	if w4beta is not None:
		scorefxn.set_weight(rosetta.core.scoring.hbond_lr_bb, w4beta)
        	scorefxn.set_weight(rosetta.core.scoring.fa_elec, w4beta)
        	scorefxn.set_weight(rosetta.core.scoring.ss_pair, w4beta)
        	scorefxn.set_weight(rosetta.core.scoring.sheet, w4beta)

        scorefxn.show(pose)

	mmap = MoveMap()
        mmap.set_bb(True)
        mmap.set_chi(True)
        mmap.set_jump(True)

        fastrelax = rosetta.protocols.relax.FastRelax(scorefxn, ncycles)
	fastrelax.max_iter(200)
	fastrelax.dualspace(True)
	fastrelax.set_movemap(mmap)
        fastrelax.apply(pose)

        scorefxn.show(pose)

        return pose


def main(argv):

	if len(argv) < 2:
		Usage()
		exit(1)

	inputFile = argv[0]
	cstFile = argv[1]

	assert inputFile.endswith('.pdb') or inputFile.endswith('.fasta') or inputFile.endswith('.seq')

	try:
                opts, args = getopt.getopt(argv[2:],"w:d:a:e:n:t:s:rbqp",[ "w4distance=", "w4dihedral=", "w4angle=", "w4beta=", "ncycles=", "tolerance=", "savefolder=", "randomPhiPsi=", "nonblist=", "quick=", "perturb="])
                print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

	w4dihedral = 1.
	w4angle = 1.
	w4distance = 1.0
	w4beta = None

	ncycles = 1000
	tolerance = 0.0001

	UseNBList = True
	DoRelax = True
	SampleByPhiPsiDistribution = True

	UsePerturbation = False

	savefolder = os.getcwd()

        for opt, arg in opts:
                if opt in ("-w", "--w4distance"):
                        w4distance = np.float32(arg)
                elif opt in ("-d", "--w4dihedral"):
                        w4dihedral = np.float32(arg)
		elif opt in ("-a", "--w4angle"):
			w4angle = np.float32(arg)
		elif opt in ("-e", "--w4beta"):
			w4beta = np.float32(arg)

		elif opt in ("-n", "--ncycles"):
			ncycles = np.int32(arg)
		elif opt in ("-t", "--tolerance"):
			tolerance = np.float32(arg)
			assert tolerance>0

		elif opt in ("-s", "--savefolder"):
			savefolder = arg

		elif opt in ("-r", "--randomPhiPsi"):
			SampleByPhiPsiDistribution = False

		elif opt in ("-b", "--nonblist"):
			UseNBList = False

		elif opt in ("-q", "--quick"):
			DoRelax = False

		elif opt in ("-p", "--perturb"):
			UsePerturbation = True

                else:
                        Usage()
                        exit(1)

	target = os.path.basename(inputFile).split('.')[0]

	PhiPsiDistribution = None

	if (not inputFile.endswith('.pdb')) and SampleByPhiPsiDistribution:
		seq = ReadFASTAFile(inputFile)
		PhiPsiDistribution = ExtractPhiPsiDistribution(seq, cstFile)

	#init()
	init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 100')
	rosetta.basic.options.set_boolean_option( 'run:nblist_autoupdate', True )
	pose = InitializePose(inputFile, PhiPsiDistribution)
	if pose is None:
		print 'ERROR: the intial pose is None'
		exit(1)

	switch = SwitchResidueTypeSetMover("centroid")
	switch.apply(pose)

	pose = Fold(pose, cstFile, tolerance=tolerance, ncycles=ncycles, UseNBList=UseNBList, UsePerturbation=UsePerturbation)
	if pose is None:
		print 'ERROR: the folded pose is None'
		exit(1)

	modelFile = os.path.join(savefolder, target + '.fold.' + str(os.getpid()) + '.pdb')
	print 'Writing one model to ', modelFile
	pose.dump_pdb(modelFile)

	if DoRelax:
		pose.remove_constraints()
		pose = Relax(pose, cstFile, w4angle=w4angle, w4dihedral=w4dihedral, w4distance=w4distance, w4beta=w4beta)
		modelFile = os.path.join(savefolder, target + '.relaxed.' + str(os.getpid()) + '.pdb')
		print 'Writing one model to ', modelFile
		pose.dump_pdb(modelFile)

if __name__ == "__main__":
        main(sys.argv[1:])

