import os
import sys
import numpy as np
import getopt
import random
import shutil
import tempfile
import socket

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta import *

from pyrosetta.rosetta.protocols.minimization_packing import MinMover

#from ScoreOneModel import GetScore
#from ScoreOneModel import PrintScore

def Usage():
	print 'python FoldNRelax.py inputFile cstFile [-I extraInput ] [ -w weight4distance ] [ -d weight4dihedral ] [-a weight4angle] [-e w4beta] [ -n ncycles4initFolding ] [-t tolerance] [-s savefolder ] [-r] [-b] [-q] [-p]'
	print '       inputFile: a file ending with .fasta, .seq, and .pdb'
	print '       cstFile: a constraint file including all kinds of constraints acceptable by Rosetta'
	print '	      -I: when this option specified, cstFile shall be interpreted as the predicted dist/orientation matrix file instead of Rosetta constraint file'
	print '	          extraInput is used to inform how to derive Rosetta constrains from predicted dist/ori matrix, e.g., potential type and some hyperparamters'
	print '	          It also provides predicted Phi/Psi angle file, e.g., -I alpha=1.61,phipsi=T1111.predictedProperties.pkl or -I alpha=1.61,T111.predictedProperties.pkl'
	print '       -n: the number of iterations for each MinMover (default 1000) at initial folding stage'
	print '	      -t: the tolerance for termination of initial folding, default 0.0001, i.e., 0.01%'
	print '	      -s: the folder for result save, default current work directory'
	print '	'
	print '       -w: weight for distance constraints used by FastRelax, default 0.2'
	print '       -d: weight for orientation dihedral used by FastRelax, default 0.2'
	print '	      -a: weight for orientation angle used by FastRelax, default 0.2'
	print '	      -e: weight for beta sheet used by FastRelax, default None'
	print ' '
	print '       -r: if specified, do not sample phi/psi angles by the predicted distribution, otherwise sample (default)'
	print '       -b: if specified, do not use neighbor list in minimization, otherwise use neighbor list (default)'
	print '	      -q: if specified, do not run FastRelax, otherwise run FastRelax (default)'
	print '       -p: if specified, use perturbation in the first stage of folding, default No'

def ReadFASTAFile(inputFile):
	with open(inputFile, 'r') as fh:  
    		sequence = fh.readlines()
  
    	# removing the trailing "\n" and any header lines
    	sequence = [line.strip() for line in sequence if not '>' in line]
    	sequence = ''.join( sequence )    # combine into a single sequence
	return sequence

## randomly check the SPLINE files in a cstfile to make sure the path information is correct
## PyRosetta itself does not check if any SPLINE files exist or not
def CheckCSTFile(cstfile, numChecks=10):
	with open(cstfile) as fh:
		content=[line.strip() for line in list(fh) if 'SPLINE' in line ]
	if len(content)<5:
		print 'ERROR: there are very few SPLINE constraints in', cstfile
		return False

	numChecks2 = min(numChecks, len(content) )
	for i in range(numChecks2):
		oneRow = random.choice(content)
		fields = oneRow.split()
		if len(fields) < 8:
			print 'ERROR: an incorrect SPLINE constraint in', cstfile, ':', oneRow
			return False
		location = fields.index('SPLINE')
		if location + 2 >= len(fields):
			print 'ERROR: an incorrect SPLINE constraint in', cstfile, ':', oneRow
			return False
		splinefile = fields[location+2]
		if not os.path.isfile(splinefile):
			print 'ERROR: invalid file for SPLINE constraint in', cstfile, ':', splinefile
			return False
	return True

## derive Rosetta constraints from pairMatixFile and propertyFile and save the constraints to saveFolder
## return the Rosetta constraint file to user
## the user shall be responsible for cleanup saveFolder after this cstfile is not needed any more
def DeriveRosettaCSTFile(seq, pairMatrixFile, propertyFile, saveFolder, param4Potential=1.61):

	import cPickle
	from Folding.GenPairwisePotentialFromPrediction import CalcDistOriPotential
	from Folding.GenPropertyPotential4Rosetta import GeneratePhiPsiPotential
	from GeneratePairPotential4Rosetta import GenerateSplinePotential, WriteSplineConstraints
	from DL4PropertyPrediction import PropertyUtils
	
	with open(pairMatrixFile, 'rb') as fh:
		name, sequence, distOriProbMatrix, contactProbMatrix, labelWeight, labelDistribution = cPickle.load(fh)[:6]

	if seq is not None and seq != sequence:
		print 'ERROR: query sequence is inconsistent with that in predicted distance/orientation file', pairMatrixFile
		exit(1)

	predData = (distOriProbMatrix, labelWeight, labelDistribution)
	pairPotential, cutoffs, validProbs, _, _ = CalcDistOriPotential(predData, param4Potential=param4Potential)
	pairConstraints = GenerateSplinePotential( (name, sequence, pairPotential, cutoffs, validProbs))

	print 'Finish calculating pairwise constraints from ', pairMatrixFile, 'and ', propertyFile

	## get a savefile, here bname could be a target name or an alignment name
	bname = os.path.basename(pairMatrixFile).split('.')[0]
	savefile = os.path.join(saveFolder, bname + '.pairPotential4Rosetta.SPLINE.txt')
	savefolder4histfile = os.path.join(saveFolder, 'SplinePotential4' + name + '/')
	WriteSplineConstraints(pairConstraints, savefile=savefile, savefolder4histfile=savefolder4histfile)

	print 'Finish writing pairwise constraints to ', savefile

	## handle PhiPsi angles
        name1, sequence1, predProperty = PropertyUtils.LoadPredictedProperties(propertyFile)[:3]

	if seq is not None and seq != sequence1:
		print 'ERROR: the query sequence is inconsistent with that in predicted property file', propertyFile
		exit(1)

	if sequence != sequence1:
		print 'ERROR: inconsistent sequences in the predicted distance/orientation file and the property file'
		print '	dist/ori file: ', pairMatrixFile
		print '	property file: ', propertyFile
		exit(1)

        if not predProperty.has_key('PhiPsi_vonMise2d4'):
                print 'ERROR: no predicted Phi/Psi in', propertyFile
                exit(1)

        PhiPsiList = predProperty['PhiPsi_vonMise2d4']
        PhiPsiConstraints = GeneratePhiPsiPotential(sequence, PhiPsiList)

        if len(PhiPsiConstraints)<1:
                print 'ERROR: failed to generate Phi/Psi constraints from', propertyFile
                exit(1)

	with open(savefile, 'a') as fh:
		fh.write('\n'.join(PhiPsiConstraints) )

	print 'The resultant Rosetta constraint is saved to ', savefile

	return savefile
	

def GetAngles(pose):
        phis = np.array([ pose.phi(i) for i in range(1, pose.total_residue() + 1 ) ])
        psis = np.array([ pose.psi(i) for i in range(1, pose.total_residue() + 1 ) ])
        return phis, psis

def SetAngles(pose, phis, psis):
        for i, phi, psi in zip(range(1, pose.total_residue() + 1), phis, psis):
                pose.set_phi(i, phi)
                pose.set_psi(i, psi)
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
	
## extract phi/psi distribution from cstFile, which has the Rosetta format with AMBER function
def ExtractPhiPsiDistribution(seq, cstFile):

	assert seq is not None
	seqLen = len(seq)	
	PhiDistribution = [None] * seqLen 
	PsiDistribution = [None] * seqLen 
	rows = None

	with open(cstFile, 'r') as f:
		content = f.readlines()
		rows = [ c for c in content if 'AMBER' in c and 'Dihedral' in c ]

	if len(rows) < 5:
		print "ERROR: there are very few Phi/Psi constraints in ", cstFile
		exit(1)

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

def InitializePose(inputFile, PhiPsiDistribution=None):
	if inputFile.endswith('.fasta') or inputFile.endswith('.seq'):
		sequence = ReadFASTAFile(inputFile)
		pose = pose_from_sequence(sequence)

		if PhiPsiDistribution is None:
			phis = np.random.uniform(-180, 180, pose.total_residue() )
                        psis = np.random.uniform(-180, 180, pose.total_residue() )
                        SetAngles(pose, phis, psis)
		else:
			PhiDistribution = PhiPsiDistribution['phi']
			phis = SampleDihedralsByAMBER(PhiDistribution)
			PsiDistribution = PhiPsiDistribution['psi']
			psis = SampleDihedralsByAMBER(PsiDistribution)

			SetAngles(pose, phis, psis)
	else:
		pose = pose_from_pdb(inputFile)

	return pose

def RemoveClash(scorefxn, mover, pose):
    for _ in range(0, 5):
        if float(scorefxn(pose)) < 10:
            break
        mover.apply(pose)

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

## pose shall already contain some constraints
def Fold(pose, ncycles=1000, tolerance=0.0001, UseNBList=True, UsePerturbation=False):
	assert pose is not None

	mmap = MoveMap()
	mmap.set_bb(True)
	mmap.set_chi(False)
	mmap.set_jump(True)
	mmap.show()

	scriptdir = os.path.dirname(os.path.realpath(__file__))
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

	## remove clash in the initial pose
	RemoveClash(sf_vdw, min_mover_vdw, pose)

    	repeat_mover = RepeatMover(min_mover, 4)
	repeat_mover.apply(pose)

	if UsePerturbation:
        	pose = MinimizeEnergyByPerturbation(pose, min_mover, sf, sigmas=[10, 7.5, 3, 2])

	min_mover_cart.apply(pose)

	RemoveClash(sf_vdw, min_mover1, pose)

	sf.show(pose)

	switch = SwitchResidueTypeSetMover("fa_standard")
	switch.apply(pose)

	return pose

def RelaxWithoutConstraints(pose, ncycles=1):
        assert pose is not None
	assert ncycles <=3 

        scorefxn = create_score_function('ref2015')
        scorefxn.show(pose)

        mmap = MoveMap()
        mmap.set_bb(True)
        mmap.set_chi(True)
        mmap.set_jump(True)

        fastrelax = rosetta.protocols.relax.FastRelax(scorefxn, ncycles)
        fastrelax.max_iter(200)
        #fastrelax.dualspace(True)
        fastrelax.set_movemap(mmap)
        fastrelax.apply(pose)
        scorefxn.show(pose)

        return pose

## pose shall already have some constraints
def Relax(pose, w4distance=0.2, w4angle=0.2, w4dihedral=0.2, w4beta=None, ncycles=5):
        assert pose is not None
	assert ncycles >= 2

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
	pairFile = argv[1]

	pairFileIsCST = True
	PhiPsiFile = None
	param4Potential = 1.61

	assert inputFile.endswith('.pdb') or inputFile.endswith('.fasta') or inputFile.endswith('.seq')
	if not os.path.isfile(inputFile):
		print 'ERROR: invalid input file for folding/relax:', inputFile
		exit(1)
	if not os.path.isfile(pairFile):
		print 'ERROR: invalid file for pairwise information:', pairFile
		exit(1)


	try:
                opts, args = getopt.getopt(argv[2:],"I:w:d:a:e:n:t:s:rbqp",["extraInput=", "w4distance=", "w4dihedral=", "w4angle=", "w4beta=", "ncycles=", "tolerance=", "savefolder=", "randomPhiPsi=", "nonblist=", "quick=", "perturbation="])
                #print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

	w4dihedral = 0.2
	w4angle = 0.2
	w4distance = 0.2
	w4beta = None

	ncycles = 1000
	tolerance = 0.0001

	UseNBList = True
	DoRelax = True
	SampleByPhiPsiDistribution = True

	UsePerturbation = False

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

                elif opt in ("-w", "--w4distance"):
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

	if UsePerturbation:
		print 'Will use Phi/Psi perturbation to minimize potential'

	target = os.path.basename(inputFile).split('.')[0]

	init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 100')
	rosetta.basic.options.set_boolean_option( 'run:nblist_autoupdate', True )

	seq = ReadFASTAFile(inputFile)
	seqLen = len(seq)

	if not pairFileIsCST:
		assert PhiPsiFile is not None
		if not os.path.isfile(PhiPsiFile):
			print 'ERROR: invalid file for predicted Phi/Psi angles: ', PhiPsiFile
			exit(1)

		if param4Potential > 10:
			param4Potential = random.uniform(1.57, 1.63)

		print 'alpha for DFIRE:', param4Potential
		assert param4Potential <= 1.63 and (param4Potential >= 1.57), 'In FoldNRelax.py, param for DFIRE distance potential is out of range'


		## create a CST folder and generate the CST file
		machine = socket.gethostname().split('.')[0]
		print 'Running FoldNRelax.py on', machine

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

	PhiPsiDistribution = None
	if (not inputFile.endswith('.pdb')) and SampleByPhiPsiDistribution:
		#seq = ReadFASTAFile(inputFile)
		PhiPsiDistribution = ExtractPhiPsiDistribution(seq, cstFile)

	pose = InitializePose(inputFile, PhiPsiDistribution)
	if pose is None:
		print 'ERROR: the intial pose is None'
		exit(1)

	switch = SwitchResidueTypeSetMover("centroid")
	switch.apply(pose)

	## read in constraints
        constraints = protocols.constraint_movers.ConstraintSetMover()
        constraints.constraint_file(cstFile)
        constraints.add_constraints(True)
	constraints.apply(pose)

	## remove the cstFile and its savefolder after adding the constraints
	if not pairFileIsCST:
		shutil.rmtree(cstfolder)

	pose = Fold(pose, tolerance=tolerance, ncycles=ncycles, UseNBList=UseNBList, UsePerturbation=UsePerturbation)
	if pose is None:
		print 'ERROR: the folded pose is None'
		exit(1)

	## some compute nodes have a short process ID. To avoid conflit, we add a random number after the process ID
	savefile = target + '.fold.' + str(os.getpid() * 1000 + np.random.randint(low=1, high=1000)) + '.pdb'

	if DoRelax:
		pose = RelaxWithoutConstraints(pose, ncycles=2)
		pose = Relax(pose, w4angle=w4angle, w4dihedral=w4dihedral, w4distance=w4distance, w4beta=w4beta)
		
		modelFile = os.path.join(savefolder, savefile[: -len('.pdb')] + '.relaxed.' + str(os.getpid()) + '.pdb')
		print 'Writing the relaxed model to ', modelFile
		pose.dump_pdb(modelFile)
	else:
		modelFile = os.path.join(savefolder, savefile)
		print 'Writing the initial model to ', modelFile
		pose.dump_pdb(modelFile)


if __name__ == "__main__":
        main(sys.argv[1:])

