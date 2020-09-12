import os
import sys
import numpy as np
import cPickle
import math
import getopt

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import quasi_newton, conjugate_gradients, actions

from Common.LoadFASTA import LoadFASTAFile

from DL4DistancePrediction3.Utils.SplineCurves4Potential import GenSplineFunc
from DL4DistancePrediction3.config import  Response2LabelName
from DL4PropertyPrediction.PropertyUtils import LoadPredictedProperties

eps = np.finfo(np.float32).eps

def Usage():
	print 'RunModeller.py pdbFile/seqFile distPotential_PKL [-p predictedProperty_PKL ] [-w weight4Dist:weight4PhiPsi] [-s minSeqSep for dist potential ] [ -a algorithm]'
	print '  -p: specify the predicted property file in PKL '
	print '  -w: if two values provided, it specifies weight for both distance and PhiPsi potentials. If only one value, then weight for distance. Default value: 1.'
	print '  -s: minimum sequence separation between two residues for which their distance potential is used, default 6'
	print '  -a: optimization method, QN (quasi-newton) or CG (conjugate gradient), default CG'

def LoadPredictedPhiPsi(file):
	if not os.path.isfile(file):
		print 'ERROR: cannot open file: ', file
		exit(1)

	content = LoadPredictedProperties(file)
        assert len(content) >=3
        name, sequence, predProperty = content[:3]

	if not predProperty.has_key('PhiPsi_vonMise2d4'):
                print 'ERROR: the property file does not contain predicted Phi and Psi information: ', file
                exit(1)

        PhiPsiList = predProperty['PhiPsi_vonMise2d4']

	Phis = []
	Psis = []
	for i, PhiPsi in enumerate(PhiPsiList):
		if i>0:
			force = 1./(eps + PhiPsi[2])
			phase = -PhiPsi[0]

			if phase < 0: ## do this conversion for MODELLER
				phase = np.pi + phase
				force = - force

			Phis.append((force, phase))
		else:
			Phis.append((None, None))

		if i < len(PhiPsiList) - 1:
			force = 1./(eps + PhiPsi[3])
			phase = -PhiPsi[1]

			if phase < 0: ## do this conversion for MODELLER
				phase = np.pi + phase
				force = - force

			Psis.append((force, phase))
		else:
			Psis.append((None, None))

	return (Phis, Psis)
		


def CreatePhiPsiConstraints(mdl, propertyFile, strength=1.0):

	Phis, Psis = LoadPredictedPhiPsi(propertyFile)
	rsrs = []
	for res in mdl.residues:
		#print res.index, res.phi, res.psi
		ind = res.index -1 
		if res.phi is not None and Phis[ind][0] is not None:
			force = Phis[ind][0] * strength
			phase = Phis[ind][1]
			period = 1
			r = forms.cosine(group=physical.phi_dihedral, feature=features.dihedral(res.phi.atoms), phase=phase, force=force, period=period)
			rsrs.append(r)
		if res.psi is not None and Psis[ind][0] is not None:
			force = Psis[ind][0] * strength
			phase = Psis[ind][1]
			period = 1
			r = forms.cosine(group=physical.psi_dihedral, feature=features.dihedral(res.psi.atoms), phase=phase, force=force, period=period)
			rsrs.append(r)

	return rsrs
		

def Dist2Energy(d, func, dfunc, deriv=True):
	f = 0.
	g = 0.
	if d < 16.:
		f = func(d)
		g = dfunc(d)

	if deriv:
		return f, g
	return f

class MyRGEnergyTerm(terms.energy_term):
    	_physical_type = physical.ca_distance
	_cutoff = 20.

    	def __init__(self, atomName=None, strength=1.0):
        	self.strength = strength
		self.atomName = atomName

	def eval(self, mdl, deriv, indats):
        	atoms = self.indices_to_atoms(mdl, indats)

		x0 = []
		## calculate center of coordinates and squared coordinates
		xsqr = ysqr = zsqr = 0.
		cx = cy = cz = 0.
		numAtoms = 0		
		for (num, a) in enumerate(atoms):
			if self.atomName is not None and a.name != self.atomName:
				continue
			x0.append(a.x)
			cx += a.x
			cy += a.y
			cz += a.z
			numAtoms += 1

			xsqr += a.x ** 2
			ysqr += a.y ** 2
			zsqr += a.z ** 2

		cx /= numAtoms
		cy /= numAtoms
		cz /= numAtoms

		xsqr /= numAtoms
		ysqr /= numAtoms
		zsqr /= numAtoms

		## calcualte radius of gyration
		R = math.sqrt( (xsqr - cx ** 2) + (ysqr - cy ** 2) + (zsqr - cz ** 2) )
		e = R * self.strength

		print "Radius of Gyration: ", R
		#print x0

		if not deriv:
			return e

        	dvx = [0.] * len(indats)
        	dvy = [0.] * len(indats)
        	dvz = [0.] * len(indats)

		## calculate derivative
		for (num, a) in enumerate(atoms):
			if self.atomName is not None and a.name != self.atomName:
				continue
			## derivate from center position
			dvx[num] -= cx / numAtoms
			dvy[num] -= cy / numAtoms
			dvz[num] -= cz / numAtoms

			## derivate from xsqr, ysqr and zsqr
			dvx[num] += a.x/ numAtoms
			dvx[num] += a.y/ numAtoms
			dvx[num] += a.z/ numAtoms

			## scale
			dvx[num] *= (self.strength / (eps + R))
			dvy[num] *= (self.strength / (eps + R))
			dvz[num] *= (self.strength / (eps + R))

		#print dvx

		return (e, dvx, dvy, dvz)

class MyPairwiseEnergyTerm(terms.energy_term):

    	_physical_type = physical.xy_distance
	_cutoff = 16.5

   	#potMatrix is a tensor with dimension L*L*numBins where L is the seqLena and numBins is the number of distance bins
	#midPoints is the middle point of each distance bin
    	def __init__(self, potMatrix=None, distCutoff=None, atomNames=None, strength=1.0, minSeqSep=6):
        	self.strength = strength
		self.atomNames = atomNames
		assert atomNames is not None

		if atomNames[0] == 'N' and atomNames[1] == 'O':
			self._physical_type = physical.n_o_distance
		elif atomNames[0] == 'CA' and atomNames[1] == 'CA':
			self._physical_type = physical.ca_distance
		
        	terms.energy_term.__init__(self)


		self.distCutoff = distCutoff
		self.minSeqSep = minSeqSep

		## add code here to generate splines (and derivative) from discrete distance potential
		#self.potMatrix = potMatrix
		self.funcs, self.dfuncs = GenSplineFunc(potMatrix, distCutoff)

    	def eval(self, mdl, deriv, indats):
        	atoms = self.indices_to_atoms(mdl, indats)
        	e = 0.
        	dvx = [0.] * len(indats)
        	dvy = [0.] * len(indats)
        	dvz = [0.] * len(indats)

		for (num1, a1) in enumerate(atoms):
			if a1.name != self.atomNames[0]:
				continue
			rsd1 = a1.residue.index - 1

			for (num2, a2) in enumerate(atoms):
				if a2.name != self.atomNames[1]:
					continue
				rsd2 = a2.residue.index -1

				if rsd2 < (rsd1 + self.minSeqSep ):
					continue


				d = math.sqrt( eps + (a1.x-a2.x)**2 + (a1.y-a2.y)**2 + (a1.z-a2.z)**2 )
				f, g = Dist2Energy(d, self.funcs[rsd1][rsd2], self.dfuncs[rsd1][rsd2], deriv=True) 

				#print 'in MyPairwiseEnergyTerm.eval():', rsd1, rsd2, 'dist, energy, deriv=', d, f, g

				e += f

				if deriv: 
					dvx[num1] += (a1.x - a2.x)/(eps + d) * g 
					dvx[num2] += (a2.x - a1.x)/(eps + d) * g 

					dvy[num1] += (a1.y - a2.y)/(eps + d) * g
					dvy[num2] += (a2.y - a1.y)/(eps + d) * g 

					dvz[num1] += (a1.z - a2.z)/(eps + d) * g 
					dvz[num2] += (a2.z - a1.z)/(eps + d) * g

		e *= self.strength

		if deriv:
			dvx = [ x * self.strength for x in dvx ]
			dvy = [ y * self.strength for y in dvy ]
			dvz = [ z * self.strength for z in dvz ]
			"""
			print dvx
			print dvy
			print dvz
			"""

        	if deriv:
            		return (e, dvx, dvy, dvz)
        	else:
            		return e

def main(argv):

	if len(argv) < 2:
		Usage()
		exit(1)

	inputFile = argv[0]
	distPotentialFile = argv[1]
	minSeqSep = 6
	propertyFile = None
	w4Dist = 1.
	w4PhiPsi = 1.
	algorithm = 'CG'
	allAlgorithms = set(['CG', 'QN'])

	if not os.path.isfile(inputFile):
		print 'ERROR: the input file does not exist: ', inputFile
		exit(1)

	if not os.path.isfile(distPotentialFile):
		print 'ERROR: the input distance potential file does not exist: ', distPotentialFile
		exit(1)

	try:
                opts, args = getopt.getopt(argv[2:],"p:s:w:a:",["propertyFile=", "minSeqSep=", "weight=", "algorithm="])
        except getopt.GetoptError:
                Usage()
                exit(1)

        for opt, arg in opts:
		if opt in ("-p", "propertyFile"):
			propertyFile = arg
			if not os.path.isfile(propertyFile):
				print 'ERROR: the file does not exist: ', propertyFile
				exit(1)

		elif opt in ("-s", "minSeqSep"):
			minSeqSep = np.int32(arg)
			if minSeqSep < 1:
				print 'ERROR: a positive integer is required for minSeqSep'
				exit(1)

		elif opt in ("-w", "weight"):
			fields = arg.split(':')
			if len(fields) < 1 or len(fields) > 2:
				print 'ERROR: wrong input format for weight factors'
				exit(1)		

			w4Dist = np.float32(fields[0])			
			if len(fields) == 2:
				w4PhiPsi = np.float32(fields[1])

		elif opt in ("-a", "algorithm"):
			algorithm = arg
			if algorithm not in allAlgorithms:
				print 'ERROR: only the following algorithms are allowed: ', allAlgorithms
				exit(1)	
		else:
			Usage()
			exit(1)

	assert w4Dist >= 0
	assert w4PhiPsi >=0 

	inputIsPDB = True

	if inputFile.endswith('.pdb'):
		inputIsPDB = True
		pdbFile = inputFile
	elif inputFile.endswith('.seq') or inputFile.endswith('.fasta'):
		inputIsPDB = False
		seqFile = inputFile
	else:
		print 'ERROR: the input file shall be either a seqeunce file (ends with .seq or .fasta) or a PDB file (ends with .pdb)'
		exit(1)

	
	distf = open(distPotentialFile, 'r')
	name, seq4potential, distPotMatrices, distCutoffs = cPickle.load(distf)
	distf.close()

	seed = max(2, os.getpid() % 50000)
	env = environ(rand_seed = -seed)
	# Read parameters (needed to build models from internal coordinates)
	env.libs.topology.read('${LIB}/top_heav.lib')
	env.libs.parameters.read('${LIB}/par.lib')

	#env.edat.contact_shell = 20.
	#env.edat.energy_terms.append( MyRGEnergyTerm(atomName='CA', strength=100.0) )
	#env.edat.energy_terms.append( MyRGEnergyTerm(strength=10.0) )

	for response, pot in distPotMatrices.items():
		apt = Response2LabelName(response)
		print 'response: ', response
		if apt == 'CaCa':
			atomNames = ('CA', 'CA')
		elif apt == 'NO':
			atomNames =('N', 'O')
		elif apt == 'CbCb':
			atomNames = ('CB', 'CB')
		else: 
			continue
		"""
		elif apt == 'CgCg':
			atomNames = ('CG', 'CG')
		elif apt == 'CaCg':
			atomNames = ('CA', 'CG')
		"""
		env.edat.energy_terms.append( MyPairwiseEnergyTerm(potMatrix=pot, atomNames=atomNames, distCutoff=distCutoffs[response], strength=4.0, minSeqSep=6) )

	#print'energy terms:',  env.edat.energy_terms, len(env.edat.energy_terms)

	if not inputIsPDB:
		sequence = LoadFASTAFile(seqFile).upper()
		assert sequence == seq4potential.upper()
		# Build a model from one-letter codes, and write to a PDB file:
		mdl = model(env)
		mdl.build_sequence(sequence)
	else:
		mdl = model(env)
		mdl = complete_pdb(env, pdbFile)

		## assure that the sequence in the PDB file is same as that in the potential file
		assert len(seq4potential) == len(mdl.residues)
		for res, mres in zip(seq4potential, mdl.residues):
			assert res == mres.code

	## add stereo restraints and dihedral restraints
	sel = selection(mdl).only_atom_types('CB') | selection(mdl).only_mainchain()
	if not inputIsPDB:
		## generate a random 3D model as initial
		dih = 'phi psi omega chi1 chi2 chi3 chi4 chi5'
		sel.rotate_dihedrals(change='RANDOMIZE', deviation=90.0, dihedrals=dih)
		#sel.randomize_xyz(deviation=(-1.5*math.sqrt(len(sequence)) ) )
		mdl.write(file=name + '-init.pdb')
	else:
		sel.randomize_xyz(deviation=0.5)
		mdl.write(file=name + '-init.pdb')

	mdl.restraints.make(sel, restraint_type='stereo', spline_on_site=False)
	if propertyFile is not None:
		rsrs = CreatePhiPsiConstraints(mdl, propertyFile, strength=2.0)
		for rsr in rsrs:
			mdl.restraints.add(rsr)

	mdl.restraints.write(name + '-rsr.txt')
	#exit(0)

	print "initial energy: ", sel.energy()

	## minimize energy
	if algorithm=='CG': 
		optimizer=conjugate_gradients()
	else:
		optimizer=quasi_newton()

	optimizer.optimize(sel, max_iterations=800, actions=actions.trace(1, name + '-opt-logfile.txt'))

	print "final energy:", sel.energy()
	mdl.write(file=name + '-final.pdb')
	##sel.write(file=name + '-final.pdb')

if __name__ == "__main__":
        main(sys.argv[1:])

