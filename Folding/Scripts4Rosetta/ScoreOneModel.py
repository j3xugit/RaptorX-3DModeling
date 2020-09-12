import os
import sys
import numpy as np
import getopt

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta import *

from FoldNRelax import CheckCSTFile

def Usage():
	print 'python ScoreOneModel.py inputFile cstFile'
	print '       inputFile: a PDB file ending with .pdb'
	print '       cstFile: a constraint file including all kinds of constraints acceptable by Rosetta'

def PrintScore(score, savefile=None, ResidueLevel=True):
	overall = [score['totalPot'], score['totalDistPot'], score['totalDihedralPot'], score['totalAnglePot'] ]
	average = [ e/len(score['distPots']) for e in overall ]

	overall = [ '{:.2f}'.format(e) for e in overall ]
	average = [ '{:.2f}'.format(e) for e in average ]

	header = '\t'.join( ['total'] + overall )
	avgStr = '\t'.join( ['avg'] + average )

	allLines = [header, avgStr]

	if ResidueLevel:
		for i, distPot, diPot, aPot in zip(range(len(score['distPots'])), score['distPots'], score['dihedralPots'], score['distPots']):
			totPot = distPot + diPot + aPot
			line ='\t'.join( [ str(i+1), '{:.2f}'.format(totPot), '{:.2f}'.format(distPot), '{:.2f}'.format(diPot), '{:.2f}'.format(aPot) ] )
			allLines.append(line)

	if savefile is None:
		print '\n'.join(allLines)
		return

	with open(savefile, 'w') as fh:
		fh.writelines('\n'.join(allLines) )

def GetScore(pose):
	score = dict()
	totalDistPot= pose.energies().total_energies()[atom_pair_constraint]
	totalAnglePot= pose.energies().total_energies()[angle_constraint]
	totalDihedralPot= pose.energies().total_energies()[dihedral_constraint]
	#print totalDistPot, totalAnglePot, totalDihedralPot

	score['totalDistPot'] = totalDistPot
	score['totalAnglePot'] = totalAnglePot
	score['totalDihedralPot'] = totalDihedralPot
	score['totalPot'] = totalDistPot + totalAnglePot + totalDihedralPot

	#print score

	distPots = [ pose.energies().residue_total_energies(i)[atom_pair_constraint] for i in np.arange(1, pose.total_residue()+1 ) ]
	anglePots = [ pose.energies().residue_total_energies(i)[angle_constraint] for i in np.arange(1, pose.total_residue()+1 ) ]
	dihedralPots = [ pose.energies().residue_total_energies(i)[dihedral_constraint] for i in np.arange(1, pose.total_residue()+1 ) ]

	score['distPots'] = distPots
	score['anglePots'] = anglePots
	score['dihedralPots'] = dihedralPots

	return score

def Score(pose=None, scorefxn=None, constraints=None):
        assert pose is not None
        constraints.apply(pose)
        scorefxn(pose)

	return GetScore(pose)

def main(argv):

	if len(argv) < 2:
		Usage()
		exit(1)

	inputFile = argv[0]
	cstFile = argv[1]

	assert inputFile.endswith('.pdb') 
	target = os.path.basename(inputFile).split('.')[0]

	init('-out:level 100')
	rosetta.basic.options.set_boolean_option( 'run:nblist_autoupdate', True )

	if not os.path.isfile(cstFile):
                print 'ERROR: invalid constraint file: ', cstFile
                exit(1)

	if not CheckCSTFile(cstFile):
		print 'ERROR: incorrect cst file: ', cstFile
		exit(1)

        constraints = protocols.constraint_movers.ConstraintSetMover()
        constraints.add_constraints(True)
        constraints.constraint_file(cstFile)

        scorefxn = create_score_function('ref2015')
        scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, 1)
        scorefxn.set_weight(rosetta.core.scoring.dihedral_constraint, 1)
        scorefxn.set_weight(rosetta.core.scoring.angle_constraint, 1)

	pose = pose_from_pdb(inputFile)
	score = Score(pose, scorefxn, constraints)
	score['modelFile'] = inputFile

	PrintScore(score)
	#print os.path.basename(s['modelFile']), s['totalPot'], s['totalDistPot'], s['totalAnglePot'], s['totalDihedralPot']

if __name__ == "__main__":
        main(sys.argv[1:])

