import os
import sys
import numpy as np

sys.path.append(os.path.join(os.environ['DL4DistancePredHome'], 'Utils'))
from CalcAtomDistMatrixFromPDB import CalcAtomDistMatrix

sys.path.append(os.path.join(os.environ['ModelingHome'], 'Common'))
from LoadFASTA import LoadFASTAFile

def Usage():
	print "python CalcNOEViolationPerResidue.py seqFile modelFile restraintFile"
	print "	      seqFile: the file for primary sequence in FASTA format"
	print "	      modelFile: predicted 3D model file in PDB format"
	print "	      restraintFile: the distance restraint file used to build the model, e.g., contact.tbl"

def LoadDistRestraintFile(restraintFile, proteinLen):

	fh = open(restraintFile)
	content = [ line.strip() for line in list(fh) ]
	fh.close()

	restraints = dict()
	for line in content:
		fields = line.split()
		resid1 = int(fields[2])
		resid2 = int(fields[7])
		atom1 = fields[5][:-1][:2]
		atom2 = fields[10][:-1][:2]
		dist = np.float32(fields[11])
		ldev = np.float32(fields[12])
		udev = np.float32(fields[13])

		if atom1.upper() == 'CA' and atom2.upper() == 'CA':
			apt = "CaCa"
		elif atom1.upper() == 'CB' and atom2.upper() == 'CB' :
			apt = "CbCb"
		elif atom1.upper() == 'N' and atom2.upper() == 'O':
			apt = "NO"
		elif atom1.upper() == 'CG' and atom2.upper() == 'CG':
			apt = "CgCg"
		elif atom1.upper() == 'CA' and atom2.upper() == 'CG':
			apt = "CaCg"
		else:
			#print "WARNING: unsupported atom pair type in line: ", line
			continue

		if not restraints.has_key( apt ):
			restraints[apt] = - np.ones( (proteinLen, proteinLen, 3), dtype=np.float32 ) 

		restraints[apt][resid1-1, resid2-1,:] = np.array([dist, ldev, udev])
		restraints[apt][resid2-1, resid1-1,:] = np.array([dist, ldev, udev])

	return restraints


if len(sys.argv)<3:
	Usage()
	exit(1)

seqFile=sys.argv[1]
modelFile=sys.argv[2]
restraintFile=sys.argv[3]

seq=LoadFASTAFile(seqFile)
proteinLen = len(seq)

## atomDistMatrix is a dictionary, atomDistMatrix['CbCb'] is a matrix of Cb-Cb distance
atomDistMatrix = CalcAtomDistMatrix(modelFile)

## restraints is a dictionary, restraints['CbCb'] is a 3D arrary with the 3rd dimension contains information for estimated distance and lower, upper deviations
restraints = LoadDistRestraintFile(restraintFile, proteinLen)

## currently we only consider CbCb atoms

dist = atomDistMatrix['CaCa']

if not restraints.has_key('CaCa'):
	exit(0)

bounds = restraints['CaCa']

violation = abs(dist - bounds[:,:,0])

#violation = np.zeros( (proteinLen, proteinLen), dtype=np.float32 )

#np.putmask(violation, diff > bounds[:,:,2], diff - bounds[:,:,2])
#np.putmask(violation, (-diff) > bounds[:,:,1], -diff - bounds[:,:,1])
np.putmask(violation, bounds[:,:,0]<0, 0)

count = np.ones( (proteinLen, proteinLen), dtype=np.int32 )
np.putmask(count, bounds[:,:,0]<0, 0)

vioPerResidue = np.sum(violation, axis=1)
countPerResidue = np.sum(count, axis=1)
vioPerResidue = vioPerResidue / (countPerResidue + 0.00001)

## save the results
modelName=os.path.basename(modelFile).split('.')[0]
savefile = modelName + ".localQuality.txt"
fh = open(savefile, 'w')
fh.write('REMARK: Ca-Ca NOE violation per residue for ' + modelFile + '\n')
fh.write('REMARK: the 1st, 2nd and 3rd columns are residue index, amino acid and NOE violation, respectively\n')
for index, aa, vio in zip(range(proteinLen), seq, vioPerResidue):
	v = '%.2f\n' % vio
	fh.write(' '.join([ str(index+1), aa, v ]) )

fh.close()
