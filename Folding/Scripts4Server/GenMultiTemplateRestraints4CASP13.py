import os
import sys
import numpy as np
import socket
from subprocess import call 

OldDistDir="Dist_CASP13_TPLEC52CSet7952DTAtom15Models/"
OldAngleDir="PhiPsiSS8_CASP13_TPLSet7952Models/"

DistDirTop20="Dist_CASP13_TPLEC52CSet10820Top20DTAtom10Models/"
DistDirTop40="Dist_CASP13_TPLEC52CSet10820Top40DTAtom10Models/"

AngDirTop40="PhiPsiSS8_CASP13_TPLSet10820Top30Models/"
AngDirTop20="PhiPsiSS8_CASP13_TPLSet10820Top20Models/"


SeqDir="/mnt/data/CASP13/SEQ/"
TPLSeqDir="/mnt/data/CASP13/TPL-SEQ/"

DistDirPrefix="/mnt/data/CASP13/PredictedDistMatrix/"
AngleDirPrefix="/mnt/data/CASP13/PredictedProperties/"

def Usage():
	print 'python GenMultiTemplateRestraint4CASP13.py targetName Top20/Top40 template1 template2 ... '
	

def GenRestraintFile4OneTarget(target, templates, model="Top20"):
	if model == 'Top20':
		DistDir = DistDirTop20
		AngleDir = AngDirTop20
	elif model == 'Top40':
		DistDir = DistDirTop40
		AngleDir = AngDirTop40
	else:
		print 'unsupported model in GenRestraintFile4OneTarget: ', model
		exit(-1)

	DistDir = os.path.join(DistDirPrefix, DistDir)
	AngleDir = os.path.join(AngleDirPrefix, AngleDir)

	cmd = ['python', os.path.join(os.environ['ModelingHome'],'DL4DistancePrediction2/TPLMergePredictedDistMatrix.py')]
	for temp in templates:
		cmd.append(os.path.join(DistDir, target+'-'+temp+'.predictedDistMatrix.pkl'))
	call(cmd)

	cmd = ['python', os.path.join(os.environ['ModelingHome'],'DL4PropertyPrediction/TPLMergePredictedProperties.py')]
	for temp in templates:
		cmd.append(os.path.join(AngleDir, target+'-'+temp+'.predictedProperties.pkl'))
	call(cmd)

	machine = socket.gethostname()
	account = ''.join(['RaptorX', '@', machine, ':'])
	cmd = ['scp', os.path.join(SeqDir, target+'.fasta'), account + os.path.join(TPLSeqDir, '-'.join([target] + templates) + '.fasta') ]
	call(cmd)
	cmd = ['ssh', account[:-1], 'chmod', 'a+r', os.path.join(TPLSeqDir, '-'.join([target] + templates) + '.fasta') ]
	call(cmd)

	DistDestDir=DistDir
	cmd = ['chmod', 'a+r', '-'.join([target] + templates) + '.predictedDistMatrix.pkl']
	call(cmd)
	print 'Copying the predicted dist matrix file to a shared directory ...'
	cmd = [ "scp", '-'.join([target] + templates) + '.predictedDistMatrix.pkl', account + DistDestDir ]
	call(cmd)

	AngleDestDir=AngleDir	
	cmd = ['chmod', 'a+r', '-'.join([target] + templates) + '.predictedProperties.pkl']
	call(cmd)
	print 'Copying the predicted property file to a shared directory ...'
	cmd = [ "scp", '-'.join([target] + templates) + '.predictedProperties.pkl', account + AngleDestDir ]
	call(cmd)


if len(sys.argv)<4:
	Usage()
	exit(-1)
	
target = sys.argv[1]
model = sys.argv[2]
templates = sys.argv[3:]

GenRestraintFile4OneTarget(target, templates, model)
