import os
import sys
import numpy as np
from subprocess import call 

DistDir="Dist_CASP12DMTPL_TPLEC52CSet7952DTAtom15Models/"
AngleDir="PhiPsiSS8_CASP12DMTPL_TPLSet7952Models/"

def Usage():
	print 'python GenMultiTemplateRestraint4OneProtein.py listFile'
	

def GenRestraintFile4OneTarget(target, templates):
	cmd = ['python', os.path.join(os.environ['ModelingHome'],'DL4DistancePrediction2/TPLMergePredictedDistMatrix.py')]
	for temp in templates:
		cmd.append(os.path.join(DistDir, target+'-'+temp+'.predictedDistMatrix.pkl'))
	call(cmd)

	cmd = ['python', os.path.join(os.environ['ModelingHome'],'DL4PropertyPrediction/TPLMergePredictedProperties.py')]
	for temp in templates:
		cmd.append(os.path.join(AngleDir, target+'-'+temp+'.predictedProperties.pkl'))
	call(cmd)

	cmd = ['cp', 'CASP12DM-SEQ/'+target+'.fasta', '-'.join([target] + templates) + '.fasta' ]
	call(cmd)

if len(sys.argv)<2:
	Usage()
	exit(-1)
	
listFile = sys.argv[1]


fh = open(listFile, 'r')
content = [ line.strip() for line in list(fh) ]
fh.close()

for line in content:
	fields = line.split()
	target = fields[0]
	templates = fields[1:]
	GenRestraintFile4OneTarget(target, templates)
