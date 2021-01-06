import sys
import os
import tempfile
import cPickle

from GenPropertyFeaturesFromHHM import LoadTestData4Properties

## this script generates input features for multiple proteins from HHM files for property prediction

def Usage():
	print 'python GenPropertyFeatures4Proteins.py proteinList hhmFolder/a3mFolder [savefile]'
	print 'This script generates a property feature file (in PKL) for multiple proteins from their MSA files'
	print '\tproteinList: a file containing a list of protein names, each in one line'
	print '\thhmFolder/a3mFolder: a folder containing the MSA files (i.e, .hhm or .a3m files); when .a3m files are used as input, please make sure that hhmake is accessible'
	print '\tif not specified, the result file is named after proteinListFileName.propertyFeatures.pkl'

def main(argv):

	if len(argv) < 2:
		Usage()
		exit(1)

	proteinNameFile = argv[0]
	hhmFolder = argv[1]
	if not os.path.isdir(hhmFolder):
		print 'ERROR: the specified folder for .hhm files does not exist: ', hhmFolder
		exit(1)

	savefile = os.path.basename(proteinNameFile) + '.propertyFeatures.pkl'
	if len(argv) > 2:
		savefile = argv[2]

	with open(proteinNameFile, 'r') as fh:
		proteins = [ line.strip() for line in list(fh) ]

	features = []
	for protein in proteins:
		hhmFileIsTemp = False 
		hhmFile = os.path.join(hhmFolder, protein + '.hhm')
		if not os.path.isfile(hhmFile):
			a3mFile = os.path.join(hhmFolder, protein + '.a3m')
			if not os.path.isfile(a3mFile):
				print 'ERROR: cannot find hhmFile or a3mFile for', protein, 'in', hhmFolder
				continue
			hhmFileIsTemp = True
			hhmFh, hhmFile = tempfile.mkstemp(prefix=protein, suffix='.hhm')
			os.close(hhmFh)
			cmdStr = 'hhmake -i ' + a3mFile + ' -o ' + hhmFile
			os.system(cmdStr)
		
		feature = LoadTestData4Properties(hhmFile, protein)
		if hhmFileIsTemp:
			os.remove(hhmFile)
		if feature is None:
			print 'ERROR: failed to load hhm file: ', hhmFile
			continue
		features.append(feature)

	with open(savefile, 'wb') as fh:
		cPickle.dump(features, fh, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
        main(sys.argv[1:])

