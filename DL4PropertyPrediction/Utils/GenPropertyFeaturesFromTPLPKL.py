import numpy as np
import sys
import os
import cPickle

from DL4PropertyPrediction import PropertyUtils

## this script generates input features and labels for property prediction, training and validation.
## it takes information from a tplpkl file

def Usage():
	print 'python GenPropertyFeaturesFromTPLPKL.py tplpklFile'
	print '	This script generates property features and labels from a tpklpkl file'
	print '	tplpklFile: the tplpkl file generated for template-based modeling'

def LoadTrainData4Properties(tplFile):
	with open(tplFile) as fh:
		tpl = cPickle.load(fh)

	protein = dict()
	protein['name'] = tpl['name']
	protein['sequence'] = tpl['sequence']
	protein['length'] = tpl['length']
	protein['NEFF'] = tpl['NEFF']

	protein['PSFM'] = tpl['PSFM']
	protein['PSSM'] = tpl['PSSM']

	protein['ACC'] = tpl['ACC']
	protein['pACC'] = tpl['pACC']
	protein['CNa'] = tpl['CNa']	
	protein['CNb'] = tpl['CNb']

	protein['SS'] = tpl['SS']
	#protein['DISO'] = Ang['DISO']
        #protein['CLE'] = Ang['CLE']
	protein['Phi'] = tpl['Phi']
	protein['Psi'] = tpl['Psi']
	protein['Theta'] = tpl['Theta']
	protein['Tau'] = tpl['Tau']

	if tpl.has_key('Omg'):
		protein['Omg'] = tpl['Omg']

	##merge Phi and Psi
	protein['PhiPsi'] = np.transpose( np.array([ protein['Phi'], protein['Psi'] ]) )

	##merge Theta and Tau
	protein['ThetaTau'] = np.transpose( np.array([ protein['Theta'], protein['Tau'] ]) )

	## the missing residues have no 3D coordinates and thus, angles and solvent accessibility
	protein['Missing'] = tpl['missing']
	protein['DISO'] = tpl['missing']
	protein['SS8'] = protein['SS']
	protein['SS3'] = ''.join([PropertyUtils.SS8Letter2SS3Letter[c] for c in protein['SS8'] ] )

	protein['ForTrain'] = True

	return protein

def main(argv):
	if len(argv) < 1:
		Usage()
		exit(1)

	tplFile = argv[0]
	if not tplFile.endswith('.tpl.pkl'):
		print 'ERROR: the input file shall end with .tpl.pkl'
		exit(1)

	protein = LoadTrainData4Properties(tplFile)

	savefile = protein['name'] + '.propertyFeatures.pkl'
	with open(savefile, 'wb') as fh:
		cPickle.dump(protein, fh, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
        main(sys.argv[1:])

