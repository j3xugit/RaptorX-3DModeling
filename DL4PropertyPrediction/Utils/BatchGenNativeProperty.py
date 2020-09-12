import os
import sys
import cPickle

from GenNativeProperty import LoadTrainData4Properties

def Usage():
        print 'python BatchGenNativeProperty.py proteinList tplpklDir/nativeDir'
        print '	tplpklDir: the folder for tpl.pkl files or native.pkl files'

def main(argv):
        if len(argv) < 2:
                Usage()
                exit(1)

	proteinList = argv[0]
	tplDir = argv[1]

	with open(proteinList, 'r') as fp:
		proteins = [ p.strip() for p in list(fp) ]

	for protein in proteins:
		tplFile = os.path.join(tplDir, protein + '.tpl.pkl')
		if not os.path.isfile(tplFile):
			tplFile = os.path.join(tplDir, protein + '.native.pkl')
		property = LoadTrainData4Properties(tplFile)

		savefile = protein + '.nativeProperties.pkl'
        	with open(savefile, 'wb') as fh:
        		cPickle.dump( property, fh, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
        main(sys.argv[1:])
