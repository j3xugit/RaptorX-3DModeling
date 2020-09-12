import os
import sys

from AlignmentUtils import ParseAlignmentFile

if __name__ == '__main__':
        if len(sys.argv) < 2:
                print 'python ParseAlignmentFile.py alnFile [savefolder]'
                exit(1)

        alnFile = sys.argv[1]
	savefolder = os.getcwd()
	if len(sys.argv) >= 3:
        	savefolder = sys.argv[2]
		if not os.path.isdir(savefolder):
			os.mkdir(savefolder)

	ParseAlignmentFile(alnFile, savefolder)

