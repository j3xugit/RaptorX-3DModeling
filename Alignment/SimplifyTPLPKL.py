import os 
import sys
import cPickle


def Usage():
	print 'python SimplifyTPLPKL.py infile_TPLPKL [outDir]'
	print '	This script simplifies a heavy-weight tpl.pkl file to a lightweight tpl.pkl file for threading'
	print '	The heavy weight file is only needed for template-based distance/orientation prediction'
	print '	infile: the input file ending with .tpl.pkl'
	print '	outDir (optional): the folder for result saving, default current work directory'
	print '	The result file has the same name as input file, so make sure that you do not overwrite the input file'

if len(sys.argv) < 2:
	Usage()
	exit(1)

infile = sys.argv[1]

outDir = os.getcwd()
if len(sys.argv) > 2:
	outDir = sys.argv[2]
	if not os.path.isdir(outDir):
		os.path.mkdir(outDir)

with open(infile, 'rb') as fh:
	tpl = cPickle.load(fh)

newtpl = dict()

for k, v in tpl.iteritems():
	if k == 'atomDistMatrix':
		newtpl['atomDistMatrix'] = dict()
		newtpl['atomDistMatrix']['CbCb'] = v['CbCb']
		continue

	elif k == 'atomOrientationMatrix':
		continue

	newtpl[k] = v

savefile = os.path.join(outDir, os.path.basename(infile) )
with open(savefile, 'wb') as fh:
	cPickle.dump( newtpl, fh, protocol=cPickle.HIGHEST_PROTOCOL)
