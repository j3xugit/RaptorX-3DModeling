import os
import sys

from Bio.PDB import PDBList

##this script may download pdb files for a list of pdb codes.

if len(sys.argv)<3:
	print 'python DownloadPDBs.py pdbcodeListFile format [ResDir]'
	print '\tpdbcodeListFile: a file of PDB IDs, in upper or lower case'
	print '\tformat: pdb or mmCif. This option specifies which PDB file format needed'
	exit(1)

listFile=sys.argv[1]
with open(listFile, 'r') as fh:
	content = [ line.strip() for line in list(fh) ]

pdbcodes = []
for c in content:
	pdbcodes.extend( c.split() )

format = sys.argv[2]

ResDir=os.getcwd()
if len(sys.argv)>3:
	ResDir=sys.argv[3]

if not os.path.isdir(ResDir):
	os.mkdir(ResDir)

pdblist = PDBList()
pdblist.download_pdb_files(pdbcodes, file_format=format, pdir=ResDir)
