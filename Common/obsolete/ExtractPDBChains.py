import os
import sys
from PDBUtils import SaveOneChain2File

if len(sys.argv)<2:
	print 'python ExtractPDBChains.py chains [resDir] [pdbDir]'
	print '	chains: a list of chain names, separated by semicolon, each having format XXXXA or XXXXBB or XXXX_A or XXXX_BB where XXXX represents PDBID and A, B, AA, BB are chain letters'
	print '	the PDBID is not case sensitive'
	print '	resDir: the folder for result saving, default current work directory'
	print '	pdbDir: the folder contains the raw pdbfiles, pdb_ori by default'
	exit(1)

chains = sys.argv[1].split(';')

resDir=os.getcwd()
if len(sys.argv) >=3:
	resDir=sys.argv[2]

pdbDir = 'pdb_ori'
if len(sys.argv)>=4:
	pdbDir = sys.argv[3]


for chain in chains:
	if len(chain) < 4:
		print 'ERROR: incorrect pdb chain name'
		continue

	pdbid = chain[:4].lower()

	pdbfile1 = os.path.join(pdbDir, pdbid + '.pdb')
	pdbfile2 = os.path.join(pdbDir, pdbid + '.cif')
	if not os.path.isfile(pdbfile1) and (not os.path.isfile(pdbfile2) ):
		print 'ERROR: pdbfile or cif file does not exist for PDBID: ', pdbid
		continue
	elif os.path.isfile(pdbfile1):
		pdbfile = pdbfile1
	else:
		pdbfile = pdbfile2

	chainName = ''
	if len(chain)>4:
		if chain[4] == '_':
			chainName = chain[5:]
		else:
			chainName = chain[4:]

	savefile = os.path.join(resDir, chain+'.pdb')
	SaveOneChain2File(pdbfile, chainName, savefile)
