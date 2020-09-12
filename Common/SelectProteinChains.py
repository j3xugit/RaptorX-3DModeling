import os
import sys
from Bio.PDB import PDBParser, MMCIFParser, Select
from Bio.PDB import PDBIO
from Bio.PDB.mmcifio import MMCIFIO

class ChainSelect(Select):
	def __init__(self, chainName):
		self.chainName = chainName

	def accept_chain(self, chain):
		if chain.get_id() == self.chainName:
			return 1
		else:
			return 0

def GetValidChains(structure, chainNames):
	available = []
	for model in structure:
		for chain in model:
			if chain.get_id() in chainNames:
				available.append(chain.get_id() )
	return available

if len(sys.argv) < 3:
	print 'python SelectProteinChains.py pdbfile/mmciffile chainNames [ResDir]'
	print '	this script may extract one or multiple protein structure chains from an input PDB or CIF file'
	print '	chainNames: one or multiple chain names separated by comma, e.g, A,B, 1,2, a, Aa,Ab, default empty'
	print '	one result file will be generated for each available chain'
	print '	The file has name PDBID_chainName.pdb or PDBID_chainName.cif where PDBID is the pdb code in upper case'
	exit(1)

infile=sys.argv[1]
if not os.path.isfile(infile):
	print 'ERROR: input structure file does not exist: ', infile
	exit(1)

if infile.endswith('.pdb'):
	parser = PDBParser()
	filesuffix = '.pdb'
elif infile.endswith('.cif'):
	parser = MMCIFParser()
	filesuffix = '.cif'
else:
	print 'ERROR: the input file shall end with either .pdb or .cif'
	exit(1)

chainNames = sys.argv[2].split(',')

ResDir=os.getcwd()
if len(sys.argv) > 3:
	ResDir = sys.argv[3]
	if not os.path.isdir(ResDir):
		os.mkdir(ResDir)

structure = parser.get_structure('test', infile)

availableChains = GetValidChains(structure, chainNames)
print '#requested chains: ', len(chainNames), ' #available chains: ', len(availableChains)

if infile.endswith('.cif'):
	io = MMCIFIO()
else:
	io = PDBIO()
io.set_structure(structure)

bname = os.path.basename(infile).split('.')[0]
for chainName in availableChains:
	savefile = bname.upper() + '_' + chainName + filesuffix
	savefile = os.path.join(ResDir, savefile)
	io.save(savefile, ChainSelect(chainName))
