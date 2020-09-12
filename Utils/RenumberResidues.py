import os
import sys

def Usage():
	print 'python RenumberResidues.py pdbfile'
	print '	This script renumbers the residues of a PDB file'

if len(sys.argv) < 2:
	Usage()
	exit(1)

pdbfile = sys.argv[1]

with open(pdbfile, 'r') as fh:
	content = [ line.strip() for line in list(fh) ]

newContent = []
resNum = 0
prevID = None

for row in content:
	if not (row.startswith('ATOM') or row.startswith('TER')):
		newContent.append(row)
		continue
	
	currID = row[17:26]
	if prevID is None or currID != prevID:
		## increase the residue number by 1
		resNum += 1
		prevID = currID
	newRow = row[:22] + "%4d" % (resNum) + row[26:]
	newContent.append(newRow)

## print out
#savefile = os.path.basename(pdbfile).replace('.pdb', '.revised.pdb')
print '\n'.join(newContent)

