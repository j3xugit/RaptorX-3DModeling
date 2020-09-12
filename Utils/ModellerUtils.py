import os
import sys
import getopt
import copy

from Bio.SeqIO.FastaIO import FastaIterator

# Homology modeling with multiple templates
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

from Common import PDBUtils

def WriteSeqToPIR(sequence, numLettersPerLine=60, end='*'):
        seq = sequence + end
        lines = []
        for i in range(0, len(seq), numLettersPerLine):
                lines.append( seq[i: min(i+60, len(seq)) ] )
        return lines

def SaveAlign2PIR(mainSeq, domainSeqs, mainSeqName, domainNames, chainIDs, savefile):
        lines = []

        ## write the mainSeq firt
        firstLine = '>P1;' + mainSeqName
        secondLine = 'sequence:' + mainSeqName + ':'*8
        lines.append(firstLine)
        lines.append(secondLine)
        lines.extend( WriteSeqToPIR(mainSeq))

        for dSeq, dName, cName in zip(domainSeqs, domainNames, chainIDs):
                lines.append(' ')
                lines.append('>P1;' + dName)
                lines.append('structureM:' + dName + '::' + cName + '::' + cName + ':'*4 )
                lines.extend( WriteSeqToPIR(dSeq) )

        with open(savefile, 'w') as fh:
                fh.write('\n'.join(lines) )

        return lines


## convert a fasta file to a pir file
## if seqName is not given, then the last record in the fasta file is treated as the query sequence and others treated as templates
def FASTA2PIR(fastafile, pirfile, PDBDir=None, seqName=None):

	querySeq = None
	queryName = seqName

	templateSeqs = []
	templateNames = []
	templateChainIDs = []

	records = []	
	with open(fastafile) as fh:
		for record in FastaIterator(fh):
			records.append( (record.name, str(record.seq) ) )

	def ProcessTemplate(r):
		templateSeqs.append(r[1])
		templateNames.append(r[0])
		id = ''
		if len(r[0])>=5:
			if r[0][4]=='_':
				id = r[0][5:]
			else:
				id = r[0][4:]

		templateChainIDs.append(id)

	if seqName is None:
		querySeq = records[-1][1]
		queryName = records[-1][0]
		for r in records[:-1]:
			ProcessTemplate(r)
	else:
		for r in records:
			if r[0] == seqName:
				queryName = seqName
				querySeq = r[1]
				continue
			ProcessTemplate(r)

	if PDBDir is not None:
		## match templateSeq to its PDB file
		revisedTemplateSeqs = []
		for tSeq, tName in zip(templateSeqs, templateNames):
			realSeq = tSeq.replace('-', '')
			seq2pdbMapping, _, _, _, _, _ = PDBUtils.MapSeq2PDB(realSeq, os.path.join(PDBDir, tName + '.pdb') )

			revisedSeq = []
			i=0
			for c in tSeq:
				if c == '-':
					revisedSeq.append('-')
					continue
				if seq2pdbMapping[i] < 0:
					revisedSeq.append('-')
				else:
					revisedSeq.append(c)

				i += 1
			revisedTemplateSeqs.append(''.join(revisedSeq))
		templateSeqs = revisedTemplateSeqs

	SaveAlign2PIR(querySeq, templateSeqs, queryName, templateNames, templateChainIDs, pirfile)
	return querySeq, templateSeqs, queryName, templateNames, templateChainIDs


def ExtractNames(PIRfile):
	## get seqName and templateNames from PIRfile
	with open(PIRfile, 'r') as fh:
		content = [ line.strip() for line in list(fh) ]

	templateNames = []
	for line in content:
		if line.startswith('sequence'):
			seqName = line.split(':')[1]
			continue
		if line.startswith('structure'):
			tempName = line.split(':')[1]
			templateNames.append(tempName)

	return seqName, templateNames

def Build3DModels(PIRfile, PDBDir, seqName=None, templateNames=None, maxCaCaDist=30, numModels=1):
	if seqName is None and templateNames is None:
		sName, tNames = ExtractNames(PIRfile)
	elif seqName is None:
		sName, _ = ExtractNames(PIRfile)
		tNames = templateNames
	elif templateNames is None:
		_, tNames = ExtractNames(PIRfile)
		sName = seqName
	else:
		sName, tNames = seqName, templateNames

	#log.verbose()    # request verbose output
	log.minimal()
	env = environ()  # create a new MODELLER environment to build this model in

	# directories for input atom files
	if isinstance(PDBDir, list):
		env.io.atom_files_directory = PDBDir
	else:
		env.io.atom_files_directory = [PDBDir]

	a = automodel(env,
              alnfile  = PIRfile, # alignment filename
              knowns   = tuple(tNames),     # codes of the templates
              sequence = sName)               # code of the target

	a.max_ca_ca_distance = maxCaCaDist           # long distance restraints
	a.starting_model= 1                 # index of the first model
	a.ending_model  = numModels                # index of the last model
	a.make()                            # do the actual homology modeling

def Usage():
	print 'python ModellerUtils.py [-q seqName] PIRfile/FASTAfile templatePDBDir'


if __name__ == "__main__":

	queryName = None

	try:
                opts, args = getopt.getopt(sys.argv[1:], "q:", ["queryName=="])
        except getopt.GetoptError as err:
                Usage()
                exit(1)

        for opt, arg in opts:
                if opt in ("-q", "--queryName"):
                        queryName = arg
                else:
                        Usage()
                        exit(1)

        if len(args) < 2:
                Usage()
                exit(1)


	inputfile = args[0]
	PDBDir = args[1]

	if inputfile.endswith('.pir'):
		Build3DModels(inputfile, PDBDir)
	else:
		pirfile = os.path.basename(inputfile).replace('.fasta', '.pir')
		querySeq, templateSeqs, queryName, templateNames, _ = FASTA2PIR(inputfile, pirfile, PDBDir)
		Build3DModels(pirfile, PDBDir, seqName=queryName, templateNames=templateNames)
	
