import os
import sys
import subprocess

import numpy as np
import getopt
import shutil

from Common import PDBUtils
from Common.SequenceUtils import LoadFASTAFile
from ModellerUtils import Build3DModels

## this script assemble domain models into whole-chain models using MODELLER

def Usage():
	print '	AssembleDomainModels.py [-d savefolder | -l lengthOfLinkerRegion ] seqFile D1ModelFile D2ModelFile ...'
	print '	This script assembles domain models into a whole-chain model using MODELLER'
	print '	seqFile: the whole chain sequence file in FASTA format'
	print '	modelFiles: the model file for each domain'
	print '	-l: a small number of ALA residues to be added between two domains to reduce steric clashes, default 0'

## here we assume that the local quality (distance deviation) occpuies columns 61-66 in the ATOM record of a PDB file
def ExtractLocalQuality(pdbfile):
	with open(pdbfile, 'r') as fh:
		CAquality = [ np.float32(line[60:66]) for line in list(fh) if line.startswith('ATOM ') and line[12:16] == ' CA ' ]
	return CAquality

def PadSeqWithGaps(seq, numHeadGaps, overallLength):
	head = '-' * numHeadGaps
	tailLength = overallLength - len(seq) - numHeadGaps
	assert tailLength >= 0, 'ERROR: inconsistent domain length'
	tail = '-' * tailLength
	newSeq = head + seq + tail

	return newSeq

def WriteSeqToPIR(sequence, numLettersPerLine=60, end='*'):
	seq = sequence + end
	lines = []
	for i in range(0, len(seq), numLettersPerLine):
		lines.append( seq[i: min(i+60, len(seq)) ] )
	return lines

def SaveAlign2PIR(mainSeq, domainSeqs, mainSeqName, domainNames, chainIDs, savefile):
	lines = []

	## write the mainSeq first
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

def AddLinkers(mainSeq, domainFlag, linkerlength):

	## find where all linker regions start
	linkerStarts = []
	for i in range(len(domainFlag)-1):
		curr = domainFlag[i]
		next = domainFlag[i+1]

		if curr>0 and next>0 and curr != next:
			## find a linker position
			linkerStarts.append(i)

	## split the mainSeq into segments
	segments = []
	prevPos = -1
	for pos in linkerStarts:
		seg = mainSeq[prevPos+1 : pos+1]
		segments.append(seg)
		prevPos = pos
	if prevPos < len(mainSeq)-1:
		segments.append( mainSeq[prevPos+1:] )

	## add linkers between two segments
	linker = 'A' * linkerlength
	newSeq = linker.join(segments)
	
	return newSeq	

if __name__ == "__main__":

	if len(sys.argv) < 3:
		Usage()
		exit(1)
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "d:l:", ["savefolder==", "linkerlength=="])
	except getopt.GetoptError as err:
		Usage()
		exit(1)

	linkerlength = 0

	savefolder = os.getcwd()
	for opt, arg in opts:
		if opt in ("-d", "--savefolder"):
			savefolder = arg
			if not os.path.isdir(savefolder):
				os.mkdir(savefolder)

		elif opt in ("-l", "--linkerlength"):
			linkerlength=np.int32(arg)
			assert linkerlength >= 0

		else:
			Usage()
			exit(1)

	if len(args) < 2:
		Usage()
		exit(1)

	seqFile = args[0]
	domainModels = args[1:]

	if not os.path.isfile(seqFile):
		print "ERROR, invalid sequence file: ", seqFile
		exit(1)
	mainSeq = LoadFASTAFile(seqFile)

	domainSeqs = []
	chainIDs = []
	localQualitys = []
	for dmodel in domainModels:
		pdbseqs, _, chains = PDBUtils.ExtractSeqFromPDBFile(dmodel)
		assert len(pdbseqs) == 1
		domainSeqs.append( pdbseqs[0] )
		#print 'chain id: ', chains[0].get_id()
		chainIDs.append( chains[0].get_id() )
		locQuality = ExtractLocalQuality(dmodel)
		assert len(locQuality) == len(pdbseqs[0])
		localQualitys.append( locQuality )

	print mainSeq
	print domainSeqs

	## initial alignment between mainSeq and domainSeqs
	domainFlag = [ 0 ] * len(mainSeq)
	for dSeq, dNum in zip(domainSeqs, range(1, 1+len(domainSeqs) ) ):
		index = mainSeq.find(dSeq)
		if index < 0:
			print 'ERROR: cannot map the domain sequence to the whole-chain sequence'
			print 'domain seq: ', dSeq
			print 'whole  seq: ', mainSeq
			exit(1)
		domainFlag[index : index+len(dSeq)] = [ dNum ] * len(dSeq)
	print domainFlag

	realMainSeq = mainSeq
	if linkerlength > 0:
		newMainSeq = AddLinkers(mainSeq, domainFlag, linkerlength)
		mainSeq = newMainSeq
	
	## align domainSeqs to possibly new mainSeq
	paddedDomainSeqs = []
	quality = [ 99.0 ] * len(mainSeq)
	for dSeq, locQuality in zip(domainSeqs, localQualitys):
		index = mainSeq.find(dSeq)
		if index < 0:
			print 'ERROR: cannot map the domain sequence to the whole-chain sequence'
			print 'domain seq: ', dSeq
			print 'whole  seq: ', mainSeq
			exit(1)
		newDomainSeq = PadSeqWithGaps(dSeq, index, len(mainSeq) )
		paddedDomainSeqs.append(newDomainSeq)
		quality[index: index + len(dSeq) ] = locQuality

	## test
	print mainSeq
	for subSeq in paddedDomainSeqs:
		print subSeq

	## get names
	mainName = os.path.basename(seqFile).split('.')[0]

	## copy the domain models to the same folder
	domainNames = []
	for dmodel in domainModels:
		bname = '.'.join(os.path.basename(dmodel).split('.')[:-1])
		domainNames.append(bname)
		newModelFile = os.path.join(savefolder, bname+'.pdb')
		if not os.path.isfile(newModelFile):
			os.link(dmodel, newModelFile )

	## output CA quality
	qualityfileBaseName = '-'.join( ['CA.quality', mainName] + domainNames ) + '.txt'
	qualityfile = os.path.join(savefolder, qualityfileBaseName)
	qualityStr = [ '{:6.2}'.format(q) for q in quality ] 
	with open(qualityfile, 'w') as fh:
		fh.write('\n'.join(qualityStr) )

	## generate a PIR file
	bname = '-'.join([mainName] + domainNames) 
	pirfile = bname + '.pir'
	pirfile2 = os.path.join(savefolder, pirfile)	
	SaveAlign2PIR(mainSeq, domainSeqs=paddedDomainSeqs, mainSeqName=bname, domainNames=domainNames, chainIDs=chainIDs, savefile=pirfile2)

	scriptdir = os.path.dirname(os.path.realpath(__file__))

	## build 3D models for the whole chain; savefolder has the template PDB files
	currDir = os.getcwd()
	os.chdir(savefolder)
	numModels = 2
	Build3DModels(pirfile, './', seqName=bname, templateNames=domainNames, maxCaCaDist=30, numModels=numModels)

	## add quality to models
	for i in range(1, 1+numModels):
		rawModelFile = bname + '.B9999XXXX.pdb'.replace('XXXX', "{:04d}".format(i))
		newModelFile = bname + '.' + "{:02d}".format(i) + '.pdb'
		addQualityCmd = [os.path.join(scriptdir, 'AddErrorEstimate2PDB'), rawModelFile, qualityfileBaseName, newModelFile]
		#subprocess.run(addQualityCmd)
		os.system(' '.join(addQualityCmd) )

	"""
	pdbfile = bname + '.pdb' 
	os.rename(bname + '.B99990001.pdb', pdbfile)
	"""

	## go back 
	os.chdir(currDir)
