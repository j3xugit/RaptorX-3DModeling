import os
import sys
import subprocess

import numpy as np
import getopt
import shutil
import copy

#from Common import PDBUtils
from Common.SequenceUtils import LoadFASTAFile
#from ModellerUtils import Build3DModels
from Alignment.AlignmentUtils import ReadAlignment, ExtractSeqFromAlignment

## this script assemble domain models into whole-chain models using MODELLER

def Usage():
	print '	AssembleDomainAlignments.py [-d savefolder ] seqFile D1AlnFile D2AlnFile ...'
	print '	This script assembles domain alignments into a whole-chain alignment'
	print '	seqFile: the whole chain sequence file in FASTA format'
	print '	aliFiles: the alignment file for each domain in FASTA format, in which template is placed before query sequence'
	print '		There shall be no overlap between two domain sequences, otherwise it may not work correctly'

def PadSeqWithGaps(seq, numHeadGaps, overallLength):
	head = '-' * numHeadGaps
	tailLength = overallLength - len(seq) - numHeadGaps
	assert tailLength >= 0, 'ERROR: inconsistent domain length'
	tail = '-' * tailLength
	newSeq = head + seq + tail

	return newSeq

def SaveAln2FASTA(queryName, querySeq, tplSeqs, savefile):
	outStrs = []
	for tplName, tplSeq in tplSeqs:
		outStrs.append('>' + tplName)
		outStrs.append(tplSeq)
	outStrs.extend( [ '>'+queryName, querySeq] )
	finalStr = '\n'.join(outStrs)
	with open(savefile, 'w') as fh:
		fh.write(finalStr)
	return finalStr

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
		opts, args = getopt.getopt(sys.argv[1:], "d:", ["savefolder=="])
	except getopt.GetoptError as err:
		Usage()
		exit(1)

	savefolder = os.getcwd()
	for opt, arg in opts:
		if opt in ("-d", "--savefolder"):
			savefolder = arg
			if not os.path.isdir(savefolder):
				os.mkdir(savefolder)

		else:
			Usage()
			exit(1)

	if len(args) < 2:
		Usage()
		exit(1)

	seqFile = args[0]

	if not os.path.isfile(seqFile):
		print "ERROR, invalid sequence file: ", seqFile
		exit(1)
	mainSeq = LoadFASTAFile(seqFile)
	newSeq = copy.deepcopy(mainSeq)

	domainModels = args[1:]
	alignments = []
	for dmodel in domainModels:
		alignment, names = ReadAlignment(dmodel, returnNames=True)
		tplSeq, tgtSeq = ExtractSeqFromAlignment(alignment)
		index = newSeq.find(tgtSeq)
		if index < 0:
			print 'ERROR: cannot map the domain alignment sequence to the whole-chain alignment sequence'
			print 'domain seq: ', tgtSeq
			print 'whole  seq: ', newSeq
			exit(1)

		## alignment[0] and alignmentp1[1] are the template and query seq in alignment, respectively
		newSeq = newSeq.replace(tgtSeq, alignment[1])

		## names[0] is the template name
		alignments.append( (alignment, names[0]) )

	print newSeq

	newTplSeqs = []
	for aln, tplName in alignments:
		querySeq = aln[1]
		index = newSeq.find(querySeq)
		if index < 0:
			print 'ERROR: cannot map the domain alignment sequence to the whole-chain alignment sequence'
			print 'domain seq: ', querySeq
			print 'whole  seq: ', newSeq
			exit(1)
		newTplSeq = PadSeqWithGaps(aln[0], index, len(newSeq) )
		newTplSeqs.append((tplName, newTplSeq))

	print newTplSeqs
	

	## get names
	mainName = os.path.basename(seqFile).split('.')[0]
	tplNames = [ a for a, b in newTplSeqs ]

	## generate a FASTA file for the new alignment
	bname = '-'.join([mainName] + tplNames) 
	fastafile = bname + '.fasta'
	fastafile2 = os.path.join(savefolder, fastafile)	
	SaveAln2FASTA(mainName, newSeq, tplSeqs=newTplSeqs, savefile=fastafile2)
