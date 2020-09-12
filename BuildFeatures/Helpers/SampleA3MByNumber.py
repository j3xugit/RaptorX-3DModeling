import os
import sys
import numpy as np

## this script samples an subMSA from from an MSA (in a3m format)
## the input MSA shall have the following format:
## >seq1_name
## AAAAAAAAAAAAAACCCCCCCCGGGGGGGGGGGGGGGGGGGGGNNNNNNNNNNNNNNNLLLLLLLLLLLLLLLL
## >seq2_name
##LLLLLLLLLLLLLLLLKKKKKKKKKKKKKKKNNNNNNNNNNNNNNNNGGGGGGGGGGGGCCCCCCCCCCCCCCCC

## That is, each sequence has two lines. The fist line is the sequence name and annotation and the 2nd line is the amino acid sequence
## further, the first sequence in the a3m file shall be the query sequence under prediction

def Usage():
	print 'python SampleA3MByNumber.py inputFile numSeqs [savefile]'
	print '	inputFile: an MSA file in a3m format'
	print '	numSeqs: the number of sequences to be sampled. If <0 or > #seqs in input MSA, then output the input MSA as the result'
	print '	savefile: the file for result saving'

def GenerateOneSample(input, numSampledSeqs):
	numLines = len(input)
	assert numLines % 2 == 0

	if numSampledSeqs<1 or numSampledSeqs>(numLines/2):
		numSeqs = numLines/2
	else:
		numSeqs = numSampledSeqs

	## always keep the first 2 lines since they are the query sequence
	newMSA = input[:2]
	numSeqs -= 1

	indices = np.arange(2, numLines, 2)
	selected = np.random.choice(indices, size=numSeqs, replace=False) 
	selected = np.sort( selected )

	for i in selected:
		if input[i].startswith('ss_pred') or input[i].startswith('ss_conf'):
			continue
		newMSA.extend(input[i:i+2])
	return newMSA
			
if len(sys.argv)<3:
	Usage()
	exit(1)

inputFile=sys.argv[1]
numSeqs=np.int32(sys.argv[2])

if numSeqs<1:
	print 'ERROR: the # of sequences to be sampled shall be at least 1'
	exit(1)

name = os.path.basename(inputFile).split('.')[0]
savefile= name + '_sampled.a3m'
if len(sys.argv)>=4:
	savefile = sys.argv[3]

## load in the input MSA
with open(inputFile, 'r') as fh:
	input = list(fh)

oneMSA = GenerateOneSample(input, numSeqs)
with open(savefile, 'w') as fh:
	fh.writelines(oneMSA)

print 'The resultant file is saved to', savefile
