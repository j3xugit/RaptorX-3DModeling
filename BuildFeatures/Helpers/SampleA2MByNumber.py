import os
import sys
import numpy as np

## this script generates a number of samples of a multiple sequence alignment (in a2m format)
## the input MSA shall have the following format:
## AAAAAAAAAAAAAACCCCCCCCGGGGGGGGGGGGGGGGGGGGGNNNNNNNNNNNNNNNLLLLLLLLLLLLLLLL
## LLLLLLLLLLLLLLLLKKKKKKKKKKKKKKKNNNNNNNNNNNNNNNNGGGGGGGGGGGGCCCCCCCCCCCCCCC
## where each line represents one sequence

## That is, each sequence has two lines. The fist line is the sequence name and annotation and the 2nd line is the amino acid sequence
## further, the first sequence in the a3m file shall be the query sequence under prediction

def Usage():
	print 'python SampleA2MByNumber.py inputFile numSeqs [savefile]'
	print '	inputFile: input MSA in a2m format'
	print '	numSeqs: the number of sequences to be sampled'
	print '	savefile: the result file name'

def GenerateOneSample(input, numSampledSeqs):
	numLines = len(input)

	if numSampledSeqs<1 or numSampledSeqs>numLines:
		numSeqs = numLines
	else:
		numSeqs = numSampledSeqs

	## always keep the first line since it is the query sequence
	newMSA = [ input[0] ]
	numSeqs -= 1

	indices = np.arange(1, numLines)
	selected = np.random.choice(indices, size=numSeqs, replace=False) 
	selected = np.sort( selected )

	for i in selected:
		newMSA.append(input[i])
	return newMSA
			
if len(sys.argv)<3:
	Usage()
	exit(1)

inputFile=sys.argv[1]
numSeqs=np.int32(sys.argv[2])

if numSeqs<1:
	print 'ERROR: #sequences to be sampled shall be at least 1'
	exit(1)

name = os.path.basename(inputFile).split('.')[0]
savefile= name + '_sampled.a2m'
if len(sys.argv)>=4:
	savefile = sys.argv[3]

## load in the input MSA
with open(inputFile, 'r') as fh:
	input = list(fh)

oneMSA = GenerateOneSample(input, numSeqs)
with open(savefile, 'w') as fh:
	fh.writelines(oneMSA)
