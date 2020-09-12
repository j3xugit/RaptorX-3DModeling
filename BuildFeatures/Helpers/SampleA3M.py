import os
import sys
import numpy as np

## this script generates a number of samples of a multiple sequence alignment (in a3m format)
## the input MSA shall have the following format:
## >seq1_name
## AAAAAAAAAAAAAACCCCCCCCGGGGGGGGGGGGGGGGGGGGGNNNNNNNNNNNNNNNLLLLLLLLLLLLLLLL
## >seq2_name
##LLLLLLLLLLLLLLLLKKKKKKKKKKKKKKKNNNNNNNNNNNNNNNNGGGGGGGGGGGGCCCCCCCCCCCCCCCC

## That is, each sequence has two lines. The fist line is the sequence name and annotation and the 2nd line is the amino acid sequence
## further, the first sequence in the a3m file shall be the query sequence under prediction

def Usage():
	print 'python SampleA3M.py inputFile ratio [numSamples]'
	print '	inputFile: input MSA in a3m format'
	print '	ratio: the ratio of MSA to be sampled'
	print '	numSamples: number of sub-MSAs to be sampled'
	print ' the result files have names XXX.a3m_SY where XXX is the protein name and Y is the sample No. '

def GenerateOneSample(input, ratio):
	numSeqs = len(input)
	assert numSeqs%2 == 0

	## here we assume that the first 2 lines are the query sequence
	newMSA = input[:2]

	for i in np.arange(2, len(input), 2):
		r = np.random.uniform()
		if r < ratio * (1 + np.random.uniform(-0.05, 0.05) ):
			newMSA.extend(input[i:i+2])
	return newMSA
			
if len(sys.argv)<3:
	Usage()
	exit(1)

inputFile=sys.argv[1]
ratio=np.float32(sys.argv[2])

if ratio<=0 or ratio>1:
	print 'ERROR: the sampling rate is not in (0, 1]'
	exit(1)

numSamples = 1
if len(sys.argv) >= 4:
	numSamples = np.int32(sys.argv[3])

if numSamples<1:
	print 'ERROR: the number of sampled MSAs shall be at least 1'
	exit(1)

## load in the input MSA
with open(inputFile, 'r') as fh:
	input = list(fh)

name = os.path.basename(inputFile)

for sample in range(numSamples):
	oneMSA = GenerateOneSample(input, ratio)
	savefile = name + '_S' + str(sample)
	with open(savefile, 'w') as fh:
		fh.writelines(oneMSA)
