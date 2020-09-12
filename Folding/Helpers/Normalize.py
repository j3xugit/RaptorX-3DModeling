import os
import sys
import numpy as np

## this script normalizes all columns except the first in a file

def Usage():
	print 'python Normalize.py textfile [savefile]'

if len(sys.argv) < 2:
	Usage()
	exit(1)

infile=sys.argv[1]

with open(infile, 'r') as fh:
	content = [ line.strip() for line in list(fh) ]

modelFiles = []
scores = []
for row in content:
	fields = row.split()
	assert len(fields) >= 3
	modelFiles.append(fields[0])
	s = [ np.float32(e) for e in fields[1:] ]
	scores.append(s)

scores = np.array(scores)

scores_m = np.mean(scores, axis=0, keepdims=True)
scores_std = np.std(scores, axis=0, keepdims=True)

#print 'std of ', infile, scores_std

scores_normalized = (scores - scores_m)/(scores_std + 0.000001)

## print
new_scores = np.concatenate( (scores_normalized, scores), axis=1)

scoreStrs = []
for m, score in zip(modelFiles, new_scores):
	scoreStr =[m] +  [ str(e) for e in score ]
	scoreStr = '\t'.join(scoreStr)
	scoreStrs.append(scoreStr)

if len(sys.argv)>=3:
	outfile = sys.argv[2]
	with open(outfile, 'w') as fh:
		fh.writelines('\n'.join(scoreStrs))
else:
	print '\n'.join(scoreStrs)
