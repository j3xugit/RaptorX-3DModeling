import numpy as np
import os

import sys

if len(sys.argv)<2:
	print 'python FindBest.py qualityFile4OneTarget'
	exit(1)

qualityFile = sys.argv[1]
names = os.path.basename(qualityFile).split('-')[:2]
name = '-'.join(names)

with open(qualityFile, 'r') as file:
	content = file.readlines()
	scores = []
	for row in content:
		fields = row.split()
		assert len(fields) == 7
		model = fields[0].replace('_result', '')
		length = fields[1]
		TM = np.float32(fields[3])
		GDT = np.float32(fields[5])
		GHA = np.float32(fields[6])
		scores.append( [model, TM, GDT, GHA] )

	if len(scores) < 1:
		exit(0)

	sortedTM = sorted(scores, key=lambda s: s[1], reverse=True)
	sortedGDT = sorted(scores, key=lambda s: s[2], reverse=True)
	sortedGHA = sorted(scores, key=lambda s: s[3], reverse=True)
	sortedTMGDT = sorted(scores, key=lambda s: (s[1] + s[2])/2, reverse=True )

	nameStr = [ str(e) for e in [os.path.dirname(qualityFile).split('/')[-1][:25], name, length] ]
	bestScores = [sortedTM[0][1], sortedGDT[0][2], (sortedTMGDT[0][1] + sortedTMGDT[0][2])/2, sortedGHA[0][3] ]
	bestStr = [ "{0:.4f}".format(e) for e in bestScores ]
	
	topK = len(scores)/10
	if topK < 1:
		topK = 1

	topTM = np.mean([ tm[1] for tm in sortedTM[:topK] ])
	topGDT = np.mean([ gdt[2] for gdt in sortedGDT[:topK] ])
	topGHA = np.mean([ gha[3] for gha in sortedGHA[:topK] ])
	topTMGDT = np.mean([ (s[1]+s[2])/2 for s in sortedTMGDT[:topK] ])

	topScores = [ topTM, topGDT, topTMGDT, topGHA ]
	topStr = [ "{0:.4f}".format(e) for e in topScores ]

	resultStr = nameStr + bestStr + ['top'+str(topK)] + topStr

	finalStr = '\t'.join(resultStr)
	print finalStr
