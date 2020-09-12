import os
import sys
import numpy as np

def Usage():
	print 'python AnalyzeScoreNQuality.py scorefile qualityfile'
	
if len(sys.argv)<3:
	Usage()
	exit(1)

scorefile=sys.argv[1]
qualityfile=sys.argv[2]

with open(scorefile, 'r') as fh:
	scores = [ line.strip() for line in list(fh) ]

modelScores = dict()

for score in scores:
	fields = score.split()
	assert len(fields) == 9
	model = os.path.basename(fields[0])
	s = [ np.float32(f) for f in fields[1:] ]
	modelScores[model] = s

with open(qualityfile, 'r') as fh:
	quality = [ line.strip() for line in list(fh) ]

modelQuality = dict()
for q in quality:
	fields = q.split()
	assert len(fields) == 7
	model = fields[0] + '.pdb'
	qq = [ np.float32(f) for f in fields[3:] ]
	modelQuality[model] = qq

scores = []
quality = []
## calculate correlation between quality and scores

for m, q in modelQuality.iteritems():
	if modelScores.has_key(m):
		quality.append(q)
		scores.append(modelScores[m][:4])

scores=np.array(scores)
quality=np.array(quality)

corr = np.corrcoef(quality, y=scores, rowvar=False)
print corr[0:4, 4:]
