import os
import sys
import numpy as np

def Usage():
	print 'python CollectTPLNOEScore.py TPLNOEfile [DeepThreaderRankFile]'
	

def LoadTemplates(tempListFile):
	fh = open(tempListFile)
	templates = [ line.strip() for line in list(fh) ]
	fh.close()
	return set(templates)

def LoadFASTAFile(seqFile):

        ##load in the seq file
        seqfh = open(seqFile, 'r')
        content = [ line.strip() for line in list(seqfh) ]
        seqfh.close()
        content2 = [ c for c in content if not c.startswith('>') and not c.startswith('#') ]
        sequence = ''.join(content2)

        return sequence

def LoadRankFile(rankFile):
	fh = open(rankFile, 'r')
	content = list(fh)
	fh.close()

	## seq lenght
	assert content[17].startswith('Query Length')
	seqLen = np.int32(content[17].split('=')[1])
	##print 'seqLen=', seqLen

	## get template and alignment length
	assert content[26].startswith('No ')
	assert content[27].startswith('1 ')

	results = dict()
	for row in content[27:1027]:
		fields = row.split()
		assert len(fields) == 13
		template = fields[1]
		alignLen = np.int32(fields[9])

		results[template] = alignLen

	return results, seqLen

def LoadNOEFile(NOEfile):
	fh = open(NOEfile, 'r')
	content = [ line.strip() for line in list(fh) ]
	fh.close()

	## parse content
	assert len(content) % 6 ==0

	#print len(content)

	data = []
	for i in range(0, len(content), 6):
		modelFile = content[i]
		modelFlag = 'EC52C'

		if not modelFlag in modelFile:
			modelFlag = 'EC25C'

		#print modelFile
		## find template
		fields = modelFile.split('/')[-2].split('_')[0].split('-')
		if len(fields) != 2:
			print "Skip ", modelFile
			continue

		template = fields[1]

		## read the 5 lines of NOE scores
		NOEs = []
		for j in range(1, 6):
			line = content[i+j]
			fields = line.split()
			assert len(fields) == 3
			avgNOE = np.float32(fields[2])
			NOEs.append(avgNOE)
		noe = np.average(NOEs)

		data.append( [template, noe, modelFlag == 'EC52C' ] )

	return data

def LoadAlignModelFile(alignModelFile):
	fh = open(alignModelFile, 'r')
	content = [ line.strip() for line in list(fh) ]
	fh.close()

	data = dict()
	GDTs = []
	for line in content:
		fields = line.split()
		assert len(fields) == 13
		template = fields[0]
		seqLen = np.int32(fields[3])
		alignLen = np.int32(fields[4])
		
		data[template] = [seqLen, alignLen]
		GDTs.append(np.float32(fields[-2]) )

	bestOf1 = GDTs[0]
	bestOf5 = max(GDTs[:5])
	bestOf20 = max(GDTs[:20])
	bestOf100 = max(GDTs[:100])
	bestOfAll = max(GDTs)

	bests = [ bestOf1, bestOf5, bestOf20, bestOf100, bestOfAll ]

	return data, bests

if len(sys.argv)<2:
	Usage()
	exit(-1)

NOEfile=sys.argv[1]

rankFile=None
ranks = None
seqLen = None
if len(sys.argv)>=3:
	rankFile = sys.argv[2]
	ranks, seqLen = LoadRankFile(rankFile)

proteinName = os.path.basename(NOEfile).split('.')[0]

##print 'BEST', proteinName, ' '.join([ "%.3f" % x for x in best ])

NOEs = LoadNOEFile(NOEfile)
data = NOEs
data.sort(key=lambda x: x[1])

#print NOEs

if ranks is not None:
	data = []
	for d in NOEs:
		template = d[0]
		alignLen = ranks[template]
		rankScore = d[1] * (1+ 0.15*alignLen*1./seqLen)/1.15 
		data.append([template, rankScore, seqLen, alignLen] + d[1:])

	data.sort(key=lambda x: x[1] )
	


for d in data:
	dstr = [d[0] ] + [ "%.3f" % d0 for d0 in d[1:] ]
	print ' '.join(dstr)
