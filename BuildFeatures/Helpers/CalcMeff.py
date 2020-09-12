import numpy as np
import sys
import os

if len(sys.argv) <2:
	print 'python CalcMeff.py a3m_file or a2m_file [seqid threshold in [0, 1], default 0.7]'
	exit(-1)

threshold = 0.7
a3mFile = sys.argv[1]

if len(sys.argv) >= 3:
	threshold = np.float32(sys.argv[2])

fh = open(a3mFile, 'r')
content = [ line.strip() for line in list(fh) ]
fh.close()

import string

content2 = [ line.translate(None, string.ascii_lowercase) for line in content if not line.startswith('>') ]

print '# loaded sequences in the a3m file: ', len(content2)
if len(content2)<1:
	print 'there is no sequence in this MSA file '
	exit(-1)

queryLen = len(content2[0])

##calculate effective lengths
effLens = []

for line in content2:
	effLen = np.sum( [ c!='-' for c in line[ : min(queryLen, len(line)) ] ] )
	effLens.append(effLen)

requiredIDs = [ np.int32( effLen * threshold +0.5 ) for effLen in effLens ]

##calculate similarity

sim = np.ones( (len(content2)), dtype=np.int32)

for line1, i in zip(content2, range(len(content2))):
	if i%100 == 0 :
		print '#processing %d-th sequence...' %i
	len1 = min(len(line1), queryLen)

	for line2, j in zip(content2[i+1:], range(i+1, len(content2))):
		len2 = min(len(line2), queryLen)
		minLen = min(len1, len2)

		ID = np.sum ( [ (c1 != '-' and (c1 !='.') and (c1==c2) ) for c1, c2 in zip(line1[: minLen], line2[: minLen]) ] )
		sim[i] += (ID >= requiredIDs[i] )
		sim[j] += (ID >= requiredIDs[j] )

Meff = np.sum( [ 1./s for s in sim ] )

print 'Meff=', Meff				
