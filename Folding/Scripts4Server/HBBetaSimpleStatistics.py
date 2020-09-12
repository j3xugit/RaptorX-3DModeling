import numpy as np
import os
import sys

rawfile = sys.argv[1]

f= open(rawfile, 'r')
content = [ line.strip() for line in list(f) ]
f.close()

threshold4Beta = 10

if len(sys.argv) > 2:
	threshold4Beta = int(sys.argv[2])

numIgnoredChains = 0
keptChains = []
for line in content:
	fields = line.split()
	numBetaPairings = int(fields[2])
	if numBetaPairings < threshold4Beta:
		numIgnoredChains += 1
		continue

	keptChains.append( map(int, fields[1:]) )

print '#ignored chains: ', numIgnoredChains

keptChains = np.array( keptChains)
summary = np.sum(keptChains, axis=0)
print summary
print summary[2:6]*1.0/summary[-4:]
print summary[7:11]*1.0/summary[-4:]
	
