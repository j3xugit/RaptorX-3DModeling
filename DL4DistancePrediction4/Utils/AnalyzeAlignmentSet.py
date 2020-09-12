import numpy as np
import sys

AlignSetFile = sys.argv[1]
SubSetFile = sys.argv[2]

f=open(SubSetFile, 'r')
proteins = [ line.strip() for line in list(f) ]
f.close()

f=open(AlignSetFile, 'r')
aligns = [ line.split() for line in list(f) ]
f.close()

## calculate the average SeqID and logP of the whole alignment set

seqIDs = [ np.float32(a[2]) for a in aligns ]
Ps = [ np.log(np.float32(a[3])) for a in aligns if a[3]!='0.0' ]

print np.average(seqIDs), np.average(Ps)

## calculate the average SeqID and logP of the subset

proteins = set(proteins)
subaligns = [ a for a in aligns if a[1] in proteins ]

seqIDs = [ np.float32(a[2]) for a in subaligns ]
Ps = [ np.log(np.float32(a[3])) for a in subaligns if a[3]!='0.0' ]

print np.average(seqIDs), np.average(Ps)
