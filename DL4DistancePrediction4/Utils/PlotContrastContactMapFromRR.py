#!/usr/bin/env python

import numpy as np
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import sys
import os


if len(sys.argv) < 5:
    print 'Usage: PlotContrastContactMapFromRR.py protein_name RRFile1 RRFile2 NativeRRFile [MethodName1 | MethodName2]'
    print '	RRFile: a list of contacts. Each line has 2 to 5 fields and the first 2 fields are the indices of one contact.'
    exit(1)

target = sys.argv[1]
predFile1 = sys.argv[2]
predFile2 = sys.argv[3]
nativeFile = sys.argv[4]

methodName1 = 'RaptorX'
methodName2 = 'CCMpred'

if len(sys.argv) >= 7:
	methodName1 = sys.argv[5]
	methodName2 = sys.argv[6]

def LoadContacts(RRFile):
        fh = open(RRFile, 'r')
        content = list(fh)
        fh.close()

        data = []
        for line in content:
                fields=line.split()
		i = np.int32(fields[0])
		j = np.int32(fields[1])

		if i < j:
                	data.append( (i, j) )
		else:
                	data.append( (j, i) )

        return data


native = LoadContacts(nativeFile)
pred1 = LoadContacts(predFile1)
pred2 = LoadContacts(predFile2)

## determine the range of residues
allpairs = native + pred1 + pred2
start = min( [ min(p) for p in allpairs ] )
end = max( [ max(p) for p in allpairs ] )
L = end - start + 1

pred1_correct = np.array( list( set(pred1).intersection( set(native) ) ) ).astype(np.int32)
pred1_wrong = np.array( list( set(pred1) - set(native) ) ).astype(np.int32)

pred2_correct = np.array( list( set(pred2).intersection( set(native) ) ) ).astype(np.int32)
pred2_wrong = np.array( list( set(pred2) - set(native) ) ).astype(np.int32)

native2 = np.array( list(native) ).astype(np.int32)

plt.figure(1, figsize=(8, 8))
#plt.figure(1, figsize=(3.42, 3.42))
plt.axis([start, end, end, start])
plt.scatter(native2[:,0], native2[:,1], marker = 'o', color = 'grey', s =3)
plt.scatter(native2[:,1], native2[:,0], marker = 'o', color = 'grey', s =3)

## draw diagonal line
plt.plot((start, end), (start,end))

plt.scatter(pred1_correct[:,0], pred1_correct[:,1], marker = '*', color = 'r', s = 10)
plt.scatter(pred1_wrong[:,0], pred1_wrong[:,1], marker = 'x', color = 'g', s = 10)
#plt.annotate(methodName1, xy=(20 + start, end - 20 ), xycoords='axes points', horizontalalignment='left', verticalalignment='bottom', fontsize=20)
plt.annotate(methodName1, xy=(start + 0.05*L, end-0.05*L ), fontsize=15)

plt.scatter(pred2_correct[:,1], pred2_correct[:,0], marker = '*', color = 'r', s = 10)
plt.scatter(pred2_wrong[:,1], pred2_wrong[:,0], marker = 'x', color = 'g', s = 10)
#plt.annotate(methodName2, xy=(20 + start, end-20), xycoords='axes points', horizontalalignment='right', verticalalignment='upper', fontsize=20)
plt.annotate(methodName2, xy=(end-0.15*L, start+0.05*L), fontsize=15)

plt.title('Predicted and native contact matrix of ' + target, fontdict = {'fontsize' : 15})

#plt.savefig(target + ".ContactMap.pdf", format='pdf', dpi=300)
plt.savefig(target + ".ContactMap.tif", format='tif', dpi=300)

