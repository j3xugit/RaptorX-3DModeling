import numpy as np
import sys
import os
import cPickle
import random

numTraces = 100000
numSteps = 600
CaCaLength = 3.8
CaCbLength = 1.53

nDims = 3

eps = sys.float_info.epsilon

def str_display(ls):
        if not isinstance(ls, (list, tuple, np.ndarray)):
                str_ls = '{0:.7f}'.format(ls)
                return str_ls

        str_ls = ['{0:.7f}'.format(v) for v in ls ]
        str_ls2 = ' '.join(str_ls)
        return str_ls2


##sample on a sphere, s=size and r=radius of sphere
def SampleOnSphere(s, r):
	#steps = np.random.normal(size=(numTraces, numSteps, nDims))
	steps = np.random.normal(size=s).astype(np.float32)
	#print steps

	## convert each step to a point on a sphere with radius=r
	steps_sqr = np.square(steps)
	steps_norm = np.sqrt( np.sum(steps_sqr, axis=-1, keepdims=True) + eps )
	steps_sphere_unit = np.divide(steps, steps_norm)
	steps_sphere = r * steps_sphere_unit

	return steps_sphere

## calculate the length of an array of vectors
def CalcNorm(traces):
	traces_sqr = np.square(traces)
	norms = np.sqrt( np.sum(traces_sqr, axis=-1) )
	return norms

## calculate the occurring frequencey of an array of vectors by its length discretized according to bins
def BinCount(traces, bins):
	norms = CalcNorm(traces)
	normsT = norms.transpose()
	labelMatrix = np.digitize(normsT, bins) - 1

	numLabels = len(bins)

	freqs = []
	for row in labelMatrix:
		count = np.bincount(row, minlength=numLabels )
		freq = count * 1. / np.sum(count)
		freqs.append(freq)
	
	return np.array(freqs)

## here bins is a list of 68 points used to discretize distance
bins=np.array ( np.linspace(0, 34.0, num=69).tolist()  ).astype(np.float32)

## initialize random seed
a = list( str(os.getpid())  +  os.urandom(8) )
random.shuffle(a)
seed = ''.join(a)
print 'random seed: ', seed
random.seed(a=seed)


freqs = []
for i in xrange(100):

	## in this RW, there are numSteps + 1 residues, each trajectory has numSteps Ca-Ca virtual bonds
	## random walk for Ca trace
	CaCaSteps = SampleOnSphere(s=(numTraces, numSteps, nDims), r=CaCaLength).astype(np.float32)
	##assume all traces start from the original point of Ca
	CaCaTraces = np.cumsum(CaCaSteps, axis=1)

	## sample Ca-Cb vectors for each residue. Each residue has one Cb atom and thus, each trajectory has numSteps + 1 Cb atoms
	CaCbSteps = SampleOnSphere(s=(numTraces, numSteps+1, nDims), r=CaCbLength).astype(np.float32)

	## add Cb atoms for residues 1, 2, ..., numSteps
	CaCbTraces = CaCaTraces + CaCbSteps[:, 1:, :]

	## add Cb atoms for residue 0
	CbCbTraces = CaCbTraces + CaCbSteps[:, 0:1, :]

	#freq = BinCount(CaCaTraces, bins)
	freq = BinCount(CbCbTraces, bins)
	freqs.append(freq)

final_freq = np.mean(freqs, axis=0)

for row in final_freq:
	print str_display(row)

##np.set_printoptions(threshold=sys.maxsize)
##print freqs

savefile = 'RWRefState-' + str(os.getpid()) + '.pkl'
fh = open(savefile, 'wb')
cPickle.dump(final_freq, fh, protocol=cPickle.HIGHEST_PROTOCOL)
fh.close()
