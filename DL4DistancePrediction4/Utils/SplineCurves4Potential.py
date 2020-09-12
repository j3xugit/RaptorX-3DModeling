import os
import sys
from scipy.interpolate import UnivariateSpline
import numpy as np

##this script generates a matrix of spline funcs and their derivatives
##potMatrix is the potential matrix and distCutoffs is the distance boundaries for discretization.
##Currently, distCutoffs = [0, 4, 4.5, ...] or [ 0, 4.5, 5.0, ...]. Otherwise the results may be incorrect
##the last entry in potMatrix[i, j] shall be 0

## return two dicts: spline and dsplines. Can access them by [i][j] where i < j
def GenSplineFunc(potMatrix, distCutoffs, deriv=True):

	assert distCutoffs[0] == 0.
	assert potMatrix[:,:,0].all() == 0.0

	x = distCutoffs[1:]
        assert x[0] == 4.0 or x[0] == 4.5

        ##here we assume that binWidth = 0.5
        binWidth = 0.5
	for i in xrange(1, len(distCutoffs)-1):
		assert binWidth == (distCutoffs[i+1] - distCutoffs[i])

        xmin = 0.
        xPrefix = np.arange(xmin, x[0], binWidth).tolist()

	## yPrefix[-1] corresponds to distance=3.0, i.e., at this point, potential = 5.
        yPrefix = np.arange(35.0, 0., -5.).tolist()

        xmax = x[-1] + binWidth
	##the potential for distance >= xmax is assumed to be 0
        xSuffix = [ xmax ]
        #ySuffix = [ 0.0 ]

	## the data points for spline		
        xk = sum([xPrefix, x.tolist(), xSuffix], [])
	xk = np.array(xk) - binWidth/2.
        #print xk

	splines = dict()
	if deriv:
		dsplines = dict()

        size = potMatrix.shape
        for i in xrange(size[0]):
		splines[i] = dict()
		if deriv:
			dsplines[i] = dict()

                for j in xrange(i+1, size[1]):
                	y = potMatrix[i, j]

                        ## set potential for distance=3.5-binWidth/2. y[0] is the potential for distance = 4.0 or 4.5 minus binWidth/2
                        E435 = (5. + y[0])/2
			"""
                        if y[0] < 0:
                        	E435 = 0.
			"""

                        yPrefix2 = [E435 ] + [ y[0] ] * (len(xPrefix) - len(yPrefix) -1)
                        #yk = yPrefix + yPrefix2 + y[:-1].tolist() + ySuffix
                        yk = yPrefix + yPrefix2 + y.tolist()

                        #print yk
                        assert len(xk) == len(yk)

                        spl = UnivariateSpline(xk, yk, s=0)
			splines[i][j] = spl

			if deriv:
				spld = spl.derivative()
				dsplines[i][j] = spld

	if deriv:
		return splines, dsplines

	return splines

