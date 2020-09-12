import numpy as np
import sys
import theano
import theano.tensor as T
from numpy import random as rng

from NN4LogReg import NN4LogReg
from NN4Normal import NN4Normal
from NN4PhiPsi import NN4PhiPsi
from ResNet4Property import ResNet
import config

from config import Response2LabelName, Response2LabelType

class ResNet4Properties:

	def __init__(self, rng, seqInput, mask_seq=None, modelSpecs = None):
		"""
		seqInput has shape (batchSize, seqLen, n_in_seq)
		mask_seq has shape (batchSize, #cols_to_be_masked)
       		"""

		self.modelSpecs = modelSpecs	
		n_in_seq = modelSpecs['n_in_seq']
		n_hiddens_seq = modelSpecs['conv1d_hiddens']
		seq_repeats = modelSpecs['conv1d_repeats']
		n_hiddens_logreg = modelSpecs['logreg_hiddens']

		hwsz_seq=modelSpecs['halfWinSize_seq']
	
		self.mask_1d = mask_seq
		self.layers = []
        
		# sequence convolution 

		if modelSpecs['network'].startswith('ResNet'):
        		seqConv = ResNet(rng, input=seqInput, n_in=n_in_seq, n_hiddens=n_hiddens_seq, n_repeats=seq_repeats, halfWinSize=hwsz_seq, mask=mask_seq, activation=modelSpecs['activation'], batchNorm=modelSpecs['batchNorm'], version=modelSpecs['network'])
		else:
			print 'Unimplemented deep network type: ', modelSpecs['network']
			exit(-1)

		self.layers.append(seqConv)

		## conv_out has shape (batchSize, seqLen, seqConv.n_out)
        	conv_out = seqConv.output

		##flatten all
        	selected = conv_out.dimshuffle(2, 0, 1).flatten(2).dimshuffle(1, 0) 
		n_in4logreg = seqConv.n_out 

		self.outputList = []
		self.output_probList = []
		self.predictors = []

		self.params4var = []
		self.paramL14var = 0
		self.paramL24var = 0

		for res in modelSpecs['responses']:
			labelType = Response2LabelType(res)
			predictor = None

			if labelType.startswith('vonMise'):
				assert (config.responseValueDims[labelType] == 2)
            			predictor = NN4PhiPsi(rng=rng, input=selected, n_in=n_in4logreg, n_variables=config.responseValueDims[labelType], n_out=config.responseProbDims[labelType], n_hiddens=n_hiddens_logreg)
				self.params4var += predictor.params4var
				self.paramL14var += predictor.paramL14var
				self.paramL24var += predictor.paramL24var

			elif labelType.startswith('Gauss'):
            			predictor = NN4Normal(rng=rng, input=selected, n_in=n_in4logreg, n_variables=config.responseValueDims[labelType], n_out=config.responseProbDims[labelType], n_hiddens=n_hiddens_logreg)
				self.params4var += predictor.params4var
				self.paramL14var += predictor.paramL14var
				self.paramL24var += predictor.paramL24var

			elif labelType.startswith('Discrete'):
				assert (config.responseValueDims[labelType] == 1)
            			predictor = NN4LogReg(rng=rng, input=selected, n_in=n_in4logreg, n_out=config.responseProbDims[labelType], n_hiddens=n_hiddens_logreg)

			else:
				print 'incorrect response name or label type: ', res
				exit(-1)
				
	    		self.layers.append(predictor)
			self.predictors.append(predictor)

            		## output 
	    		y_pred = predictor.y_pred.reshape( (conv_out.shape[0], conv_out.shape[1], config.responseValueDims[labelType]) )
	    		output_prob = predictor.output.reshape( (conv_out.shape[0], conv_out.shape[1], config.responseProbDims[labelType]) )

	    		self.outputList.append( y_pred )
	    		self.output_probList.append( output_prob )

		## y_pred is the predicted target value
		## output_prob contains information for probability distribution of a target value
		self.y_pred = T.concatenate( self.outputList, axis=2 )
		self.output4prob = T.concatenate( self.output_probList, axis=2 ) 

        	self.params = []
		self.paramL1 = 0
		self.paramL2 = 0

		for layer in self.layers:
	    		self.params  += layer.params
	    		self.paramL1 += layer.paramL1
	    		self.paramL2 += layer.paramL2


    	## w has the same shape as z, it contains weight for each element in z
	## this function returns a vector
    	def loss(self, zList, useMeanOnly=False, weightList=None):
		losses = []
		if weightList is not None and len(weightList)>0:
			for predictor, z, w in zip(self.predictors, zList, weightList):
				zflat = z.dimshuffle(2, 0, 1).flatten(2).dimshuffle(1, 0)
				wflat = w.dimshuffle(2, 0, 1).flatten(2).dimshuffle(1, 0)
				losses.append(  predictor.loss(zflat, useMeanOnly, sampleWeight=wflat)   )
		else:
			for predictor, z in zip(self.predictors, zList):
				zflat = z.dimshuffle(2, 0, 1).flatten(2).dimshuffle(1, 0)
				losses.append(  predictor.loss(zflat, useMeanOnly )  )

		return T.stack(losses)


    	## the overall weighted errors, this function returns a vector
    	def errors(self, zList, weightList=None):
		errs = []
		if weightList is not None and len(weightList)>0:
			for predictor, z, w in zip(self.predictors, zList, weightList):
				zflat = z.dimshuffle(2, 0, 1).flatten(2).dimshuffle(1, 0)
				wflat = w.dimshuffle(2, 0, 1).flatten(2).dimshuffle(1, 0)
				e = predictor.errors(zflat, wflat)
				errs.append(e)
		else:
			for predictor, z in zip(self.predictors, zList):
				zflat = z.dimshuffle(2, 0, 1).flatten(2).dimshuffle(1, 0)
				e = predictor.errors(zflat )
				errs.append(e)

		return T.concatenate(errs)

	"""
    	## calculate the confusion matrix
    	def confusionMatrix(self, zList):

	    	confMs = []
		for z in zList:
			if z does not have type Discrete:
				## return an empty matrix
				m = T.ones( (1, 1) )
				confMs.append(m)
				continue

			selMatrix = T.ones_like(zList[0])
	    		if self.mask_1d is not None:
	        		mshape = self.mask_1d.shape
	        		selMatrix = T.set_subtensor(selMatrix[:, :mshape[1], :], self.mask_1d)
	        		selMatrix = T.set_subtensor(selMatrix[:, :, :mshape[1]], self.mask_1d.dimshuffle(0,2,1) )

            		dataLen = truth.shape[1]

	    		pred_truth = truth * numLabels + pred
	        	m = T.bincount( pred_truth, minlength=numLabels * numLabels ).reshape((numLabels, numLabels))	
			confMs.append( m )

	    	return T.stacklists( confMs )
	"""

	    

def BuildModel(modelSpecs, forTrain=True):
        rng = np.random.RandomState()

        ## x is for sequential features
        x = T.tensor3('x')

        ## mask for x 
        xmask = T.bmatrix('xmask')
        propertyPredictor = ResNet4Properties( rng, seqInput=x, mask_seq=xmask, modelSpecs=modelSpecs )

        ## labelList is a list of label matrices, each with shape (batchSize, seqLen, numLabels)
        labelList = []
        if forTrain:
                ## when this model is used for training. We need to define the label variable
		labelList = []
		for res in modelSpecs['responses']:
			labelType = Response2LabelType(res)
			if labelType.startswith('Discrete'):
                		labelList.append( T.itensor3('label4' + res ) )
			else:
                		labelList.append( T.tensor3('label4' + res ) )

        ## weightList is a list of label weight matices, each with shape (batchSize, seqLen, 1)
	## we always use weight to deal with residues without 3D coordinates
        weightList = []
        if len(labelList)>0:
                weightList = [ T.tensor3('weight4' + res ) for res in modelSpecs['responses'] ]

	if len(labelList)>0:
        	return propertyPredictor, x, xmask, labelList, weightList
	else:
        	return propertyPredictor, x, xmask

## need revision
def TestResNet4LocalProperties():

    x = T.tensor3('x')
    y = T.tensor4('y')
    xmask = T.bmatrix('xmask')
    ymask = T.btensor3('ymask')
    selection = T.wtensor3('selection')

    import cPickle
    fh = open('seqDataset4HF.pkl')
    data = cPickle.load(fh)
    fh.close()

    distancePredictor = ResNet4DistMatrix( rng=np.random.RandomState(), seqInput=x, matrixInput=y, n_in_seq=data[0][0].shape[2], n_in_matrix=data[1][0].shape[3], n_hiddens_seq=[3, 5], n_hiddens_matrix=[2], hwsz_seq=4, hwsz_matrix=4, mask_seq=xmask, mask_matrix=ymask)

    """
    f = theano.function([x, y, xmask, ymask], distancePredictor.output_1d)
    g = theano.function([x, y, xmask, ymask], distancePredictor.output_2d)
    """

    dataLen = 300
    batchSize = 60
    a = np.random.uniform(0, 1, (batchSize, dataLen, 20)).astype(np.float32)
    b = np.random.uniform(0, 1, (batchSize, dataLen, dataLen, 3)).astype(np.float32)
    amask = np.zeros((batchSize, 0)).astype(np.int8)
    bmask = np.zeros((batchSize, 0, dataLen)).astype(np.int8)
    sel = np.ones((batchSize, dataLen, dataLen)).astype(np.int8)
    #print a
    #print b
    c = np.random.uniform(0, 3, (batchSize, dataLen, dataLen)).round().astype(np.int8)
    np.putmask(c, c>=2, 2)
    """
    c[0, 1, 13]=1
    c[0, 2, 15]=1
    c[0, 4, 16]=1

    c[0, 1, 27]=1
    c[0, 2, 28]=1
    c[0, 4, 29]=1

    c[1, 0, 13]=2
    c[1, 1, 15]=2
    c[1, 3, 16]=2

    c[2, 0, 23]=2
    c[2, 1, 25]=2
    c[2, 3, 26]=2
    """
    #sel = c

    #out1d = f(a, b, amask, bmask)
    #out2d = g(a, b, amask, bmask)

    #print out1d
    #print out2d

    z = T.btensor3('z')
    loss = distancePredictor.loss(z, selection)
    errs = distancePredictor.ErrorsByRange(z)
    accs = distancePredictor.TopAccuracyByRange(z)
    confM = distancePredictor.confusionMatrix(z)

    h = theano.function([x, y, xmask, ymask, selection, z], confM, on_unused_input='ignore')
    #l, e, accu = h(a, b, amask, bmask, sel, c)

    cms = []
    for i in np.arange(5):
    	cm = h(data[0][i], data[1][i], data[2][i], data[3][i], data[4][i], data[5][i])
        print cm
	cms.append(cm)

    sumofcms = np.sum(cms, axis=0)*1.

    for i in xrange(sumofcms.shape[0]):
        sumofcms[i] /= np.sum(sumofcms[i])

    confusions = sumofcms
    print confusions
    print np.sum(confusions[0])
    print np.sum(confusions[1])
    print np.sum(confusions[2])

    """
    params = contactPredictor.params
    cost = loss + 0.1 * contactPredictor.paramL2
    grads = [ T.grad(cost, p) for p in params ]
    fg = theano.function([x, y, xmask, ymask, selection, z], grads)
    g0 = fg(a, b, amask, bmask, sel, c)

    param_shapes = [ p.get_value().shape for p in params ]
    print g0
    """

if __name__ == "__main__":

    TestResNet4LocalProperties()
    
