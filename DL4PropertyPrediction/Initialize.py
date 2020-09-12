import theano.tensor as T

def InitializeModelSpecs():
        modelSpecs = dict()
        modelSpecs['trainFile'] = None
        modelSpecs['validFile'] = None
        modelSpecs['predFile'] = None
        modelSpecs['checkpointFile'] = None

        modelSpecs['network'] = 'ResNet1D'
        modelSpecs['responseStr'] = 'PhiPsi_vonMise'
        modelSpecs['responses'] = ['PhiPsi_vonMise']
        modelSpecs['w4responses'] = [1.]

        modelSpecs['algorithm'] = 'Adam'
        modelSpecs['numEpochs'] = [19, 2]
        modelSpecs['lrs'] = [np.float32(0.0002),  np.float32(0.0002)/10 ]

        modelSpecs['algorithm4var'] = 'Adam'
        modelSpecs['numEpochs4var'] = modelSpecs['numEpochs']
        modelSpecs['lrs4var'] = modelSpecs['lrs']

        modelSpecs['algStr'] = 'Adam:21:0.00022'

        modelSpecs['minibatchSize'] = 200
        modelSpecs['maxbatchSize'] = 1500

        modelSpecs['validation_frequency'] = 100
        modelSpecs['patience'] = 5

        ##default number of hidden units at 1d convolutional layer
        modelSpecs['conv1d_hiddens'] = [80, 100, 120, 140]
        modelSpecs['conv1d_repeats'] = [ 0,  0,  0,  0]

        ## for the logistic regression at the final stage
        modelSpecs['logreg_hiddens'] = [ ]

        modelSpecs['halfWinSize_seq'] = 5
        modelSpecs['activation'] = T.nnet.relu

        modelSpecs['L2reg'] = np.float32(0.0001)
        modelSpecs['w4diso'] = 3.

        modelSpecs['batchNorm'] = True
        modelSpecs['UseSampleWeight'] = True
        modelSpecs['UsePSFM'] = False
        modelSpecs['UseTemplate'] = False
        modelSpecs['UseSequenceEmbedding'] = False
        modelSpecs['UseOneHotEncoding'] = False

        return modelSpecs

