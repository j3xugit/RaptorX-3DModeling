import numpy as np
import cPickle
import theano.tensor as T
import theano
import os
import os.path
import sys
import time
import datetime
import random
import gzip

import config
from config import Response2LabelType, Response2LabelName
import getopt

def Usage():
    	print 'python TrainPropertyPredictor.py -n network_type -y response -c conv1d_hiddens -e logreg_hiddens -w halfWinSize -t trainfile -v validfile -p predfile -a trainAlgorithm -g regfactors -r restart_file -i numEpochs -s minibatchSize -k keyword_value'

    	print '-n: specify a network type, e.g., SeqResNet1D (default). We may use 1d CNN or LSTM or their combination for 1d sequence.' 

	print '-y: specify response type (e.g., PhiPsi_vonMise2d4 or SS3_Discrete3C or SS8_Discrete8C:1;PhiPsi_vonMise2d2:2;pACC_Gauss1d:1) where SS8 is property name, Discrete8C is property value type and 1 is the weight factor for this property'
	print '	   Gauss1d indicates to model the reponse variables using single-variable Normal distribution'
	print '	   Discrete3C indicates using a finite number of classes, here 3C indicates 3 classes'

    	print '-c: the number of hidden units for 1d convolutional ResNet, e.g. 80,200:2,1 where the first part indicates the number of hidden neurons and the 2nd is the number of repeats'
    	print '-e: the number of hidden units for the multi-layer logistic regression layers (default empty), e.g. 30,30'
    	print '-w: the half window size used by convolutional layers (default 3), e.g. 5 '

	print '-a: training algorithm and learning rate (if needed), e.g., Adam:19+0.0002:2+0.00002 (default), SGDM:18+0.01:10+0.001'
        print '         In each stage (e.g., 19+0.0002), the first integer is the number of epochs to be run for this stage and the 2nd value is the learning rate '

    	print '-g: the L2 and L1 regularization factors, respectively, e.g., 0.0001 (default for L2) and 0 (default for L1)'
    	print '    where the first one is for L2 and the 2nd one (if available) for L1. By default the first one is 0.0001 and no L1 factor'
    	print '-i: the total number of epochs used by the optimization algorithm, default 20'
    	print '-s: the smallest and largest number of residues in a minibatch, e.g. 600,1500'

    	print '-t: specify one or more training files'
    	print '-v: specify one or more validation files'
    	print '-p: specify one or more files containing data to be predicted'

    	print '-r: specify a checkpoint file for restart'
    	print '-k: specify a set of keywords:value pair, e.g., UseTemplate:yes;UseSequenceEmbedding:yes'

def ParseArguments(argv, modelSpecs):

    	try:
        	opts, args = getopt.getopt(argv,"n:y:t:v:p:c:e:w:a:r:g:i:s:k:",["network=","response=","trainfile=","validfile=","predictfile=","conv1d=","logreg_hiddens=","halfWinSize=","algorithm=","restartfile=","regfactor=", "numEpochs=", "minibatchSize=","kvpairs="])
        	print opts, args
    	except getopt.GetoptError:
       		Usage()
        	exit(-1)

	##we need at least a training file and a validation file
    	if len(opts) < 2:
       		Usage()
        	exit(-1)

    	for opt, arg in opts:
       		if opt in ("-n", "--network"):

			## for training, we set ResNet1D to ResNet1DV21, which is almost same as ResNet1D
			if arg == 'ResNet1D':
				modelSpecs['network'] = 'ResNet1DV21'
			else:
				modelSpecs['network'] = arg

			if modelSpecs['network'] not in config.allNetworks:
				print 'Currently only support the network types in ', config.allNetworks
				exit(-1)

		elif opt in ("-y", "--response"):
			modelSpecs['responseStr'] = arg
			responseSet = []
			weightSet = []

			fields = arg.split(';')
			for f in fields:
				words = f.split(':')
				if len(words) < 1:
					print 'Error: wrong format for the response argument: ', arg
					exit(-1)

				responseSet.append(words[0])
				if len(words) >= 2:
					weightSet.append(max(0, np.float32(words[1]) ) )
				else:
					weightSet.append(1.)

			modelSpecs['responses'] = responseSet
			modelSpecs['w4responses'] = np.array(weightSet).astype(theano.config.floatX)
			print 'responses=', modelSpecs['responses'], ' their weights= ', modelSpecs['w4responses']

			## correctness check
			for res in modelSpecs['responses']:
	    			fields = res.split('_')
				assert ( len(fields) == 2 )

				if fields[0] not in config.allLabelNames:
                    			print 'Please specify a correct property name. It must be one of the following: ', config.allLabelNames
                    			exit(-1)
				if fields[1] not in config.allLabelTypes:
                    			print 'Please specify a correct property value type. It must be one of the following: ', config.allLabelTypes
                    			exit(-1)


        	elif opt in ("-a", "--algorithm"):
			modelSpecs['algStr'] = arg
			mean_var = arg.split(';')
                        fields = mean_var[0].split(':')

                        if fields[0] not in config.allAlgorithms:
                                print 'currently only the following algorithms are supported: ', config.allAlgorithms
                                exit(-1)
                        modelSpecs['algorithm'] = fields[0]

                        if len(fields) > 1:
                                numEpochs = []
                                lrs = []
                                for index, f in zip( xrange(len(fields)-1 ), fields[1: ]):
                                        f2 = f.split('+')
                                        assert ( len(f2) == 2)
                                        numEpochs.append(np.int32(f2[0]) )
                                        lrs.append(np.float32(f2[1]) )
                                modelSpecs['numEpochs'] = numEpochs
                                modelSpecs['lrs' ] = lrs

                        if len(mean_var) > 1:
                                fields = mean_var[1].split(':')
                        	if fields[0] not in config.allAlgorithms:
                                	print 'currently only the following algorithms are supported: ', config.allAlgorithms
                                	exit(-1)
                                modelSpecs['algorithm4var'] = fields[0]
                                if len(fields) > 1:
                                        numEpochs = []
                                        lrs = []
                                        for index, f in zip( xrange(len(fields)-1 ), fields[1: ]):
                                                f2 = f.split('+')
                                                assert ( len(f2) == 2)
                                                numEpochs.append(np.int32(f2[0]) )
                                                lrs.append(np.float32(f2[1]) )
                                        modelSpecs['numEpochs4var'] = numEpochs
                                        modelSpecs['lrs4var' ] = lrs
			else:
                                modelSpecs['algorithm4var'] = modelSpecs['algorithm']
                                modelSpecs['numEpochs4var'] = modelSpecs['numEpochs']
                                modelSpecs['lrs4var' ] = modelSpecs['lrs']
				


        	elif opt in ("-t", "--trainfile"):
	    		modelSpecs['trainFile'] = [ f.strip() for f in arg.split(';') ]
        	elif opt in ("-v", "--validfile"):
	    		modelSpecs['validFile'] = [ f.strip() for f in arg.split(';') ]
        	elif opt in ("-p", "--predictfile"):
	    		modelSpecs['predFile'] = [ f.strip() for f in arg.split(';') ]

        	elif opt in ("-c", "--conv1d_hiddens"):
            		fields = arg.split(':')
	    		modelSpecs['conv1d_hiddens'] = map(int, fields[0].split(','))

	    		if len(fields) == 2:
	    			modelSpecs['conv1d_repeats'] = map(int, fields[1].split(','))
			else:
				modelSpecs['conv1d_repeats'] = [0] * len(modelSpecs['conv1d_hiddens'])

			assert  len(modelSpecs['conv1d_repeats']) == len(modelSpecs['conv1d_hiddens'] )

		elif opt in ("-e", "--logreg_hiddens"):
	    		modelSpecs['logreg_hiddens']  = map(int, arg.split(','))

        	elif opt in ("-w", "--halfWinSize"):
            		halfWinSize = np.int32(arg)
	        	modelSpecs['halfWinSize_seq'] = max(0, halfWinSize)

        	elif opt in ("-g", "--regfactor"):
            		regs = map(np.float32, arg.split(','))
            		if len(regs)>0 and regs[0]>0:
				modelSpecs['L2reg'] = regs[0]

            		if len(regs)>1 and regs[1]>0:
	        		modelSpecs['L1reg'] = regs[1]

        	elif opt in ("-r", "--restart"):
			modelSpecs['checkpointFile'] = arg

		elif opt in ("-s", "--minibatchSize"):
	    		fields = arg.split(',')
            		minibatchSize = max(10, np.int32(fields[0]) )
            		modelSpecs['minibatchSize'] = minibatchSize
            		if len(fields) > 1:
                		modelSpecs['maxbatchSize'] = max(minibatchSize, np.int32(fields[1]) )

		elif opt in ("-k", "--kvpairs"):
	     		items = [ i.strip() for i in arg.split(';') ]
	     		for item in items:
				k, v = item.split(':')
				if k not in config.allKeywords:
					print 'Error: the input keyword ', k, ' not in the allowed list: ', config.allKeywords
					exit(-1)

				if k.lower() == 'activation':
		    			if v.upper() == 'RELU':
						modelSpecs['activation'] = T.nnet.relu
		    			elif v.upper() == 'TANH':
						modelSpecs['activation'] = T.tanh
		    			else:
						print 'unsupported activation function'
						sys.exit(-1)
				elif k.lower() == 'w4disorder':
					modelSpecs['w4diso'] = np.float32(v)

				elif v.upper() in ['True'.upper(), 'Yes'.upper()]:
		    			modelSpecs[k] = True
				elif v.upper() in ['False'.upper(), 'No'.upper()]:
		    			modelSpecs[k] = False
				else:
					modelSpecs[k] = v

        	else:
            		print Usage()
            		exit(-1)

    	if modelSpecs['trainFile'] is None:
       		print "Please provide one or more training files ending with .pkl and separated by ;"
        	exit(-1)

    	if modelSpecs['validFile'] is None:
       		print "Please provide one or more validation files ending with .pkl and separated by ;"
        	exit(-1)

    	if modelSpecs['maxbatchSize'] < modelSpecs['minibatchSize']:
		print 'The largest number of data points in a batch is smaller than the smallest number. Please reset them.'
		exit(-1)

	return modelSpecs

