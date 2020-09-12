import os
import numpy as np
import getopt

import config
from utilsNoT import IsNumber

def Usage():
    	print 'python TrainDistancePredictor.py -N ModelID -n network_type -y response -c conv1d_hiddens -d conv2d_hiddens -l dilation_factors -e logreg_hiddens -w halfWinSize -x method_for_seq2matrix -t trainMetaData_json -v validMetaData_json -p predfile_PKL -a trainAlgorithm -g regfactors -r restart_file -s minibatchSize -k keyword_value'
	print ' '
	print '-N: optional, specify a model ID (default None), e.g., A, B, 1, 2'
    	print '-n: specify a network type, e.g., ResNet2D or DilatedResNet2D (default). ' 
	print ' '
	print '-y: specify responses, e.g., CaCa+CbCb+CgCg:25C:2;Beta:2C:1. Several responses are separated by semicolon. Each response consists of label names, response label type and optionally weight factor separated by :. '
	print '    When several response names have the same label type and weight factor, you may separate them by + '
	print '    Response names can be atom pair types such as  CaCa, CbCb, NO, CgCg, CaCg and All where All includes all of them'
	print '    Response names can also be Beta (i.e., beta-pairing), HB (hydrogen bonding) and 5 different orientations (Ca1Cb1Cb2Ca2,N1Ca1Cb1Cb2,Cb1Cb2Ca2C2,C1Ca1Cb1Cb2,Cb1Cb2Ca2N2)'
    	print '    label types include 2C (binary classification), 3C (three labels), 12C, 25C, 34C and so on. See config.py for more definition of label types'
	print ' '
    	print '-c: the number of hidden units for 1d convolutional ResNet, e.g. 40,40:0,0 where the first part indicates the number of hidden neurons and the 2nd is the number of repeats'
    	print '-d: the number of hidden units for 2d convolutional ResNet, e.g. 30,30:1,2'
    	print '-e: the number of hidden units for the multi-layer logistic regression layers, e.g. 30,30'
    	print '-w: the half window size for both 1d and 2d convolutional layers, e.g. 10,5 for 1d and 2d, respectively'
	print '-l: the dilation factors, e.g., 2 or [1, 2, 2, 1]. It shall be a scalar or have the same length as conv2d_hiddens '
	print '-z: the half window size for Dilated network,  e.g., 7+2 or 7+[2, 1, 2, 1]. The first part (before +) is for 1d and the 2nd part for 2d, which shall have the same length as conv2d_hiddens '
	print ' '
    	print '-x: indicate how to convert 1d sequence to 2d matrix. Three strategies: SeqOnly for embedding amino acid (AA) pairs, '
	print '    Seq+SS for embedding the combination of AA and SS, OuterCat for outer concatenation of convoluted sequential features.' 
    	print '    You can use a combination of the above strategies, but SeqOnly and Seq+SS are exclusive to each other.' 
    	print '    Example: SeqOnly:4,6,12;OuterCat:80,35 or Seq+SS:4,6,12;OuterCat:80,35 where the numbers following one strategy are network parameters'
	print ' '
    	print '-a: training algorithm and learning rate (if needed), e.g., Adam:19+0.0002:2+0.00002 (default), SGDM:18+0.01:10+0.001'
	print '		In each stage (e.g., 19+0.0002), the first integer is the number of epochs to be run for this stage and the 2nd value is the learning rate '
	print ' ' 
    	print '-g: the L2 and L1 regularization factors, respectively, e.g., 0.0001 (default for L2) and 0 (default for L1)'
    	print '    where the first one is for L2 and the 2nd one (if available) for L1. By default the first one is 0.0001 and no L1 factor'
    	print '-s: the lower and upper bounds of a minibatch size in terms of the number of residue pairs, e.g. 90000,160000'
	print '	   The upper bound specifies the size of the max submatrix sampled from the original contact/distance matrix, e.g., 160000 implies the max size is 400*400'
	print '    Some minibatches may have fewer residue pairs than the lower bound, but all minibatches have at least half of the specified residue pairs, e.g., 90000 implies that a minibatch has at least 90000/2 residue pairs'
	print ' '
	#print '-D: a list of PKL files for the set of real protein data for train, validation and even prediction'
    	print '-t: a meta data file (json format) for training proteins. The real data is specified by the meta data'
	print '		The real data is saved as a PKL file for each protein'
    	print '-v: a meta data file (json format) for validation proteins. The real data is specified by the meta data'
    	print '-p: one or several PKL files for test proteins. The files have all the information needed for test proteins '
	print ' '
    	print '-r: specify a checkpoint file for restart'
    	print '-k: specify a set of keywords:value pair, e.g., UsePSICOV:no;UseTemplate:yes;TPLMemorySave:yes'

def ParseArguments(argv, modelSpecs):

    	try:
        	opts, args = getopt.getopt(argv,"N:n:y:t:v:p:c:d:e:w:x:a:r:g:s:k:l:z:",["modelID=", "network=","response=", "trainfile=","validfile=","predictfile=","conv1d=","conv2d=","logreg_hiddens=","halfWinSize=","seq2matrix=","algorithm=","restartfile=","regfactor=", "minibatchSize=","kvpairs=","dilation=","hwsz="])
        	print opts, args
    	except getopt.GetoptError:
       		Usage()
        	exit(1)

	## at least train and validation data files shall be provided
    	if len(opts) < 4:
       		Usage()
        	exit(1)

    	for opt, arg in opts:
		if opt in ("-N", "--modelID"):
			modelSpecs['ModelID'] = arg

       		elif opt in ("-n", "--network"):
			modelSpecs['network'] = arg
			if modelSpecs['network'] not in config.allNetworks:
				print 'Currently only support the network types in ', config.allNetworks
				exit(1)

		elif opt in ("-y", "--response"):
			modelSpecs['responseStr'] = arg
			responseSet = []
			weightSet = []

	    		fields = arg.split(';')
			for f in fields:
				## each f represents one or multiple responses, consisting of one or multiple response names, one response label type and optionally one weight factor, separated by :
				words = f.split(':')
				if len(words) < 2:
					print 'ERROR: this response has a wrong format : ', arg
                                        exit(1)

				## label name
				names = config.ParseLabelNames(words[0])
		
				## label type
				labelType = words[1]

				## weight factor for all responses encoded by names and labelType
				w = 1.
				if len(words) > 2:
					w = np.float32(words[2])
				assert w>0 

				if labelType[0].isdigit():
					labelType = 'Discrete' + words[1]

				if labelType not in config.allLabelTypes:
					print labelType, ' is not a correct label type. It must be one of the following: ', config.allLabelTypes
                                        exit(1)
	
				for n in names:
					#print n, labelType
					if n not in config.allLabelNames:
                                        	print n, ' is not an allowed response name. It must be one of the following: ', config.allLabelNames
                                        	exit(1)

					response = n + '_' + labelType
					responseSet.append(response)
					weightSet.append( w )

			print 'train a model for response: ', responseSet, ' with weight: ', weightSet
			modelSpecs['responses'] = responseSet
			modelSpecs['w4responses'] = np.array(weightSet).astype(np.float32)


        	elif opt in ("-a", "--algorithm"):
			modelSpecs['algStr'] = arg

	    		mean_var = arg.split(';')
			fields = mean_var[0].split(':')	

	    		if fields[0] not in config.allAlgorithms:
				print 'currently only the following algorithms are supported: ', config.allAlgorithms
	        		exit(1)
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
				modelSpecs['lrs4var'] = modelSpecs['lrs']
			

		elif opt in ("-x", "--seq2matrix"):

	    	## arg shall look like 'SeqOnly:3,4,5;OuterCat:80,35' or 'Seq+SS:3,4,5;OuterCat:80,35' where the numbers are parameters, e.g., for SeqOnly, 
		## the following numbers are the length of embedding vectors for short-, medium- and long-range residue pairs, respectively
	    	## if the vector length is set to 5 if not specified
	    	## if only a single length provided, then all ranges are set to this length
	    		options = {}
	    		fields = arg.split(';')
	    		for f in fields:
				columns = f.strip().split(':')
	        		if len(columns) !=2 :
		    			print 'unsupported format in the seq2matrix option: ', arg
		    			exit(1)

				if columns[0] not in config.allEmbeddingModes:
					print 'unsupported seq2matrix strategy: ', arg
                                        print 'allowed options: ', config.allEmbeddingModes
                                        exit(1)

	        		sizes = [ np.int32(x) for x in columns[1].split(',') ]
	        		options[ columns[0] ] = sizes
			
			##keep only one embedding method. Seq+SS always overrides SeqOnly
	    		if options.has_key('SeqOnly') and options.has_key('Seq+SS'):
				del options['SeqOnly']

	    		modelSpecs['seq2matrixMode'] = options

		#elif opt in ("-D", "--dataset"):
		#	modelSpecs['dataset'] = [ f.strip() for f in arg.split(';') ]
        	elif opt in ("-t", "--trainfile"):
	    		#modelSpecs['trainFile'] = [ f.strip() for f in arg.split(';') ]
	    		modelSpecs['trainFile'] = arg
        	elif opt in ("-v", "--validfile"):
	    		#modelSpecs['validFile'] = [ f.strip() for f in arg.split(';') ]
	    		modelSpecs['validFile'] = arg
        	elif opt in ("-p", "--predictfile"):
	    		#modelSpecs['predFile'] = [ f.strip() for f in arg.split(';') ]
	    		modelSpecs['predFile'] = arg

        	elif opt in ("-c", "--conv1d_hiddens"):
            		fields = arg.split(':')
	    		modelSpecs['conv1d_hiddens'] = map(int, fields[0].split(','))

	    		if len(fields) == 2:
	    			modelSpecs['conv1d_repeats'] = map(int, fields[1].split(','))

	    		assert len(modelSpecs['conv1d_hiddens']) == len(modelSpecs['conv1d_repeats'])

        	elif opt in ("-d", "--conv2d_hiddens"):
            		fields = arg.split(':')
	    		modelSpecs['conv2d_hiddens'] = map(int, fields[0].split(','))

	    		if len(fields) == 2:
	    			modelSpecs['conv2d_repeats'] = map(int, fields[1].split(','))

	    		assert len(modelSpecs['conv2d_hiddens']) == len(modelSpecs['conv2d_repeats']) 

		elif opt in ("-e", "--logreg_hiddens"):
	    		modelSpecs['logreg_hiddens']  = map(int, arg.split(','))

        	elif opt in ("-w", "--halfWinSize"):
            		halfWinSize = map(int, arg.split(','))
            		if len(halfWinSize) >= 1:
	        		modelSpecs['halfWinSize_seq'] = max(0, halfWinSize[0])

            		if len(halfWinSize) >= 2:
	        		modelSpecs['halfWinSize_matrix'] = max(0, halfWinSize[1])	

		elif opt in ("-l", "--dilation"):
			dilation = map(np.int32, arg.split(','))
			modelSpecs['conv2d_dilations'] = [ max(d, 1) for d in dilation ]

		elif opt in ("-z", "-hwsz"):
			fields = arg.split('+')
			modelSpecs['conv1d_hwsz'] = max(0, np.int32(fields[0]) )
			modelSpecs['conv2d_hwszs'] = [ max(0, w) for w in map(np.int32, fields[1].split(',') ) ]

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
            		minibatchSize = max(1000, int(fields[0]) )
            		modelSpecs['minibatchSize'] = minibatchSize
            		if len(fields) > 1:
                		modelSpecs['maxbatchSize'] = max(minibatchSize, int(fields[1]) )
			if len(fields) > 2:
				modelSpecs['boundingBoxOffset'] = max(0, int(fields[2]) )

		elif opt in ("-k", "--kvpairs"):
	     		items = [ i.strip() for i in arg.split(';') ]
	     		for item in items:
				k, v = item.split(':')
				if k not in modelSpecs.keys():
					print 'FATAL ERROR: not implemented options: ', k
					exit(1)

				modelSpecs[k] = v
				if IsNumber(v):
					modelSpecs[k] = np.float32(v)

				if v.upper() in ['True'.upper(), 'Yes'.upper()]:
		    			modelSpecs[k] = True

				if v.upper() in ['False'.upper(), 'No'.upper()]:
		    			modelSpecs[k] = False

				if k.lower() == 'activation':
		    			if v.upper() in config.allActivations:
						modelSpecs['activation'] = v.upper()
		    			else:
						print 'allowed activation functions: ', config.allActivations
						exit(1)

				print 'Extra model specfication: ', k, modelSpecs[k]

        	else:
            		Usage()
            		exit(1)

	## check if conv2d_hwszs and conv2d_dilations have the same length as conv2d_repeats
	if modelSpecs['network'].startswith('DilatedResNet'):
		if len(modelSpecs['conv2d_dilations']) == 1:
			modelSpecs['conv2d_dilations'] = modelSpecs['conv2d_dilations'] * len(modelSpecs['conv2d_repeats'])
		elif len(modelSpecs['conv2d_dilations']) != len(modelSpecs['conv2d_repeats']):
			print 'ERROR: conv2d_dilations and conv2d_repeats do not have the same length.'
			exit(1)

		if len(modelSpecs['conv2d_hwszs']) == 1:
			modelSpecs['conv2d_hwszs'] = modelSpecs['conv2d_hwszs'] * len(modelSpecs['conv2d_repeats'])
		elif len(modelSpecs['conv2d_hwszs']) != len(modelSpecs['conv2d_repeats']):
			print 'ERROR: conv2d_hwszs and conv2d_repeats do not have the same length.'
			exit(1)

    	## if ExSS is set to True (i.e., SS is not used), then revise the embedding mode
    	if  modelSpecs['seq2matrixMode'].has_key('Seq+SS') and modelSpecs.has_key('ExSS') and (modelSpecs['ExSS'] is True):
		modelSpecs['seq2matrixMode']['SeqOnly'] = modelSpecs['seq2matrixMode']['Seq+SS']
		modelSpecs['seq2matrixMode'].pop('Seq+SS')

	"""
    	if modelSpecs['dataset'] is None:
       		print "Please provide one or more PKL files for the overall dataset for training, validation and even prediction"
        	exit(1)
	"""

    	if modelSpecs['trainFile'] is None or not os.path.isfile(modelSpecs['trainFile']):
       		print "Please provide one file ending with .json for training meta data"
        	exit(1)

    	if modelSpecs['validFile'] is None or not os.path.isfile(modelSpecs['validFile']):
       		print "Please provide one file ending with .json for validation meta data"
        	exit(1)


    	if modelSpecs['maxbatchSize'] < modelSpecs['minibatchSize']:
		print 'The largest number of data points in a batch is smaller than the smallest number. Please reset them.'
		exit(1)

	ratioSet = []
	for res in modelSpecs['responses']:
		ratioSet.append( config.topRatios[config.Response2LabelName(res) ] )
	modelSpecs['topRatios'] = ratioSet

	return modelSpecs

