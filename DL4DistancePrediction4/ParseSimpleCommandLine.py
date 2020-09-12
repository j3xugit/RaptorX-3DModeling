import numpy as np
import cPickle
import os
import os.path
import sys
import random

import config
from config import Response2LabelName, Response2LabelType

import getopt

def Usage():
    	print 'python CalcEmpiricalRefState.py -y response -D trainData -s minibatchSize -k keyword_value'
	print '	-y: specify responses, e.g., CaCa+CbCb+CgCg:25C:2;Beta:2C:1. Several responses are separated by semicolon. Each response consists of response names, response label type and optionally weight factor separated by :. '
	print '    When several response names have the same label type and weight factor, you may separate them by + '
	print '    Response names can be atom pair types such as  CaCa, CbCb, NO, CgCg, CaCg and All where All includes all of them'
	print '    Response names can also be Beta (i.e., beta-pairing) and HB (hydrogen bonding)'
    	print '    2C for two labels; 3C for three labels; 12C for 12 labels; 25C for 25 labels. See config.py for more definition of distance labels'
	print ' '
    	print '	-s: the smallest and largest number of residue pairs in a minibatch, e.g. 90000,160000'
	print '	-D: a list of PKL files for the set of real protein data for train, validation and even prediction'
    	print '	-k: specify a set of keywords:value pair'

def ParseArguments(argv, modelSpecs):

    	try:
        	opts, args = getopt.getopt(argv,"y:D:s:k:",["response=","dataset=", "minibatchSize=", "kvpairs="])
        	print opts
    	except getopt.GetoptError:
       		Usage()
        	exit(1)

	## at least a data file shall be provided
    	if len(opts) < 3:
       		Usage()
        	exit(1)

    	for opt, arg in opts:
		if opt in ("-y", "--response"):
			modelSpecs['responseStr'] = arg
			responseSet = []
			weightSet = []

			## we examine the accuracy of the top seqLen * ratio predicted contacts where ratio is an element in ratioSet
			ratioSet = []

	    		fields = arg.split(';')
			for f in fields:
				## each f represents one response, consisting of response name, response label type and optionally weight factor
				words = f.split(':')
				if len(words) < 2:
					print 'Error: wrong format for the response argument: ', arg
                                        exit(1)

				## label name
				if words[0].upper() == 'All'.upper():
					names = config.allAtomPairTypes
				else:
					names = words[0].split('+')
		
				## label type
				labelType = words[1]

				if len(words) == 2:
					w = 1.
				else:
					w = np.float32(words[2])

				if labelType[0].isdigit():
					labelType = 'Discrete' + words[1]

				if labelType not in config.allLabelTypes:
					print labelType, 'is not a correct label type. It must be one of the following: ', config.allLabelTypes
                                        exit(1)
	
				for n in names:
					if n not in config.allLabelNames:
                                        	print n, 'is not an allowed response name. It must be one of the following: ', config.allLabelNames
                                        	exit(1)

					response = n + '_' + labelType
					responseSet.append(response)
					weightSet.append( w )
					##ratioSet.append( config.topRatios[n] )

			print responseSet, weightSet, ratioSet
			modelSpecs['responses'] = responseSet
			#modelSpecs['w4responses'] = np.array(weightSet).astype(theano.config.floatX)
			modelSpecs['w4responses'] = np.array(weightSet).astype(np.float32)
			##modelSpecs['topRatios'] = ratioSet


		elif opt in ("-D", "--dataset"):
			modelSpecs['dataset'] = [ f.strip() for f in arg.split(';') ]

		elif opt in ("-s", "--minibatchSize"):
	    		fields = arg.split(',')
            		minibatchSize = max(1000, int(fields[0]) )
            		modelSpecs['minibatchSize'] = minibatchSize
            		if len(fields) > 1:
                		modelSpecs['maxbatchSize'] = max(minibatchSize, int(fields[1]) )

		elif opt in ("-k", "--kvpairs"):
	     		items = [ i.strip() for i in arg.split(';') ]
	     		for item in items:
				k, v = item.split(':')
				modelSpecs[k] = v
				if v.upper() in ['True'.upper(), 'Yes'.upper()]:
		    			modelSpecs[k] = True
				if v.upper() in ['False'.upper(), 'No'.upper()]:
		    			modelSpecs[k] = False

				print 'Extra model specfication: ', k, modelSpecs[k]

        	else:
            		print Usage()
            		exit(1)


    	if modelSpecs['dataset'] is None:
       		print "Please provide one or more PKL files for the overall dataset fr training, validation and even prediction"
        	exit(1)

    	if modelSpecs['maxbatchSize'] < modelSpecs['minibatchSize']:
		print 'The largest number of data points in a batch is smaller than the smallest number. Please reset them.'
		exit(1)

	return modelSpecs

