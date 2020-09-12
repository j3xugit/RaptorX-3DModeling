import numpy as np
import os
import sys
import cPickle

from config import Response2LabelName, Response2LabelType

SS3Code2Letter = ['H', 'E', 'L' ]
SS3Letter2Code = { 'H':0, 'E':1, 'L':2, 'C':2}

SS8Code2Letter = ['H', 'G', 'I', 'E', 'B', 'T', 'S', 'L' ]
SS8Letter2Code = {'H':0, 'G':1, 'I':2, 'E':3, 'B':4, 'T':5, 'S':6, 'L':7, 'C':7 }

##convert 8-state to 3-state. 
##different conversion rules may have different prediction accuracy
SS8Letter2SS3Code = {'H':0, 'G':0, 'I':0, 'E':1, 'B':1, 'T':2, 'S':2, 'L':2, 'C':2 }
SS8Letter2SS3Letter = {'H':'H', 'G':'H', 'I':'H', 'E':'E', 'B':'E', 'T':'L', 'S':'L', 'L':'L', 'C':'L' }
SS8Code2SS3Code = [0, 0, 0, 1, 1, 2, 2, 2]

##vector coding of SS
SS8Coding=dict()
SS8Coding['H']=np.array([1, 0, 0, 0, 0, 0, 0, 0 ]).astype(np.int32)
SS8Coding['G']=np.array([0, 1, 0, 0, 0, 0, 0, 0 ]).astype(np.int32)
SS8Coding['I']=np.array([0, 0, 1, 0, 0, 0, 0, 0 ]).astype(np.int32)
SS8Coding['E']=np.array([0, 0, 0, 1, 0, 0, 0, 0 ]).astype(np.int32)
SS8Coding['B']=np.array([0, 0, 0, 0, 1, 0, 0, 0 ]).astype(np.int32)
SS8Coding['S']=np.array([0, 0, 0, 0, 0, 1, 0, 0 ]).astype(np.int32)
SS8Coding['T']=np.array([0, 0, 0, 0, 0, 0, 1, 0 ]).astype(np.int32)
SS8Coding['L']=np.array([0, 0, 0, 0, 0, 0, 0, 1 ]).astype(np.int32)
SS8Coding['C']=np.array([0, 0, 0, 0, 0, 0, 0, 1 ]).astype(np.int32)

SS3Coding=dict()
SS3Coding['H']=np.array([1, 0, 0]).astype(np.int32)
SS3Coding['E']=np.array([0, 1, 0]).astype(np.int32)
SS3Coding['L']=np.array([0, 0, 1]).astype(np.int32)
SS3Coding['C']=np.array([0, 0, 1]).astype(np.int32)

SS3Coding['G'] = SS3Coding['H']
SS3Coding['I'] = SS3Coding['H']
SS3Coding['B'] = SS3Coding['E']
SS3Coding['S'] = SS3Coding['L']
SS3Coding['T'] = SS3Coding['L']


## the pACC cutoff for B and M are 10% and 42%, respectively
ACCCode2Letter = ['B', 'M', 'E']
ACCLetter2Code = {'B':0, 'M':1, 'E':2 }

## the CLE coding
CLEsize = 18 ## CLE has letters from A to R inclusive
CLECode2Letter = [ chr(c) for c in (np.array( range(CLEsize) ) + ord('A') ) ]
CLELetter2Code = dict()
for letter in CLECode2Letter:
	CLELetter2Code[letter] = ord(letter) - ord('A')

CLECoding = dict()
CLEeye = np.eye(CLEsize, M=CLEsize, dtype=np.int32)
for letter in CLECode2Letter:
	CLECoding[letter] = CLEeye[ord(letter) - ord('A'), : ]
	
##convert the coding representation (0, 1, 2) of a property (e.g., secondary structure) to meaningful representation (e.g., H, E, L)
##coding is an integer vector or a 2d matrix with shape (seqLen, 1)
def Coding2String(coding, response):

	if coding.ndim == 2:
		code = coding[:,0]
	else:
		code = coding

	labelType = Response2LabelType(response)

	if response.startswith('SS'):
		if labelType.endswith('3C'):
			str = ''.join( [ SS3Code2Letter[c] for c in code ] ) 
		elif labelType.endswith('8C'):
			str = ''.join( [ SS8Code2Letter[c] for c in code ] ) 
		else:
			print 'ERROR: unsupported response and labelType: ', response
			exit(1)
		return str

	if response.startswith('ACC'):
		assert ( labelType.endswith('3C') )
		str = ''.join( [ ACCCode2Letter[c] for c in code ] )
		return str

	if response.startswith('CLE'):
		assert ( labelType.endswith('18C') )
		str = ''.join( [ CLECode2Letter[c] for c in code ] )
		return str
	
	print 'ERROR: unsupported response: ', response
	exit(1)

##convert the meaningful representation (H, E, L) of a property (e.g., secondary structure) to coding representation (0, 1, 2)
## return a matrix with shape (seqLen, 1)
def String2Coding(str, response):

	labelType = Response2LabelType(response)

	if response.startswith('SS'):
		if labelType.endswith('3C'):
			code = [ SS3Letter2Code[c] for c in str ]  
			code = np.array(code).astype(np.int32).reshape( (len(str), 1) )
		elif labelType.endswith('8C'):
			code = [ SS8Letter2Code[c] for c in str ] 
			code = np.array(code).astype(np.int32).reshape( (len(str), 1) )
		else:
			print 'ERROR: unsupported response and labelType: ', response
			exit(1)
		return code

	if response.startswith('ACC'):
		assert ( labelType.endswith('3C') )
		code = [ ACCLetter2Code[c] for c in str ] 
		code = np.array(code).astype(np.int32).reshape( (len(str), 1) )
		return code

	if response.startswith('CLE'):
		assert ( labelType.endswith('18C') )
		code = [ CLELetter2Code[c] for c in str ]
		code = np.array(code).astype(np.int32).reshape( (len(str), 1) )
		return code

	print 'ERROR: unsupported response: ', response
	exit(1)

def LoadPredictedProperties(propertyFile):
	f = open(propertyFile, 'rb')
	data = cPickle.load(f)
	f.close()
	return data

##evaluate the prediction accuracy for a single protein. prediction is a dictionary with response as the key
def EvaluateSinglePropertyPrediction(prediction, nativeLabelFile):

	from DataProcessor import LoadNativeLabelsFromFile

	errors = dict()
	nativeLabels = LoadNativeLabelsFromFile(nativeLabelFile)

	for response, pred in prediction.iteritems():
		native = nativeLabels[Response2LabelName(response)]
		missing = nativeLabels['Missing']

		if response.startswith('DISO'):
			numResidues = len(pred)
			totalError = sum( [ p!=t for p, t in zip(pred, native) ] )
			tmpError = np.array([numResidues, totalError])

		elif response.startswith('SS'):
			numResidues = sum( [ m==0 for m in missing ] )
			totalError = sum( [ p!=t for p, t, m in zip(pred, native, missing) if m==0 ] )
			tmpError = np.array([numResidues, totalError])

		elif 'Phi' in response or 'Psi' in response:
			invalidResidues = [0] * len(missing)
			for i in xrange( len(missing) ):
				if missing[i]==1:
					invalidResidues[i]=1 
					if i > 0:
						invalidResidues[i-1]=1 
					if i < len(missing)-1:
						invalidResidues[i+1]=1

			invalidResidues[0] = 1
			invalidResidues[len(missing)-1] = 1
 
			numResidues = sum( [ m==0 for m in invalidResidues ] )
			err1 = abs( pred - native )
                	err2 = np.float32(2 * np.pi) - err1
                	err = np.minimum(err1, err2)
			totalError = np.sum( [ e for e, m in zip(err, invalidResidues ) if m==0 ], axis=0 )
			tmpError = np.array([numResidues] + list(totalError) )
		else:
			print 'The Evaluate function not implemented for response: ', response
			exit(1)

		if errors.has_key(response):
			errors[response].append( tmpError)
		else:
			errors[response] = [ tmpError ]

	## calculate average error
	avgerrors = dict()
	for response, err in errors.iteritems():
		avgerrors[response] = np.average(err)

	return avgerrors


##evaluate the prediction accuracy. predictions is a dictionary with protein names as the key
def EvaluatePropertyPrediction(predictions, nativefolder):

	from DataProcessor import LoadNativeLabels

	errors = dict()
	names = []
	for name, preds in predictions.iteritems():
		#print 'name=', name
		nativeLabels = LoadNativeLabels(name, nativefolder, preds.keys() )
		if nativeLabels is None:
			continue

		names.append(name)

		for response, pred in preds.iteritems():
			native = nativeLabels[Response2LabelName(response)]
			missing = nativeLabels['Missing']

			if response.startswith('DISO'):
				numResidues = len(pred)
				totalError = sum( [ p!=t for p, t in zip(pred, native) ] )
				tmpError = np.array([numResidues, totalError])

			elif response.startswith('SS'):
				numResidues = sum( [ m==0 for m in missing ] )
				totalError = sum( [ p!=t for p, t, m in zip(pred, native, missing) if m==0 ] )
				tmpError = np.array([numResidues, totalError])

			elif 'Phi' in response or 'Psi' in response:
				invalidResidues = [0] * len(missing)
				for i in xrange( len(missing) ):
					if missing[i]==1:
						invalidResidues[i]=1 
						if i > 0:
							invalidResidues[i-1]=1 
						if i < len(missing)-1:
							invalidResidues[i+1]=1

				invalidResidues[0] = 1
				invalidResidues[len(missing)-1] = 1
 
				numResidues = sum( [ m==0 for m in invalidResidues ] )
				err1 = abs( pred - native )
                		err2 = np.float32(2 * np.pi) - err1
                		err = np.minimum(err1, err2)
				totalError = np.sum( [ e for e, m in zip(err, invalidResidues ) if m==0 ], axis=0 )
				tmpError = np.array([numResidues] + list(totalError) )
			else:
				print 'The Evaluate function not implemented for response: ', response
				exit(1)

			if errors.has_key(response):
				errors[response].append( tmpError)
			else:
				errors[response] = [ tmpError ]

	## calculate average error
	avgErrPerTarget = dict()
	avgErrPerResidue = dict()
	allerrors = dict()
	for response, e in errors.iteritems():
		
		err = np.array(e)
		err_avg = np.average(err, axis=0)
		err2 = err_avg[1:]*1./err_avg[0]

		ind_err = np.divide(err[:, 1:]*1.0, err[:, 0:1])
		err1 = np.average( ind_err, axis=0)
		
		avgErrPerTarget[response] = err1
		avgErrPerResidue[response] = err2
		"""
		print '*********************Error for response ', response, '************************'
		print 'avg by target: ', err1, ' avg by residue: ', err2
		print '                            '
		print '*********************Individual Error for response ', response, '************************'
		"""
		allerrors[response] = dict()
		for name, e0 in zip(names, ind_err):
			##print name, e0	
			allerrors[response][name] = e0

	return avgErrPerTarget, avgErrPerResidue, allerrors


## convert predicted secondary structures, return a SS3 string
def SS8String2SS3(ss8):
	return ss8.replace('G','H').replace('I','H').replace('B','E').replace('T','L').replace('S','L')

## convert predicted probability, return both prob and a SS3 string
## ss8prob is a matrix of dimension L*8
def SS8Prob2SS3(ss8prob):
	assert ss8prob.shape[1] == 8

	len = ss8prob.shape[0]
	ss3prob = np.zeros( (len, 3), dtype=ss8prob.dtype )

	for ss8_index, ss3_index in zip(range(8), SS8Code2SS3Code):
		ss3prob[:,ss3_index] += ss8prob[:, ss8_index]

	ss3code = np.argmax(ss3prob, axis=1)

	ss3letters = [ SS3Code2Letter[c] for c in ss3code ]
	ss3str = ''.join(ss3letters)

	return ss3str, ss3prob
