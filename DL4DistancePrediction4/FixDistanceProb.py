import numpy as np
import os
import sys
import getopt

import ContactUtils
import DistanceUtils

def Usage():
	print 'python FixDistanceProb.py -i pkl_files'
	print '   You may provide a list of predictedDistMatrix.pkl as input. They shall be separated by ;'



def main(argv):

	try:
		opts, args = getopt.getopt(argv,"i:",["input="])
                print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)


        if len(opts) < 1:
                Usage()
                exit(1)

	inputFiles = None

        for opt, arg in opts:
                if opt in ("-i", "--input"):
			inputFiles = arg.split(" ;")
		else:
			Usage()
			exit(1)

	if inputFiles is None:
		print 'Please provide at least one input file'
		exit(1)

	for inputFile in inputFiles:
		if not os.path.isfile(inputFile):
			print 'The input file does not exist: ', inputFile
                	exit(1)

		targetName = os.path.basename(inputFile).split('.')[0]

        	content = DistanceUtils.LoadRawDistProbFile(inputFile)

        	if len(content) != 7:
			print 'WARNING: the input file format is not what expected. Nothing is done'
			continue

                name, sequence, predictedDistProb, predictedContactProb, labelWeight, labelDistribution, refDistProb = content

                fixedProb = dict()
		contactProb = dict()
                for apt in predictedDistProb.keys():
			labelName = Response2LabelName(apt)
			labelType = Response2LabelType(apt)

			if labelName.startswith('HB') or labelName.startswith('Beta') or (not labelType.startswith('Discrete')):
				fixedProb[apt] = predictedDistProb[apt]
				contactProb[apt] = predictedContactProb[apt]
				print 'WARNING: No bias correction is done on the response: ', apt

			elif labelType.startswith('Discrete'):
				subType = labelType[len('Discrete'): ]
                        	#print 'shapes: ', predictedDistProb[apt].shape, np.array(labelWeight[apt]).shape, np.array(labelDistribution[apt]).shape
                        	fixedProb[apt] = DistanceUtils.FixDistProb( predictedDistProb[apt], labelWeight[apt], labelDistribution[apt])
				labelOf8 = DistanceUtils.LabelsOfOneDistance(config.ContactDefinition, config.distCutoffs[subType])
				contactProb[labelName] = ContactUtils.Distance2Contact(fixedProb, labelOf8)	


		## write the results
		savefile = targetName + '.fixedDistMatrix.pkl'
		fh = open(savefile, 'w')
		cPickle.dump( (name, sequence, fixedProb, contactProb, None, labelDistribution, refDistProb), fh, protocol = cPickle.HIGHEST_PROTOCOL)
		fh.close()


if __name__ == "__main__":
	main(sys.argv[1:])
