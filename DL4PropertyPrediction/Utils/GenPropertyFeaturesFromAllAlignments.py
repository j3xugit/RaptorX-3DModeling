import cPickle
import os
import sys
import getopt

from GenPropertyFeaturesFromSingleAlignment import GenPropertyFeaturesFromSingleAlignment
from Alignment.LoadAlignments import LoadAlignments

TGTDir = '/mnt/data/RaptorXCommon/TrainTestData/Distance_Contact_TrainData/7952_Old_PDB25_Training_Data/training_tgt'
TPLDir = '/mnt/data/RaptorXCommon/Server_Database/TPL_BC100'

def load_protein_list(list_file):
	fin = open(list_file, 'r')
	lines = [line.strip() for line in fin.readlines()]
	fin.close()
	return lines

def Usage():
	#print 'python GenPropertyFeaturesFromAllAlignments.py -i alignment_file [ -q tgt_dir ] [ -t tpl_dir ] [ -a tpl_angle_file_dir ]'
	print 'python GenPropertyFeaturesFromAllAlignments.py -i alignment_file [ -q tgt_dir ] [ -t tpl_dir ]'
	print '	This script generates template-based features for property prediction from a set of alignments'
	print '\r-i: an alignment file in FASTA format, which may contain multiple pairwise alignments and in each alignment, template is placed before query'
	print '\r-q: the tgt folder, each file shall end with .tgt or .tgt.pkl, default ' + TGTDir
	print '\r-t: the tpl folder, each file shall end with .tpl.pkl, default ' + TPLDir
	print '\r The resultant file is named XXX.tplPropertyFeatures.pkl where XXX is the basename of alignment_file'
	print '\r tplPropertyFeatures.pkl is a list of tuples. Each tuple has 4 items: tempName, queryName, similarityFeatures, tempProperties'

def main(argv):

	alignment_file = None
	tgt_dir = TGTDir
	tpl_dir = TPLDir
	#ang_dir = AngleDir
	
	try:
        	opts, args = getopt.getopt(argv,"i:q:t:a:",["input=", "sequenceDir=","templateDir=", "tplAngleDir="])
        	opts, args = getopt.getopt(argv,"i:q:t:",["input=", "sequenceDir=","templateDir="])
        	print opts, args
    	except getopt.GetoptError:
        	Usage()
        	exit(1)

    	if len(opts) < 1:
        	Usage()
        	exit(1)

    	for opt, arg in opts:
		if opt in ("-i", "--input"):
			alignment_file = arg
		elif opt in ("-q", "--sequenceDir"):
			tgt_dir = arg
		elif opt in ("-t", "--templateDir"):
			tpl_dir = arg
		"""
		elif opt in ("-a", "--tplAngleDir"):
			ang_dir = arg
		"""
		else:
			Usage()
			exit(1)

	if alignment_file is None or not os.path.isfile(alignment_file):
		print 'ERROR: invalid alignment file:', alignment_file
		exit(1)

	alignments, tplPool, tgtPool = LoadAlignments(alignment_file, tpl_dir, tgt_dir)

	#invalidTemplates = []
	features = []
	for aln in alignments:
		## aln has 4 items: tpl_seq, tgt_seq, tpl_name, tgt_name
		feature, xyresidues, tplProperties = GenPropertyFeaturesFromSingleAlignment(aln, tplPool[aln[2]], tgtPool[aln[3]])

		## aln[3] and aln[2] are query and template names, respectively
		features.append( (aln[3], aln[2], feature, tplProperties))

	##write template-based features to a pkl file
	resultfile = os.path.basename(alignment_file).split('.')[0] + '.tplPropertyFeatures.pkl'

	with open(resultfile, 'wb') as fh:
		cPickle.dump( features, fh, protocol=cPickle.HIGHEST_PROTOCOL)

	#print '#invalid templates: ', len(invalidTemplates)
	#print '\n'.join(invalidTemplates)

if __name__ == '__main__':
	main(sys.argv[1:])
	
