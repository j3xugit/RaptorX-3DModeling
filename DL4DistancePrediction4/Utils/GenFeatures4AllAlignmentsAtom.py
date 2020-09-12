import cPickle
import os
import sys
import getopt

"""
This script generates feature files for a batch of alignments
It calls GenFeatures4SingleAlignmentAtom.py
"""
##import SelectTemplatesFromRankFile
import GenFeatures4SingleAlignmentAtom

sys.path.append('../../Common/')
import LoadTPLTGT

sys.path.append('../../Alignment/Scripts')
from LoadAlignments import LoadAlignments

TGTDir = '/mnt/data/RaptorXCommon/TrainTestData/Distance_Contact_TrainData/7952_Old_PDB25_Training_Data/training_tgt'
TPLDir = '/mnt/data/RaptorXCommon/Server_Database/TPL_BC100'
PDBDir = '/mnt/data/RaptorXCommon/TrainTestData/Distance_Contact_TrainData/Jinbo_Folder/data4HBBeta/TPLBC100AtomCoordinates/'
RANKDir = '/mnt/data/RaptorXCommon/TrainTestData/Ranking_Project/raw_data/7952_against_bc40_rank'

"""
TGTDir = '/mnt/data/RaptorXCommon/TrainTestData/Ranking_Project/17w_data/17W_single_tgt/'
TPLDir = '/mnt/data/RaptorXCommon/TrainTestData/Ranking_Project/17w_data/17W_single_tpl/'
"""

def load_protein_list(list_file):
	fin = open(list_file, 'r')
	lines = [line.strip() for line in fin.readlines()]
	fin.close()
	return lines

"""
def LoadAlignments(alignment_file, tpl_dir, tgt_dir):

	alignments = []

        fin = open(alignment_file, 'r')
        content = [ line.strip() for line in list(fin) ]
        fin.close()

	##remove empty lines
	alignment_result = [ c for c in content if c ]

	if len(alignment_result)%4 !=0 :
		print 'the number of lines in the alignment file is incorrect: ', alignment_file
		exit(-1)


	templates = set()
	querys = set()

        for i in range(0, len(alignment_result), 4):
                tpl_name = alignment_result[i][1:]
                tpl_seq = alignment_result[i+1]
                tgt_name = alignment_result[i+2][1:]
                tgt_seq = alignment_result[i+3]

		##check to see if the tpl file exists or not
		tpl_file = os.path.join(tpl_dir, tpl_name+'.tpl')
		if not os.path.isfile(tpl_file):
			print 'WARNING: the tpl file %s does not exist!' % tpl_file
			continue

		##check to see if the tgt file exists or not
		tgt_file = os.path.join(tgt_dir, tgt_name+'.tgt')
		if not os.path.isfile(tgt_file):
			print 'WARNING: the tgt file %s does not exist!' % tgt_file
			continue

                alignments.append((tpl_seq, tgt_seq, tpl_name, tgt_name))

		if tpl_name not in templates:
			templates.add(tpl_name)

		if tgt_name not in querys:
			querys.add(tgt_name)

	print 'In total loaded ', len(alignments), ' alignments involving ', len(templates), ' templates and ', len(querys), ' query sequences'

	## load tgt and tpl
	tgtPool = {}
	tplPool = {}


	for tpl in templates:
		tpl_file = os.path.join(tpl_dir, tpl+'.tpl')
		if not os.path.isfile(tpl_file):
			print 'the tpl file %s does not exist!' % tpl_file
			exit(-1)

		tplPool[tpl] = LoadTPLTGT.load_tpl(tpl_file)

	for tgt in querys:
		tgt_file = os.path.join(tgt_dir, tgt+'.tgt')
		if not os.path.isfile(tgt_file):
			print 'the tgt file %s does not exist!' % tgt_file
			exit(-1)

		tgtPool[tgt] = LoadTPLTGT.load_tgt(tgt_file)
	
	return alignments, tplPool, tgtPool	
"""

###here we always assume that we use the pdb file fixed by Sheng Wang's tool.
###that is, the pdb file shall be consistent with the tpl file, otherwise this procedure may generate wrong results
###we only read the coordinates of Ca, Cb, Cg, N, and O
def LoadAtomCoordinates(atomCoord_file):
	fh = open(atomCoord_file, 'r')
	data = cPickle.load(fh)
	fh.close()
	coords = []
	for d in data:
		coord = dict()
		for k, v in d.items():
			if k != 'valid':
				if v is not None:
					coord[k] = v.get_array()
				else:
					coord[k] = None
		coords.append(coord)

	return coords

def Usage():
	print 'python GenFeatures4AllAlignmentsAtom.py -i alignment_file [ -q tgt_dir ] [ -t tpl_dir ] [ -p atom_coord_dir]'
	print '\r-i: alignment file in fasta format, template is placed before query '
	print '\r-q: the tgt folder, default: ' + TGTDir
	print '\r-t: the tpl folder, default: ' + TPLDir
	print '\r-p: the folder containing all the tpl atom coordinate files extracted by ExtractCoordinatesFromTPLPDB.py, default: ' + PDBDir

def main(argv):

	alignment_file = None
	tgt_dir = TGTDir
	tpl_dir = TPLDir
	pdb_dir = PDBDir
	
	try:
        	opts, args = getopt.getopt(argv,"i:q:t:p:",["input=", "sequenceDir=","templateDir=", "tplpdbDir="])
        	print opts, args
    	except getopt.GetoptError:
        	Usage()
        	exit(-1)


    	if len(opts) < 1:
        	Usage()
        	exit(-1)

    	for opt, arg in opts:
		if opt in ("-i", "--input"):
			alignment_file = arg
		elif opt in ("-q", "--sequenceDir"):
			tgt_dir = arg
		elif opt in ("-t", "--templateDir"):
			tpl_dir = arg
		elif opt in ("-p", "--tplpdbDir"):
			pdb_dir = arg
		else:
			print Usage()
			exit(-1)

	if alignment_file is None or not os.path.isfile(alignment_file):
		print 'please provide a valid alignment file for loading'
		exit(-1)

	alignments, tplPool, tgtPool = LoadAlignments(alignment_file, tpl_dir, tgt_dir)

	invalidTemplates = []
	features = []
	for aln in alignments:
		## aln has 4 items: tpl_seq, tgt_seq, tpl_name, tgt_name
		atomCoordFile = os.path.join(pdb_dir, aln[2]+'.atomCoordinates.pkl')
		if not os.path.isfile(atomCoordFile):
			print 'WARNING: the template coordinate file does not exist: ', atomCoordFile
			invalidTemplates.append(aln[2])
			continue

		tplAtomCoordinates = LoadAtomCoordinates( atomCoordFile )
		feature, distMatrix, xyresidues, coordinates = GenFeatures4SingleAlignmentAtom.GenFeature4Alignment(aln, tplPool[aln[2]], tgtPool[aln[3]], tplAtomCoordinates)
		#features.append( (aln[3], aln[2], feature, distMatrix, xyresidues, coordinates))
		features.append( (aln[3], aln[2], feature, distMatrix))

	##write template-based features to the pkl file
	resultfile = os.path.basename(alignment_file).split('.')[0] + '.tplSimFeaturesAtom.pkl'

	fh = open(resultfile, 'wb')
	cPickle.dump( features, fh, protocol=cPickle.HIGHEST_PROTOCOL)
	fh.close()

	print '#invalid templates: ', len(invalidTemplates)
	print '\n'.join(invalidTemplates)


if __name__ == '__main__':
	main(sys.argv[1:])
	
