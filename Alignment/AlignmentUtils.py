import os
import sys
import numpy as np
import scipy
import cPickle

import SimilarityScore

from Common import SequenceUtils
from Common import PDBUtils
from Common.PDBUtils import ExtractCoordinatesBySeq, CalcDistMatrix, CalcTwoROriMatrix

from Common.LoadTPLTGT import load_hhm as LoadHHM
from Common.LoadTPLTGT import load_tgt as LoadTGT
from Common.LoadTPLTGT import load_tpl as LoadTPL

"""
This script generates a set of features for a query-template alignment and also a distance/orientation matrix for target by copying coordinates from template
"""

## this function returns an alignment, which is a tuple. The first entry in the tuple is the template sequence in alignment including gaps;
## the second entry is the query sequence in alignment including gaps.
## if returnNames is True, this function also returns a tuple (templateName, queryName)

## tpl is the template object and tgt is the query object
## tgt and tpl are Python dict() and shall have the following keys: length, sequence, name, PSSM, PSFM and optionally SS

## the template sequence in alignment (after removing gaps) shall be a substring of the tpl sequence
## the target sequence in alignment (after removing gaps) shall be exactly same as the tgt sequence

## read one alignment. templateFirst=True indicates template sequence is placed before query sequence
## seqTemplatePair specifies the pair of query seq and template names
## queryName specifies the query sequence name
## a user just needs to specify one of templateFirst, seqTemplatePair and queryName
def ReadAlignment(alnfile, templateFirst=True, seqTemplatePair=None, queryName=None, returnNames=False):
	with open(alnfile, 'r') as fh:
		content = [ line.strip() for line in list(fh) ]
	content = [ line for line in content if (not line.startswith('#') and len(line)>1) ]

	if len(content) < 4:
                print "ERROR: the pairwise alignment file %s shall have at least 4 lines" % fasta
                exit(1)

	line = content[0]
        ##the first non-empty line shall starst with > followed by a protein name
        if line[0] != '>' :
                print 'ERROR: incorrect format in the fasta file at line: ', line
                exit(1)

	firstName = line[1:].split()[0]
        firstSeq = ""
	
        ##read the sequence for this first protein. the sequence could be in several lines
        i = 1
        line = content[i]
        while not line.startswith('>'):
                firstSeq += line
                i += 1
                line = content[i]

        ## get the name of the 2nd protein
        secondName = line[1:].split()[0]

        ## get the sequence of the 2nd protein
        secondSeq = ''.join(content[i+1:])

        if len(firstSeq) != len(secondSeq):
                print "ERROR: inconsistent query and template sequence length (including gaps) in alignment file", alnfile
                exit(1)

	if seqTemplatePair is None:
		if queryName is not None:
			if firstName == queryName:
				if returnNames:
					return (secondSeq, firstSeq), (secondName, firstName)
				return (secondSeq, firstSeq)

			elif secondName == queryName:
				if returnNames:
					return (firstSeq, secondSeq), (firstName, secondName)
				return (firstSeq, secondSeq)
			else:
				print 'ERROR: targetName does not match any names in alignment file: ', queryName, alnfile
				exit(1)

		elif templateFirst:
			if returnNames:
				return (firstSeq, secondSeq), (firstName, secondName)
			return (firstSeq, secondSeq)
		else:
			if returnNames:
				return (secondSeq, firstSeq), (secondName, firstName)
			return (secondSeq, firstSeq)
	else:
		queryName = seqTemplatePair[0]
		templateName = seqTemplatePair[1]
		if queryName in firstName and templateName in secondName:
			return (secondSeq, firstSeq), (templateName, queryName)
		else:
			return (firstSeq, secondSeq), (templateName, queryName)

	"""
 	## remove gaps from the aligned sequences
        firstRealSeq = firstSeq.translate(None,'-')
        secondRealSeq = secondSeq.translate(None, '-')


        firstBits = [ c != '-' for c in firstSeq ]
        secondBits = [ c != '-' for c in secondSeq ]

        firstIndex = np.cumsum(firstBits)
        secondIndex = np.cumsum(secondBits)

        ##build the index mapping, index starts from 1 and for gaps it becomes negative
        alignment = []
        for c1, c2, i, j in zip(tmp, seq, tmp_index, seq_index):
                if c1 == '-':
                        index1 = -i
                else:
                        index1 = i

                if c2 == '-':
                        index2 = -j
                else:
                        index2 = j

                alignment.append( (index1, index2) )

        #print alignment

        return alignment, nam1_content, nam2_content, nam1_full, nam2_full, nam1, nam2
	"""

##remove gaps in the alignment and return the template and query sequences without gaps
def ExtractSeqFromAlignment(alignment):
	template = alignment[0].translate(None, '-')
	query = alignment[1].translate(None, '-')
	return (template, query)
	
## alignment[0] and alignment[1] corresponds to template and query sequence, respectively
## tplAtomCoordinates shall have the 3D coordinate information for the whole template sequence
## tplAtomCoordinates have the same format as coordinates used in the above function
def CopyTemplateCoordinates(alignment, tpl, tgt, tplAtomCoordinates):	
	##check consistency between alignment, tpl and tgt
	if len(alignment[0]) != len(alignment[1]):
		print 'ERROR: the length of query and template in alignment is inconsistent'
		exit(1)

	template, query = ExtractSeqFromAlignment(alignment)

	## template and query shall be substrings of tpl and tgt, respectively
	tpl_start = tpl['sequence'].find(template)
	if tpl_start == -1:
		print 'ERROR: the template sequence in alignment is not a substring of sequence in template', tpl['name']
		exit(1)

	tgt_start = tgt['sequence'].find(query)
	if tgt_start == -1:
		print 'ERROR: the query sequence in alignment is not a substring of sequence in target', tgt['name']
		exit(1)

	##wrong if tgt_start is not 0, here we require that the query sequence in alignment is exactly same as that in tgt
	assert (tgt_start == 0)

	##index for tgt and tpl, respectively
	tgt_pos = tgt_start
	tpl_pos = tpl_start

	##effective coordinates copied from aligned template positions
	copiedCoordinates = [ None ] * tgt['length']

	for al_pos in range(len(alignment[0])):

		if alignment[0][al_pos] == '-' and alignment[1][al_pos] == '-':
			print 'WARNING: there shall not be two gaps at any aligned positions'
			continue

		## there is a gap in template, i.e., an insertion in query
		if alignment[0][al_pos] == '-':
			tgt_pos += 1
			continue

		## if there is a gap in query, just skip it
		if alignment[1][al_pos] == '-':
			tpl_pos += 1
			continue

		## copy 3D coordinates from template. It is possible that the template coordinate = None
		copiedCoordinates[tgt_pos] = tplAtomCoordinates[tpl_pos]
			
		tpl_pos += 1
		tgt_pos += 1

	return copiedCoordinates

## create matrix for query sequence by copying from templateMatrix based upon seq2template mapping
## templateMatrix: python dict() for a set of template matrices
## seq2templateMapping is a tuple of two entries. Each is a list of residue indices. The first is for seq and the 2nd for template. These two lists shall have same length
## for a sequence matrix, if it does not have corresponding entry in template, then its value is set to an invalid value

## matrix type could be Orientation or Distance
def CopyTemplateMatrix(seqLen, seq2templateMapping, templateMatrices):
	seqIndices, templateIndices = seq2templateMapping
	assert len(seqIndices) == len(templateIndices)
	assert len(seqIndices) <= seqLen

	seqMatrix_x = [ e for e in seqIndices for _ in range(len(seqIndices) ) ]
	seqMatrix_y = seqIndices * len(seqIndices)

	tempMatrix_x = [ e for e in templateIndices for _ in range(len(templateIndices) ) ]
	tempMatrix_y = templateIndices * len(templateIndices)

	seqDistMatrices = dict()
	seqOriMatrices = dict()

	tempDistMatrices, tempOriMatrices = templateMatrices

	for k, v in tempDistMatrices.iteritems():
		#seqMatrix = np.ones( (seqLen, seqLen), dtype=np.float16) * PDBUtils.InvalidDistance
		seqMatrix = np.full( (seqLen, seqLen), PDBUtils.InvalidDistance, dtype=np.float16)
		seqMatrix[ (seqMatrix_x, seqMatrix_y) ] = v[ (tempMatrix_x, tempMatrix_y) ]
		seqDistMatrices[k] = seqMatrix

	for k, v in tempOriMatrices.iteritems():
		#seqMatrix = np.ones( (seqLen, seqLen), dtype=np.float16) * PDBUtils.InvalidDegree
		seqMatrix = np.full( (seqLen, seqLen), PDBUtils.InvalidDegree, dtype=np.float16)
		seqMatrix[ (seqMatrix_x, seqMatrix_y) ] = v[ (tempMatrix_x, tempMatrix_y) ]
		seqOriMatrices[k] = seqMatrix

	return (seqDistMatrices, seqOriMatrices)

## score an alignment and also copy coordinates from template
def GetTemplateMatrixByAlignment(alignment, tpl, tgt, tplMatrices):	
	##check consistency between alignment, tpl and tgt
	if len(alignment[0]) != len(alignment[1]):
		print 'ERROR: the query and template (including gaps) in alignment have inconsistent length'
		exit(1)

	template, query = ExtractSeqFromAlignment(alignment)

	## template and query shall be substrings of tpl and tgt, respectively
	tpl_start = tpl['sequence'].find(template)
	if tpl_start == -1:
		print 'ERROR: the template sequence in alignment is not a substring of the sequence in tpl'
		exit(1)

	tgt_start = tgt['sequence'].find(query)
	if tgt_start == -1:
		print 'ERROR: the query sequence in alignment is not a substring of the sequence in tgt'
		exit(1)

	##wrong if tgt_start is not 0, here we require that the query sequence in alignment is exactly same as that in tgt
	assert (tgt_start == 0)

	##index for tgt and tpl, respectively
	tgt_pos = tgt_start
	tpl_pos = tpl_start

	seq2templateMapping = [ [], [] ]

	## alignment[0] represents the template sequence with gaps
	for al_pos in range(len(alignment[0])):

		if alignment[0][al_pos] == '-' and alignment[1][al_pos] == '-':
			print 'WARNING: there shall not be two gaps at any aligned positions'
			continue

		## there is a gap in template, i.e., an insertion in query
		if alignment[0][al_pos] == '-':
			tgt_pos += 1
			continue

		## if there is a gap in query, just skip it
		if alignment[1][al_pos] == '-':
			tpl_pos += 1
			continue

		seq2templateMapping[0].append(tgt_pos)
		seq2templateMapping[1].append(tpl_pos)

		tpl_pos += 1
		tgt_pos += 1

	copiedMatrix = CopyTemplateMatrix(tgt['length'], seq2templateMapping, tplMatrices)
	return copiedMatrix


## score an alignment and also copy coordinates from template
def ScoreAlignment(alignment, tpl, tgt, tplAtomCoordinates=None, tplMatrices=None, UseSS3=False, UseSS8=False, UseACC=False, UseEnv=False, debug=False):	
	##check consistency between alignment, tpl and tgt
	if len(alignment[0]) != len(alignment[1]):
		print 'ERROR: the query and template (including gaps) in alignment have inconsistent length'
		exit(1)

	template, query = ExtractSeqFromAlignment(alignment)

	## template and query shall be substrings of tpl and tgt, respectively
	tpl_start = tpl['sequence'].find(template)
	if tpl_start == -1:
		print 'ERROR: the template sequence in alignment is not a substring of the sequence in template', tpl['name']
		print 'TPLstring:', tpl['sequence']
		print 'TPLinAli :', template
		exit(1)

	tgt_start = tgt['sequence'].find(query)
	if tgt_start == -1:
		print 'ERROR: the query sequence in alignment is not a substring of the sequence in query', tgt['name']
		exit(1)

	##wrong if tgt_start is not 0, here we require that the query sequence in alignment is exactly same as that in tgt
	assert (tgt_start == 0)

	##index for tgt and tpl, respectively
	tgt_pos = tgt_start
	tpl_pos = tpl_start

	## a flag vector indicating insertion in query, indicated by a flag 1
	insertX = np.zeros( (tgt['length'], 1), dtype=np.uint8)

	## At each query position, there are numFeatures basic features, in which 10 are purely derived from sequence (profiles)
	numFeatures = 10
	if UseSS3:
		numFeatures += 1
	if UseSS8:
		numFeatures += 1
	if UseACC:
		numFeatures += 1
	if UseEnv:
		numFeatures += 2

	sequentialFeatures = np.zeros( (tgt['length'], numFeatures), dtype=np.float16)

	##effective coordinates copied from aligned template positions
	if tplAtomCoordinates is not None:
		copiedCoordinates = [ None ] * tgt['length']
	else:
		copiedCoordinates = None

	##for debug only
	if debug:
		XYresidues = - np.ones( (tgt['length'], 4), dtype=np.int16)
	else:
		XYresidues = None

	tgtAAIndex = SequenceUtils.Seq2OrderOf3LetterCode(tgt['sequence'])
	tplAAIndex = SequenceUtils.Seq2OrderOf3LetterCode(tpl['sequence'])

	seq2templateMapping = [ [], [] ]

	## alignment[0] represents the template sequence with gaps
	for al_pos in range(len(alignment[0])):

		if alignment[0][al_pos] == '-' and alignment[1][al_pos] == '-':
			print 'WARNING: there shall not be two gaps at any aligned positions'
			continue

		## there is a gap in template, i.e., an insertion in query
		if alignment[0][al_pos] == '-':
			## need to generate some flag features for insertion in query
			insertX[tgt_pos] = 1

			if debug:
				XYresidues[tgt_pos][0]= tgt_pos
				XYresidues[tgt_pos][1] = ord(tgt['sequence'][tgt_pos]) - ord('A')

			tgt_pos += 1
			continue

		## if there is a gap in query, just skip it
		if alignment[1][al_pos] == '-':
			tpl_pos += 1
			## no need to generate flag features for insertion in template
			continue

		if debug:
			XYresidues[tgt_pos][0]= tgt_pos
			XYresidues[tgt_pos][1] = ord(tgt['sequence'][tgt_pos]) - ord('A')
			XYresidues[tgt_pos][2]= tpl_pos
			XYresidues[tgt_pos][3] = ord(tpl['sequence'][tpl_pos]) - ord('A')

		seq2templateMapping[0].append(tgt_pos)
		seq2templateMapping[1].append(tpl_pos)

		tAA = tplAAIndex[tpl_pos]
		sAA = tgtAAIndex[tgt_pos]
		
		seq_Id = int(tAA == sAA)
		blosum80 = SimilarityScore.newBLOSUM80[tAA, sAA]
		blosum62 = SimilarityScore.newBLOSUM62[tAA, sAA]
		blosum45 = SimilarityScore.newBLOSUM45[tAA, sAA]
		cc = SimilarityScore.newCC50[tAA, sAA]      
		hdsm = SimilarityScore.newHDSM[tAA, sAA]

		x, y = tpl_pos, tgt_pos
		spScore = SimilarityScore.MutationOf2Pos6(x, y, tpl, tgt)
		spScore_ST = SimilarityScore.MutationOf2Pos6_ST(x, y, tpl, tgt)
		pmScore = SimilarityScore.MutationOf2Pos5(x, y, tpl, tgt)
		pmScore_ST = SimilarityScore.MutationOf2Pos5_ST(x, y, tpl, tgt)

		point_feature = [seq_Id, blosum80, blosum62, blosum45, spScore, spScore_ST, pmScore, pmScore_ST, cc, hdsm]

		if UseSS3:
			SS3 = SimilarityScore.SSMutationScore_3State(x, y, tpl, tgt)
			point_feature.append(SS3)

		if UseSS8:
			SS8 = SimilarityScore.SSMutationScore_6State(x, y, tpl, tgt)
			point_feature.append(SS8)

		if UseACC:
			ACC = SimilarityScore.ACCMutationScore_3State(x, y, tpl, tgt)
			point_feature.append(ACC)

		if UseEnv:
			envScore = SimilarityScore.SingletonScore_ProfileBased(x, y, tpl, tgt)
			wsEnvScore = SimilarityScore.SingletonScore_WS(x, y, tpl, tgt)
			point_feature.append(envScore)
			point_feature.append(wsEnvScore)

		sequentialFeatures[tgt_pos]=np.array(point_feature).astype(np.float32)

		## copy 3D coordinates from template. It is possible that the template coordinate = None
		if tplAtomCoordinates is not None:
			copiedCoordinates[tgt_pos] = tplAtomCoordinates[tpl_pos]

		tpl_pos += 1
		tgt_pos += 1

	if tplMatrices is not None:
		copiedMatrix = CopyTemplateMatrix(tgt['length'], seq2templateMapping, tplMatrices)
		return sequentialFeatures, insertX, copiedMatrix, XYresidues

	return sequentialFeatures, insertX, copiedCoordinates, XYresidues

##this function generates features for those query positions with effective aligned template positions
##do we need to deal with HIS tags?
def GenFeature4Alignment(alignment, tpl, tgt, tplAtomCoordinates=None, tplMatrices=None):
	sequentialFeatures, insertX, copied, XYresidues = ScoreAlignment(alignment, tpl, tgt, tplAtomCoordinates=tplAtomCoordinates, tplMatrices=tplMatrices)

	if tplMatrices is not None:
		copiedDistMatrix, copiedOriMatrix = copied
		return np.hstack([insertX, sequentialFeatures]), copiedDistMatrix, copiedOriMatrix

	if tplAtomCoordinates is not None:
		copiedDistMatrix = CalcDistMatrix(copied)
		copiedOriMatrix = CalcTwoROriMatrix(copied)
		#copiedOriMatrix = None
		return np.hstack([insertX, sequentialFeatures]), copiedDistMatrix, copiedOriMatrix

	return np.hstack([insertX, sequentialFeatures])

def GenNullAlignmentFeatures(seqLen, UseSS3=False, UseSS8=False, UseACC=False, UseEnv=False):
	## a flag vector indicating insertion in query, indicated by a flag 1
        insertX = np.zeros( (seqLen, 1), dtype=np.uint8)

        ## At each query position, there a set of numFeatures basic features, n which 10 are purely derived from sequence and sequence profiles.
        numFeatures = 10
        if UseSS3:
                numFeatures += 1
        if UseSS8:
                numFeatures += 1
        if UseACC:
                numFeatures += 1
        if UseEnv:
                numFeatures += 2

        sequentialFeatures = np.zeros( (seqLen, numFeatures), dtype=np.float16)
	sequentialFeatures = np.hstack([insertX, sequentialFeatures])

	feature = dict()
	feature['tplSimScore'] = sequentailFeatures
	feature['tplDistMatrix'] = None
	feature['tplOriMatrix'] = None
	feature['template'] = None

	return feature
		

def GenerateAlignmentFeatures(seqTemplatePair, aliDir=None, tgtDir=None, tplDir=None, pdbDir=None, modelSpecs=None):

	query = seqTemplatePair[0]
	template = seqTemplatePair[1]

	alnfile = template + '-' + query + '.fasta' 
	alnfile = os.path.join(aliDir, alnfile)
	alignment = ReadAlignment(alnfile, seqTemplatePair)

	#print 'alnfile: ', alnfile
	#print alignment

	tgtFile = os.path.join( tgtDir, query + '.hhm')
	if not os.path.isfile(tgtFile):
		tgtFile = os.path.join( tgtDir, query + '.tgt')
		if not os.path.isfile(tgtFile):
			print 'ERROR: cannot find .hhm or .tgt file for the query protein: ', query, ' in folder ', tgtDir
			exit(1)
		else:
			tgt = LoadTGT(tgtFile)
	else:
		tgt = LoadHHM(tgtFile)

	#print 'tgtFile: ', tgtFile

	tplFile = os.path.join( tplDir, template + '.hhm')
	if not os.path.isfile(tplFile):
		tplFile = os.path.join( tplDir, query + '.tpl')
		if not os.path.isfile(tplFile):
			print 'ERROR: cannot find .hhm or .tpl file for the template protein: ', template, ' in folder ', tplDir
			exit(1)
		else:
			tpl = LoadTPL(tplFile)
	else:
		tpl = LoadHHM(tplFile)

	#print 'tplFile: ', tplFile
	
	pdbfile = os.path.join(pdbDir, template + '.atomCoordinates.pkl')	
	if not os.path.isfile(pdbfile):
		print 'ERROR: cannot find the PDB file for template: ', template, ' in folder ', pdbDir
		exit(1)
	with open(pdbfile, 'rb') as fh:
		tplAtomCoordinates = cPickle.load(fh)['coordinates']

	simScore, distMatrix, oriMatrix = GenFeature4Alignment(alignment, tpl, tgt, tplAtomCoordinates)
	feature = dict()
	feature['tplSimScore'] = simScore
	feature['tplDistMatrix'] = distMatrix
	feature['tplOriMatrix'] = oriMatrix
	feature['template'] = tpl['name']

	return feature

def GetTemplateMatrix(seqTemplatePair, queryData, aliDir, tplDir, modelSpecs=None):

	query, template = seqTemplatePair[:2]

	alnfile = template + '-' + query + '.fasta' 
	alnfile = os.path.join(aliDir, alnfile)
	#print 'alnfile: ', alnfile

	alignment = ReadAlignment(alnfile, seqTemplatePair)
	#print alignment

	tplFile = os.path.join(tplDir, template + '.tpl.pkl')
	if not os.path.isfile(tplFile):
		print 'ERROR: cannot find .hhm or .tpl file for the template protein: ', template, ' in folder ', tplDir
		exit(1)

	#print 'tplFile: ', tplFile
	with open(tplFile, 'rb') as fh:
		tpl = cPickle.load(fh)

	templateMatrices = (tpl['atomDistMatrix'], tpl['atomOrientationMatrix'])	
	distMatrix, oriMatrix = GetTemplateMatrixByAlignment(alignment, tpl, tgt=queryData, tplMatrices=templateMatrices)
	feature = dict()
	#feature['tplSimScore'] = simScore
	feature['tplDistMatrix'] = distMatrix
	feature['tplOriMatrix'] = oriMatrix
	feature['template'] = tpl['name']

	return feature

def GenerateAlignmentFeatures2(seqTemplatePair, queryData, aliDir=None, tplDir=None, modelSpecs=None):

	query, template = seqTemplatePair[:2]

	alnfile = template + '-' + query + '.fasta' 
	alnfile = os.path.join(aliDir, alnfile)
	#print 'alnfile: ', alnfile

	alignment = ReadAlignment(alnfile, seqTemplatePair)
	#print alignment

	tplFile = os.path.join( tplDir, template + '.tpl.pkl')
	if not os.path.isfile(tplFile):
		print 'ERROR: cannot find .hhm or .tpl file for the template protein: ', template, ' in folder ', tplDir
		exit(1)

	#print 'tplFile: ', tplFile
	with open(tplFile, 'rb') as fh:
		tpl = cPickle.load(fh)

	templateMatrices = (tpl['atomDistMatrix'], tpl['atomOrientationMatrix'])	
	simScore, distMatrix, oriMatrix = GenFeature4Alignment(alignment, tpl, tgt=queryData, tplMatrices=templateMatrices)
	feature = dict()
	feature['tplSimScore'] = simScore
	feature['tplDistMatrix'] = distMatrix
	feature['tplOriMatrix'] = oriMatrix
	feature['template'] = tpl['name']

	return feature

def GenerateAlignmentFeatures3(queryData, aliFile, tplFile, queryName=None, modelSpecs=None):

	if queryName is not None:
		targetName = queryName
	else:
		targetName = queryData['name']

	alignment = ReadAlignment(aliFile, queryName=targetName)
	#print alignment

	if not os.path.isfile(tplFile):
		print 'ERROR: invalid template file ', tplFile
		exit(1)

	with open(tplFile, 'rb') as fh:
		tpl = cPickle.load(fh)

	templateMatrices = (tpl['atomDistMatrix'], tpl['atomOrientationMatrix'])	
	simScore, distMatrix, oriMatrix = GenFeature4Alignment(alignment, tpl, tgt=queryData, tplMatrices=templateMatrices)
	feature = dict()
	feature['tplSimScore'] = simScore
	feature['tplDistMatrix'] = distMatrix
	feature['tplOriMatrix'] = oriMatrix
	feature['template'] = tpl['name']

	return feature

def GenerateAlignmentFeatures4(queryData, aliFile, tplFolder, queryName=None, modelSpecs=None):

	if queryName is not None:
                targetName = queryName
        else:
                targetName = queryData['name']

        alignment, names = ReadAlignment(aliFile, queryName=targetName, returnNames=True)
	template = names[0]

        if not os.path.isdir(tplFolder):
                print 'ERROR: invalid template folder', tplFolder
                exit(1)

	tplFile = os.path.join(tplFolder, template + '.tpl.pkl')
	if not os.path.isfile(tplFile):
		print 'ERROR: invalid template file ', tplFile
		return None

        with open(tplFile, 'rb') as fh:
                tpl = cPickle.load(fh)

        templateMatrices = (tpl['atomDistMatrix'], tpl['atomOrientationMatrix'])
        simScore, distMatrix, oriMatrix = GenFeature4Alignment(alignment, tpl, tgt=queryData, tplMatrices=templateMatrices)
        feature = dict()
        feature['tplSimScore'] = simScore
        feature['tplDistMatrix'] = distMatrix
        feature['tplOriMatrix'] = oriMatrix
        feature['template'] = tpl['name']

        return feature


## parse an alignment file containing a set of pairwise alignments into many small files. Each file has only a single pairwise alignment
## each pairwise alignment has 4 consecutive lines, in which template information is placed before query sequence
def ParseAlignmentFile(alignmentFile, savefolder=os.getcwd()):

        with open(alignmentFile, 'r') as fin:
        	content = [ line.strip() for line in list(fin) ]
        alignments = [ c for c in content if c ]

        if len(alignments)%4 !=0 :
                print 'ERROR: the number of valid lines in the alignment file is not a multiple of 4:', alignment_file
                exit(1)

        for i in range(0, len(alignments), 4):
		assert alignments[i][0] == '>'
		assert alignments[i+2][0] == '>'
                tplname = alignments[i][1:].split()[0]
                tgtname = alignments[i+2][1:].split()[0]
		savefile = tplname + '-' + tgtname + '.fasta'
		savefile = os.path.join(savefolder, savefile)
		alignStr = '\n'.join(alignments[i:i+4])
		with open(savefile, 'w') as fh:
			fh.write(alignStr)


## mainly for test for alignment feature generation
if __name__ == '__main__':

	if len(sys.argv) < 4:
		print 'python AlignmentUtils.py alnFile tgtFile tplFolder'
		print '	this entrance is mainly used for test some feature generation methods'
		print '	alnFile: a pairwise alignment file in FASTA format'
		print '	tgtFile: a tgt file'
		print '	tplFolder: the folder for tpl.pkl files, in which the template in alnFile shall have a tpl.pkl file'
		exit(1)

	alnfile = sys.argv[1]
	tgtfile = sys.argv[2]
	tplfolder = sys.argv[3]

	tgt = LoadTGT(tgtfile)

	#queryName = sys.argv[1]
	#templateName = sys.argv[2]

	#aliDir = '/mnt/local/jinbo/CathData/Align4CathS35ToPDB70/'
	#tgtDir = '/mnt/local/jinbo/CathData/SampledMSAs/MSA_2017_E001_S0/'
	#tplDir = '/mnt/local/jinbo/TPL-PDB70_13June18/'
	#pdbDir = '/mnt/local/jinbo/TPLCoordinates-PDB70_13June18/'	

	#feature = GenerateAlignmentFeatures( (queryName, templateName), aliDir=aliDir, tgtDir=tgtDir, tplDir=tplDir, pdbDir=pdbDir)
	#feature = GenerateAlignmentFeatures2( (queryName, templateName), aliDir=aliDir, tplDir=tplDir)
	feature = GenerateAlignmentFeatures4(tgt, aliFile=alnfile, tplFolder=tplfolder, queryName=tgt['name'])
	print feature
