import sys
import os
import numpy as np
import cPickle

import SimilarityScore
from LoadAlignments import LoadAlignments

def Usage():
	print "python ScoreAlignments.py alignmentFile tplDir tgtDir [seqBounds]"
	print "	This script scores one or multiple alignments in a file (FASTA format)"
	print "	tplDir: the folder for tpl or tpl.pkl files"
	print "	tgtDir: the folder for tgt or tgt.pkl files"
	print " seqBounds: the sequence segment for which you would like to calculate alignment score" 
	print "		e.g., 1-100, both ends are inclusive. Residue index starts from 1 to sequence length"
	

##remove gaps in the alignment and return the real template and query sequences
def ExtractSeqFromAlignment(alignment):
	template = alignment[0].translate(None, '-')
	query = alignment[1].translate(None, '-')
	return (template, query)
	
## alignment is a tuple, the first entry is the template sequence in alignment including gaps
## the second entry is the query sequence in alignment including gaps
## tpl is the template object and tgt is the query object

##this function generates features for those query positions with effective aligned template positions
##do we need to deal with HIS tags?
## when seqBounds is not None, then only calculate the alignment score for the target segment defined by seqBounds
def ScoreOneAlignment(alignment, tpl, tgt, bounds=None):
	seqBounds = None
	if bounds is not None:
		assert bounds[0]>0
		if len(bounds) == 1:
			seqBounds = (bounds[0], tgt['length'])
		elif len(bounds) == 2:
			seqBounds = bounds
		else:
			print 'ERROR: incorrect information for seq bounds: ', bounds
			exit(1)
	
		assert seqBounds[1]> seqBounds[0]
		assert seqBounds[1]<= tgt['length']

	##check consistency between alignment, tpl and tgt
	if len(alignment[0]) != len(alignment[1]):
		print 'ERROR: the length of query and template in alignment is inconsistent'
		exit(1)

	template, query = ExtractSeqFromAlignment(alignment)

	## template and query shall be the substring of tpl and tgt, respectively
	tpl_start = tpl['sequence'].find(template)
	if tpl_start == -1:
		print 'ERROR: the template sequence in alignment is not a substring of the sequence in tpl', tpl['name']
		exit(1)

	tgt_start = tgt['sequence'].find(query)
	if tgt_start == -1:
		print 'ERROR: the query sequence in alignment is not a substring of the sequence in tgt', tgt['name']
		exit(1)

	""" 
	## if you want exact match, then use the below code
	if template != tpl['sequence']:
		print 'the template sequence in alignment is inconsistent with that in tpl'
		exit(1)
	if query != tgt['sequence']:
		print 'the query sequence in alignment is inconsistent with that in tgt'
		exit(1)
	"""


	##index for tgt and tpl, respectively
	tgt_pos = tgt_start
	tpl_pos = tpl_start

	## the below binary vector indicates insertions in query sequence, indicated by a flag 1
	insertX = np.ones( (tgt['length'], 1), dtype=np.int32)

	## it is set to true if the aligned template position has missing coordinates
	missingX = np.zeros( (tgt['length'], 1), dtype=np.int32)

	##for each query position, there is a vector of features
	numFeatures = 11
	localScores = np.zeros( (tgt['length'], numFeatures), dtype=np.float32)

	##for debug only
	XYresidues = - np.ones( (tgt['length'], 4), dtype=np.int16)

	for al_pos in range(len(alignment[0])):

		## there is a gap in template, i.e., an insertion in query
		if alignment[0][al_pos] == '-':
			## need to generate some flag features for insertion in query
			insertX[tgt_pos] = 1

			##write down the information at query sequence for debug only
			XYresidues[tgt_pos][0]= tgt_pos
			XYresidues[tgt_pos][1] = ord(tgt['sequence'][tgt_pos]) - ord('A')

			tgt_pos += 1
			continue


		## if there is a gap in query, just skip it (shall we change this strategy later?)
		if alignment[1][al_pos] == '-':
			tpl_pos += 1
			## no need to generate flag features for insertion in template
			continue

		## for match state
		insertX[tgt_pos] = 0

		##write down the information at both query and template for debug only
		XYresidues[tgt_pos][0]= tgt_pos
		XYresidues[tgt_pos][1] = ord(tgt['sequence'][tgt_pos]) - ord('A')
		XYresidues[tgt_pos][2]= tpl_pos
		XYresidues[tgt_pos][3] = ord(tpl['sequence'][tpl_pos]) - ord('A')

		## sAA and tAA the amino acids at query and template, respectively
		tAA = tpl['sequence'][tpl_pos]
		sAA = tgt['sequence'][tgt_pos]
		
		seq_Id = int(tAA == sAA)
		
		blosum80 = SimilarityScore.BLOSUM80(tAA, sAA)
		blosum62 = SimilarityScore.BLOSUM62(tAA, sAA)
		blosum45 = SimilarityScore.BLOSUM45(tAA, sAA)

		cc = SimilarityScore.CC50(tAA, sAA)      
		hdsm = SimilarityScore.HDSM(tAA, sAA)

		x, y = tpl_pos, tgt_pos
		spScore = SimilarityScore.MutationOf2Pos6(x, y, tpl, tgt)
		spScore_ST = SimilarityScore.MutationOf2Pos6_ST(x, y, tpl, tgt)
		pmScore = SimilarityScore.MutationOf2Pos5(x, y, tpl, tgt)
		pmScore_ST = SimilarityScore.MutationOf2Pos5_ST(x, y, tpl, tgt)

		SS3 = SimilarityScore.SSMutationScore_3State(x, y, tpl, tgt)

		"""
		## here we ignore the other structure-related scores
		#//secondary structure score
		SS8 = SimilarityScore.SS8Of2Pos(x, y, tpl, tgt)
		#//acc
		ACC = SimilarityScore.ACC_New_Score2(x, y, tpl, tgt)
		envScore = SimilarityScore.SingleOf2Pos_old(x, y, tpl, tgt)
		wsEnvScore = SimilarityScore.SingleOf2Pos_WS(x, y, tpl, tgt)
		"""

		localScore = [seq_Id, blosum80, blosum62, blosum45, spScore, spScore_ST, pmScore, pmScore_ST, cc, hdsm, SS3]
			
		localScores[tgt_pos]=np.array(localScore).astype(np.float32)

		## if one template position has no 3D coordinates
		if tpl['missing'][tpl_pos]:
			missingX[tgt_pos] = 1

		tpl_pos += 1
		tgt_pos += 1


	##calculate globalScore
	if seqBounds is not None:
		globalScore = np.sum(np.multiply(localScores, 1-insertX)[seqBounds[0]-1:seqBounds[1] ], axis=0)
	else:
		globalScore = np.sum(np.multiply(localScores, 1-insertX), axis=0)

	##calculate the number of gap openings and gaps
	if seqBounds is not None:
		start = seqBounds[0]-1
		end = seqBounds[1]
	else:
		start = 0
		end = len(insertX)

	while start<end and insertX[start]:
		start += 1

	while end>start and insertX[end-1]:
		end -= 1

	if start >= end:
		numGaps = 0
		numGapOpenings = 0
	else:
		numGaps = np.sum(insertX[start: end] )
		gapOpens = [ (insertX[i]==0 and insertX[i+1]==1) for i in range(start, end-1) ]
		numGapOpenings = np.sum(gapOpens)

	if seqBounds is not None:
		numAligned = np.sum( 1 - insertX[seqBounds[0]-1:seqBounds[1] ] )
		seqLen = seqBounds[1] - seqBounds[0] + 1
	else:
		numAligned = np.sum( 1 - insertX )
		seqLen = tgt['length']

	globalScore = [seqLen, numAligned] + list(globalScore) + [numGapOpenings, numGaps]

	if seqBounds is not None:
		start = seqBounds[0]-1
		end = seqBounds[1]
		return globalScore, localScores[start:end], insertX[start:end], missingX[start:end], XYresidues[start:end]

	return globalScore, localScores, insertX, missingX, XYresidues

## score multiple alignments. In alignmentFile, for each alignment, template is placed before query sequene
def ScoreAlignments(alignmentFile, tplDir, tgtDir, seqBounds=None):
	alignments, tpls, tgts = LoadAlignments(alignmentFile, tplDir, tgtDir)

	scores = []
	for alignment in alignments:
		tplSeq, tgtSeq, tplName, tgtName = alignment
		tpl = tpls[tplName]
		tgt = tgts[tgtName]

		if seqBounds is not None:
			globalScore, localScore, insertX, missingX, _ = ScoreOneAlignment(alignment, tpl, tgt, seqBounds)
		else:
			globalScore, localScore, insertX, missingX, _ = ScoreOneAlignment(alignment, tpl, tgt)
		scores.append( (tplName, tgtName, globalScore, localScore, insertX, missingX) )

	return scores


def str_display(ls):
        if not isinstance(ls, (list, tuple, np.ndarray)):
                str_ls = '{0:.2f}'.format(ls)
                return str_ls

        str_ls = ['{0:.2f}'.format(v) for v in ls ]
        str_ls2 = '\t'.join(str_ls)
        return str_ls2

def main(argv):
	if len(argv) < 3 or len(argv) > 4:
		Usage()
		exit(1)

	alignFile = argv[0]
	tplDir = argv[1]
	tgtDir = argv[2]

	seqBounds = None

	if len(argv)==4:
		seqBounds = []
		fields = argv[3].split('-')
		assert len(fields) >= 1
		seqBounds.append( np.int32(fields[0]) )
		if len(fields)==2:
			seqBounds.append(np.int32(fields[1]) )

	if seqBounds is not None:
		scores = ScoreAlignments(alignFile, tplDir, tgtDir, tuple(seqBounds) )
	else:
		scores = ScoreAlignments(alignFile, tplDir, tgtDir)

	scores.sort(key=lambda x: x[2][6], reverse=True)

	for s in scores:
		print "%s\t%s\t%s" % (s[0], s[1], str_display(s[2]) )

if __name__ == '__main__':
	main(sys.argv[1:])
