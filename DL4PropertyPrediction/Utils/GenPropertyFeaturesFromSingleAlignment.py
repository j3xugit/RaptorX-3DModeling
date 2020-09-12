import cPickle
import sys
import os
import numpy as np

from PropertyUtils import SS8Coding, SS3Coding
from Alignment import SimilarityScore

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
#def GenPropertyFeaturesFromSingleAlignment(alignment, tpl, tgt, tplAngleInfo):
def GenPropertyFeaturesFromSingleAlignment(alignment, tpl, tgt):
	##check consistency between alignment, tpl and tgt
	if len(alignment[0]) != len(alignment[1]):
		print 'ERROR: inconsistent length of query and template in the below alignment:'
		print alignment
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

	##here we require that the query sequence in alignment is exactly same as that in tgt
	assert (tgt_start == 0)


	##index for tgt and tpl, respectively
	tgt_pos = tgt_start
	tpl_pos = tpl_start

	## the below binary vector indicates insertions in query sequence, indicated by a flag 1
	insertX = np.ones( (tgt['length'], 1), dtype=np.int32)

	## it is set to true if the aligned template position has missing coordinates
	missingX = np.zeros( (tgt['length'], 1), dtype=np.int32)

	##for each query position, there is a vector of features
	numFeatures = 10
	sequentialFeatures = np.zeros( (tgt['length'], numFeatures), dtype=np.float32)

	tplPropertySize = 16
	effTemplateProperties = np.zeros( (tgt['length'], tplPropertySize), dtype=np.float32)

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

		"""
		## here we ignore all structure-related scores
		#//secondary structure score
		SS3 = SimilarityScore.SSOf2Pos1(x, y, tpl, tgt)
		SS8 = SimilarityScore.SS8Of2Pos(x, y, tpl, tgt)
		#//acc
		ACC = SimilarityScore.ACC_New_Score2(x, y, tpl, tgt)
		envScore = SimilarityScore.SingleOf2Pos_old(x, y, tpl, tgt)
		wsEnvScore = SimilarityScore.SingleOf2Pos_WS(x, y, tpl, tgt)
		"""

		point_feature = [seq_Id, blosum80, blosum62, blosum45, spScore, spScore_ST, pmScore, pmScore_ST, cc, hdsm]
			
		sequentialFeatures[tgt_pos]=np.array(point_feature).astype(np.float32)

		## template properties
		tplProperties = []

		## if one template position has no 3D coordinates
		if tpl['missing'][tpl_pos]:
			missingX[tgt_pos] = 1

			##for 8-state secondary structure, teat it as loop
			tplProperties.extend( list(SS8Coding['L'] ) )

			##for pACC, treat it as exposed
			tplProperties.append(1.)

			## for CNa, CNb, Phi, Psi, Theta, Tau
			tplProperties.extend( [0.0]* (tplPropertySize - 9 -1 ) )

			#for Omg
			tplProperties.append(np.pi)
		else:

			tplProperties.extend( list(SS8Coding[tpl['SS_str'][tpl_pos] ]) )
			tplProperties.append(tpl['pACC'][tpl_pos]/100.)
			tplProperties.append(tpl['CNa'][tpl_pos]/10.)
			tplProperties.append(tpl['CNb'][tpl_pos]/10.)

			tplAngleInfo = tpl
			tplProperties.append(tplAngleInfo['Phi'][tpl_pos])
			tplProperties.append(tplAngleInfo['Psi'][tpl_pos])
			tplProperties.append(tplAngleInfo['Theta'][tpl_pos])
			tplProperties.append(tplAngleInfo['Tau'][tpl_pos])
		
			if tplAngleInfo.has_key('Omg'):	
				tplProperties.append(tplAngleInfo['Omg'][tpl_pos])
			else:
				tplProperties.append(np.pi)

		effTemplateProperties[tgt_pos] = np.array(tplProperties).astype(np.float32)
			
		tpl_pos += 1
		tgt_pos += 1

	return np.hstack([insertX, sequentialFeatures]), XYresidues, effTemplateProperties

"""
if __name__ == '__main__':
	fin = open('17w_tpl_tgt_alignment.pkl', 'rb')
	rank_result, tpl_dic, tgt_dic  = cPickle.load(fin)
	fin.close()
	print len(rank_result)
	#template_name, query_name, tpl_start, query_start, tpl_seq, query_seq = rank_result[0][0]

	print len(tpl_dic)
	print len(tgt_dic)	

	for ind in range(len(rank_result)):
		template_name,  tpl_seq, query_name,query_seq = rank_result[ind]

		tpl = tpl_dic[template_name]
		tgt = tgt_dic[query_name]

		print "tpl name"
		print template_name
		print "tpl seq in tpl file"
		print tpl['sequence']
		print "alignment in tpl"
		print tpl_seq
		print "tpl SS"
		print tpl['SS_str']
		print "tpl ACC"
		print tpl['ACC']
		print "tpl SS3"
		print tpl['SS3']
		print "tpl SS8"
		print tpl['SS8']

		print "tgt name"
		print query_name
		print "tgt seq"
		print tgt['sequence']
		print "alignment in tgt"
		print query_seq
		print "tpl ACC"
		print tgt["ACC"]
		print "tgt ss3"
		print tgt['SS3']
		print "tgt ss8"
		print tgt['SS8']
		
		result = sim_feature_for_17w(template_name, query_name, tpl_seq, query_seq)
		
		for feature_pos in result:
			print feature_pos

"""
