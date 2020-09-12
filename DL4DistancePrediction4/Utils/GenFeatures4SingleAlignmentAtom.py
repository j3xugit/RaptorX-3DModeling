import cPickle
import sys
import numpy as np

sys.path.append('../../Alignment/Scripts')
import SimilarityScore

"""
This scripts generate 17 features describing a query-template alignment
The query sequence in the alignment shall be exactly same as that in the tgt file
the template sequence in the alignment shall be a subsequence of SEQRES in the tpl file  
"""

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
def GenFeature4Alignment(alignment, tpl, tgt, tplAtomCoordinates):
	

	print 'Generating similarity features for template ', tpl['name'], ' and query ', tgt['name']

	##check consistency between alignment, tpl and tgt
	if len(alignment[0]) != len(alignment[1]):
		print 'the length of query and template in alignment is inconsistent'
		exit(-1)

	template, query = ExtractSeqFromAlignment(alignment)

	## template and query shall be the substring of tpl and tgt, respectively

	tpl_start = tpl['sequence'].find(template)
	if tpl_start == -1:
		print 'the template sequence in alignment is not a substring of the sequence in tpl'
		exit(-1)

	tgt_start = tgt['sequence'].find(query)
	if tgt_start == -1:
		print 'the query sequence in alignment is not a substring of the sequence in tgt'
		exit(-1)

	##wrong if tgt_start is not 0, here we require that the query sequence in alignment is exactly same as that in tgt
	assert (tgt_start == 0)


	##index for tgt and tpl, respectively
	tgt_pos = tgt_start
	tpl_pos = tpl_start

	## insertion in query, indicated by a flag 1
	insertX = np.ones( (tgt['length'], 1), dtype=np.uint8)

	## it is set to true if the aligned template position has missing corrdinates
	missingX = np.zeros( (tgt['length'], 1), dtype=np.uint8)

	##effective coordinates copied from aligned template positions
	#effTemplateCoordinates = np.ones( (tgt['length'], 3), dtype=np.float32) * 999.0
	effTemplateCoordinates = [ dict() ] * tgt['length']

	##for each query position, there is a vector of features
	sequentialFeatures = np.zeros( (tgt['length'], 15), dtype=np.float32)

	##for debug
	XYresidues = - np.ones( (tgt['length'], 4), dtype=np.int16)

	for al_pos in range(len(alignment[0])):

		## there is a gap in template, i.e., an insertion in query
		if alignment[0][al_pos] == '-':
			## need to generate some flag features for insertion in query
			insertX[tgt_pos] = 1

			XYresidues[tgt_pos][0]= tgt_pos
			XYresidues[tgt_pos][1] = ord(tgt['sequence'][tgt_pos]) - ord('A')

			tgt_pos += 1
			continue


		## if there is a gap in query, just skip it
		if alignment[1][al_pos] == '-':
			tpl_pos += 1
			## no need to generate flag features for insertion in template
			continue

		insertX[tgt_pos] = 0

		XYresidues[tgt_pos][0]= tgt_pos
		XYresidues[tgt_pos][1] = ord(tgt['sequence'][tgt_pos]) - ord('A')
		XYresidues[tgt_pos][2]= tpl_pos
		XYresidues[tgt_pos][3] = ord(tpl['sequence'][tpl_pos]) - ord('A')

		tAA = tpl['sequence'][tpl_pos]
		sAA = tgt['sequence'][tgt_pos]
		
		seq_Id = int(tAA == sAA)
		
		blosum80 = SimilarityScore.BLOSUM80(tAA, sAA)
		blosum62 = SimilarityScore.BLOSUM62(tAA, sAA)
		blosum45 = SimilarityScore.BLOSUM45(tAA, sAA)

		x, y = tpl_pos, tgt_pos
		spScore = SimilarityScore.MutationOf2Pos6(x, y, tpl, tgt)
		spScore_ST = SimilarityScore.MutationOf2Pos6_ST(x, y, tpl, tgt)
		pmScore = SimilarityScore.MutationOf2Pos5(x, y, tpl, tgt)
		pmScore_ST = SimilarityScore.MutationOf2Pos5_ST(x, y, tpl, tgt)

		#//secondary structure score
		SS3 = SimilarityScore.SSMutationScore_3State(x, y, tpl, tgt)
		SS8 = SimilarityScore.SSMutationScore_6State(x, y, tpl, tgt)
		#//acc
		ACC = SimilarityScore.ACCMutationScore_3State(x, y, tpl, tgt)
		envScore = SimilarityScore.SingletonScore_ProfileBased(x, y, tpl, tgt)
		wsEnvScore = SimilarityScore.SingletonScore_WS(x, y, tpl, tgt)

		cc = SimilarityScore.CC50(tAA, sAA)       #   //offset by 0.5 to make the expected value 0
		hdsm = SimilarityScore.HDSM(tAA, sAA)

		point_feature = [seq_Id, blosum80, blosum62, blosum45, spScore, spScore_ST, pmScore, pmScore_ST, cc, hdsm, SS3, SS8, ACC, envScore, wsEnvScore]
			
		sequentialFeatures[tgt_pos]=np.array(point_feature).astype(np.float32)

		## if one template position has no 3D coordinates
		if tpl['missing'][tpl_pos]:
			missingX[tgt_pos] = 1
		else:
			#effTemplateCoordinates[tgt_pos] = tpl['Cb'][tpl_pos]
			effTemplateCoordinates[tgt_pos] = tplAtomCoordinates[tpl_pos]
			
		tpl_pos += 1
		tgt_pos += 1

	##calculate the distance matrix from templates
	templateDistanceMatrix = dict()
	templateDistanceMatrix['CbCb'] = - np.ones( (tgt['length'], tgt['length']), dtype=np.float16)
	templateDistanceMatrix['CaCa'] = - np.ones( (tgt['length'], tgt['length']), dtype=np.float16)
	templateDistanceMatrix['CgCg'] = - np.ones( (tgt['length'], tgt['length']), dtype=np.float16)
	templateDistanceMatrix['CaCg'] = - np.ones( (tgt['length'], tgt['length']), dtype=np.float16)
	templateDistanceMatrix['NO'] = - np.ones( (tgt['length'], tgt['length']), dtype=np.float16)

	for i in range( tgt['length'] ):
		if insertX[i] or missingX[i]:
			continue

		templateDistanceMatrix['CbCb'][i, i]=0
		templateDistanceMatrix['CaCa'][i, i]=0
		templateDistanceMatrix['CgCg'][i, i]=0
		templateDistanceMatrix['CaCg'][i, i]=0
		templateDistanceMatrix['NO'][i, i]=0

		for j in range(i+1, tgt['length'] ):
			if insertX[j] or missingX[j]:
				continue

			if (effTemplateCoordinates[i]['Cb'] is not None) and (effTemplateCoordinates[j]['Cb']  is not None):
				templateDistanceMatrix['CbCb'][i, j] = np.linalg.norm( effTemplateCoordinates[i]['Cb'] - effTemplateCoordinates[j]['Cb'] ) 
				templateDistanceMatrix['CbCb'][j, i] = templateDistanceMatrix['CbCb'][i, j]

			if (effTemplateCoordinates[i]['Ca'] is not None) and (effTemplateCoordinates[j]['Ca']  is not None):
				templateDistanceMatrix['CaCa'][i, j] = np.linalg.norm( effTemplateCoordinates[i]['Ca'] - effTemplateCoordinates[j]['Ca'] ) 
				templateDistanceMatrix['CaCa'][j, i] = templateDistanceMatrix['CaCa'][i, j]

			if (effTemplateCoordinates[i]['Cg'] is not None) and (effTemplateCoordinates[j]['Cg']  is not None):
				templateDistanceMatrix['CgCg'][i, j] = np.linalg.norm( effTemplateCoordinates[i]['Cg'] - effTemplateCoordinates[j]['Cg'] ) 
				templateDistanceMatrix['CgCg'][j, i] = templateDistanceMatrix['CgCg'][i, j]

			if (effTemplateCoordinates[j]['Cg'] is not None) and (effTemplateCoordinates[i]['Ca'] is not None):
				templateDistanceMatrix['CaCg'][i, j] = np.linalg.norm( effTemplateCoordinates[i]['Ca'] - effTemplateCoordinates[j]['Cg'] ) 

			if (effTemplateCoordinates[i]['Cg'] is not None) and (effTemplateCoordinates[j]['Ca'] is not None):
				templateDistanceMatrix['CaCg'][j, i] = np.linalg.norm( effTemplateCoordinates[j]['Ca'] - effTemplateCoordinates[i]['Cg'] )

			if (effTemplateCoordinates[i]['N'] is not None) and (effTemplateCoordinates[j]['O']  is not None):
				templateDistanceMatrix['NO'][i, j] = np.linalg.norm( effTemplateCoordinates[i]['N'] - effTemplateCoordinates[j]['O'] ) 

			if (effTemplateCoordinates[i]['O'] is not None) and (effTemplateCoordinates[j]['N']  is not None):
				templateDistanceMatrix['NO'][j, i] = np.linalg.norm( effTemplateCoordinates[j]['N'] - effTemplateCoordinates[i]['O'] )

	"""
	print insertX.shape
	print sequentialFeatures.shape
	"""
	
	#return np.hstack([insertX, XYresidues, sequentialFeatures, effTemplateCoordinates]), templateDistanceMatrix
	return np.hstack([insertX, sequentialFeatures]), templateDistanceMatrix, XYresidues, effTemplateCoordinates

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
