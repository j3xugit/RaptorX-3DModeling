import cPickle

import load_tpl_tgt
import load_rank

TGTDir = '/mnt/data/RaptorXCommon/TrainTestData/Distance_Contact_TrainData/7952_Old_PDB25_Training_Data/training_tgt'
TPLDir = '/mnt/data/RaptorXCommon/Server_Database/TPL_BC100'
RANKDir = '/mnt/data/RaptorXCommon/TrainTestData/Ranking_Project/raw_data/7952_against_bc40_rank'

def load_protein_list(list_file):
	fin = open(list_file, 'r')
	lines = [line.strip() for line in fin.readlines()]
	fin.close()
	return lines

if __name__ == '__main__':

	
	tgt_list = load_protein_list('pdb25_protein_list')
	tpl_list = load_protein_list('TPL_BC100')
	rank_list = load_protein_list('rank_file_list')
	tgt_dir = '/mnt/data/RaptorXCommon/TrainTestData/Distance_Contact_TrainData/7952_Old_PDB25_Training_Data/training_tgt'
	tpl_dir = '/mnt/data/RaptorXCommon/Server_Database/TPL_BC100'
	rank_dir = '/mnt/data/RaptorXCommon/TrainTestData/Ranking_Project/raw_data/7952_against_bc40_rank'
	
	'''
	tgt_list = load_protein_list('17w_tgt_list')
	tpl_list = load_protein_list('17w_tpl_list')
	#rank_list = load_protein_list('rank_file_list')

	tgt_dir = '/mnt/data/CASP11/QualityAssess/data/17w/17w_tgt'
	tpl_dir = '/mnt/data/CASP11/QualityAssess/data/17w/17w_tpl'
	alignment_file = '/mnt/data/CASP11/QualityAssess/data/17w/17w.fasta'
	'''

#	rank_dir = 'rank_3000'
	
	tpl_dic = {}
	tgt_dic = {}
	rank_result = []

	'''	
	alignment = []
	
	fin = open(alignment_file, 'r')
	alignment_result = fin.readlines()
	fin.close()
	
	for i in range(len(alignment_result)/4):
		tpl_name = alignment_result[i*4].strip()[1:]
		tpl_seq = alignment_result[i*4+1].strip()
		tgt_name = alignment_result[i*4+2].strip()[1:]
		tgt_seq = alignment_result[i*4+3].strip()
		
		alignment.append((tpl_name, tpl_seq, tgt_name, tgt_seq))
	
		
	'''

	count = 0
	for protein in rank_list:
		rank_result.append(load_rank.load_rank(rank_dir + '/' + protein + '.rank'))
		count += 1
		if count % 100 == 0:
			print 'processing rank file ' + str(count)
	
	
	count = 0	
	for protein in tgt_list:
#		tpl_dic[protein] = load_tpl_tgt.load_tpl(tpl_tgt_dir + '/training_tpl/' + protein + '.tpl')
		tgt_dic[protein] = load_tpl_tgt.load_tgt(tgt_dir + '/' + protein + '.tgt')
#		dis_dic[protein] = load_tpl_tgt.load_dis(tpl_tgt_dir + '/training_distcb/' + protein + '.distcb')
		count += 1
		if count % 100 == 0:
			print 'processing tgt file ' + str(count)
	count = 0	
	for protein in tpl_list:
#		print protein
		tpl_dic[protein] = load_tpl_tgt.load_tpl(tpl_dir + '/' + protein + '.tpl')
#		tgt_dic[protein] = load_tpl_tgt.load_tgt(tpl_tgt_dir + '/training_tgt/' + protein + '.tgt')
#		dis_dic[protein] = load_tpl_tgt.load_dis(tpl_tgt_dir + '/training_distcb/' + protein + '.distcb')
		count += 1
		if count % 100 == 0:
			print 'processing tpl file ' + str(count)
	
	resultfile = 'pdb25_tpl_tgt_alignment_top40.pkl'
	resultfh = open(resultfile,'wb')
	cPickle.dump( (rank_result, tpl_dic, tgt_dic), resultfh, protocol=cPickle.HIGHEST_PROTOCOL)
	resultfh.close()
