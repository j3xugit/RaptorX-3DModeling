'''
fin = open('/mnt/data/CASP11/QualityAssess/data/17w/17w.fasta', 'r')
alignment_result = fin.readlines()
fin.close()

for i in range(len(alignment_result)/4):
	tpl_name = alignment_result[i*4].strip()[1:]
	tpl_seq = alignment_result[i*4 + 1].strip()
	tgt_name = alignment_result[i*4 + 2].strip()[1:]
	tgt_seq = alignment_result[i*4 + 3].strip()

	if tpl_seq[0] == '-':
		print "line number = %d" % (i*4 + 1)
		print tpl_name
		print tpl_seq
		print tgt_name
		print tgt_seq

	if not (len(tgt_seq) == len(tpl_seq)):
		print 'not equal'

'''

import cPickle

fin = open('pdb25_all_feature_5_from_20.pkl', 'r')
rank_result = cPickle.load(fin)
fin.close()

print len(rank_result)

