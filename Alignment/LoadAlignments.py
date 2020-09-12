import numpy as np
import os
import sys

from Common import LoadTPLTGT

## this function loads a file containing some alignments. Each alignment has two proteins. The first one is assumed to be a template and the 2nd one is a query sequence.
## in addition to loading alignments, this function also loads the related template file (*.tpl) and query sequence file (*.tgt)

## alignment_file: a text file contains a set of alignments. The number of lines in this file shall be multiple of 4 since each alignment has 4 lines.
## one alignment example is as follows
## >XXXXA
## AAAAAAAAAAAAAACCCCCCCDDDDDDDDDDDDGGGGGGGGGGGHH--HHHHHHHHH
## >query
## AA--AAAAAAAAACCCCC-CCDD-DDDD-DDDGGGGGGGGGGGGHHHHHHHH-----

## tpl_dir: the folder containing the template files related to the alignment file
## tgt_dir: the folder containing the query sequence files related to the alignment file

## this function returns 3 objects: a list of alignments, a dictionary of templates, and a dictionary of query sequences
## each alignment is a tuple of 4 elements: template_sequence in the alignment (including gaps), query_sequence in the alignment (including gaps), template name, query name
## if LoadTGTTPL is False, then return only the template name set and the query name set instead of their contents

## KeepTPLSimple is valid only when LoadTGTTPLcontent is True. When it is True, do not store distance/orientation information for a template to save memory
## Otherwise, store the distance/orientation information of a template.
def LoadAlignments(alignment_file, tpl_dir, tgt_dir, LoadTGTTPLcontent=True, KeepTPLSimple=False):

        alignments = []

        fin = open(alignment_file, 'r')
        content = [ line.strip() for line in list(fin) ]
        fin.close()

        ##remove empty lines
        alignment_result = [ c for c in content if c ]

        if len(alignment_result)%4 !=0 :
                print 'ERROR: the number of lines in the alignment file is incorrect: ', alignment_file
                exit(1)

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
                	tpl_file = os.path.join(tpl_dir, tpl_name+'.tpl.pkl')
			if not os.path.isfile(tpl_file):
                        	print 'WARNING: invalid tpl file %s' % tpl_file
                        	continue

                ##check to see if the tgt file exists or not
                tgt_file = os.path.join(tgt_dir, tgt_name+'.tgt')
                if not os.path.isfile(tgt_file):
                	tgt_file = os.path.join(tgt_dir, tgt_name + '.tgt.pkl')
                	if not os.path.isfile(tgt_file):
                        	print 'WARNING: invalid tgt file %s' % tgt_file
                        	continue

                alignments.append((tpl_seq, tgt_seq, tpl_name, tgt_name))

                templates.add(tpl_name)
                querys.add(tgt_name)

        print 'In total loaded ', len(alignments), ' alignments involving ', len(templates), ' templates and ', len(querys), ' query sequences'
	if not LoadTGTTPLcontent:
		return alignments, templates, querys

        ## load tgt and tpl
        tplPool = {}
        for tplName in templates:
                tpl_file = os.path.join(tpl_dir, tplName + '.tpl')
                if not os.path.isfile(tpl_file):
                	tpl_file = os.path.join(tpl_dir, tplName + '.tpl.pkl')
			if not os.path.isfile(tpl_file):
                        	print 'ERROR: invalid tpl file %s' % tpl_file
                        	exit(1)

                tpl = LoadTPLTGT.load_tpl(tpl_file)
		if KeepTPLSimple:
			tpl.pop('atomDistMatrix', None)
			tpl.pop('atomOrientationMatrix', None)
		tplPool[tplName] = tpl
				

        tgtPool = {}
        for tgtName in querys:
                tgt_file = os.path.join(tgt_dir, tgtName + '.tgt')
                if not os.path.isfile(tgt_file):
                	tgt_file = os.path.join(tgt_dir, tgtName + '.tgt.pkl')
			if not os.path.isfile(tgt_file):
                        	print 'ERROR: invalid tgt file %s' % tgt_file
                        	exit(1)

                tgtPool[tgtName] = LoadTPLTGT.load_tgt(tgt_file)

        return alignments, tplPool, tgtPool
