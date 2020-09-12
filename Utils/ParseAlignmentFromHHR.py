#!/usr/bin/env python

"""
parse alignments from the hhr result files created with hhblits|hhsearch|hhalign -o <hhr_file>
Input: .hhr file, query.hhm and folder for the template .hhm files, E-value cutoff (default 1), topK (default 20)
this script is adapted from one file in the HHsuite package
"""

import os
import sys
import getopt

import math

import numpy as np

from collections import namedtuple

from Common.ExtractSeqFromHHM import ExtractSeqFromHHM


__author__ = 'Markus Meier (markus.meier@mpibpc.mpg.de), Jinbo Xu (jinboxu@gmail.com) '
__version__ = '1.0'
__license__ = "GPL-3"


hhr_alignment = namedtuple('hhr_alignment', ['query_id', 'query_length', 'query_neff',
                                             'template_id', 'template_length', 'template_info',
                                             'template_neff', 'query_ali', 'template_ali',
                                             'start', 'end', 'probability', 'evalue', 'score',
                                             'aligned_cols', 'identity', 'similarity', 'sum_probs'])


class HHRFormatError(Exception):
    def __init__(self, value):
        self.value = "ERROR: "+value

    def __str__(self):
        return repr(self.value)


def get_sequence_name(header):
    name = header.replace(">", "").split()[0]
    return name


def parse_result(lines):
    results = []

    query_id = None
    query_length = None
    query_neff = None
    query_seq = []
    template_id = None
    template_length = None
    template_seq = []
    template_info = None
    query_start = None
    query_end = None
    template_start = None
    template_end = None
    probability = None
    evalue = None
    score = None
    identity = None
    similarity = None
    template_neff = None
    sum_probs = None
    aligned_cols = None

    skipped_ali_tags = ["ss_dssp", "ss_pred", "Consensus"]

    is_alignment_section = False

    for line in lines:
        if(line.startswith("Query")):
            query_id = line.split()[1]
        elif(line.startswith("Match_columns")):
            query_length = int(line.split()[1])
        elif(line.startswith("Neff")):
            query_neff = float(line.split()[1])
        elif(is_alignment_section and (line.startswith("No") or line.startswith("Done!"))):
            if query_start is not None:
                result = hhr_alignment(query_id, query_length, query_neff,
                                       template_id, template_length, template_info, template_neff,
                                       query_seq, template_seq, (query_start, template_start),
                                       (query_end, template_end), probability, evalue, score,
                                       aligned_cols, identity, similarity, sum_probs)
                results.append(result)
            template_id = None
            template_info = None
            query_seq = []
            template_seq = []

            query_start = None
            query_end = None
            template_start = None
            template_end = None
        elif(line.startswith("Probab")):
            tokens = line.split()
            probability = float(tokens[0].split("=")[1])
            evalue = float(tokens[1].split("=")[1])
            score = float(tokens[2].split("=")[1])
            aligned_cols = int(tokens[3].split("=")[1])
            identity = float(tokens[4].split("=")[1].replace("%", "")) / 100.0
            similarity = float(tokens[5].split("=")[1])
            sum_probs = float(tokens[6].split("=")[1])
            if(len(tokens) > 7):
                template_neff = float(tokens[7].split("=")[1])
            continue
        elif(line.startswith(">")):
            is_alignment_section = True
            template_id = line[1:].split()[0]
            template_info = line
        elif(line.startswith("Q")):
            tokens = line.split()
            if(tokens[1] in skipped_ali_tags):
                continue

            try:
                token_2 = tokens[2].replace("(", "").replace(")", "")
                token_2 = int(token_2)
            except:
                raise HHRFormatError(("Converting failure of start index ({}) "
                                      "of query alignment").format(tokens[2]))

            if query_start is None:
                query_start = token_2
            query_start = min(query_start, token_2)

            try:
                token_4 = tokens[4].replace("(", "").replace(")", "")
                token_4 = int(token_4)
            except:
                raise HHRFormatError(("Converting failure of end index ({}) "
                                      "of query alignment").format(tokens[4]))

            if query_end is None:
                query_end = token_4
            query_end = max(query_end, token_4)
            query_seq.append(tokens[3])
        elif(line.startswith("T")):
            tokens = line.split()
            if(tokens[1] in skipped_ali_tags):
                continue
            template_seq.append(tokens[3])

            try:
                token_2 = tokens[2].replace("(", "").replace(")", "")
                token_2 = int(token_2)
            except:
                raise HHRFormatError(("Converting failure of start index ({}) "
                                      "of template alignment").format(tokens[2]))

            if template_start is None:
                template_start = token_2
            template_start = min(template_start, token_2)

            try:
                token_4 = tokens[4].replace("(", "").replace(")", "")
                token_4 = int(token_4)
            except:
                raise HHRFormatError(("Converting failure of end index ({}) "
                                      "of template alignment").format(tokens[4]))

            if template_end is None:
                template_end = token_4
            template_end = max(template_end, token_4)

            try:
                token_5 = tokens[4].replace("(", "").replace(")", "")
                token_5 = int(token_5)
            except:
                raise HHRFormatError(("Converting failure of template length ({}) "
                                      "in template alignment").format(tokens[5]))
            template_length = token_5


    if(template_id is not None and query_start is not None):
        result = hhr_alignment(query_id, query_length, query_neff,
                               template_id, template_length, template_info, template_neff,
                               "".join(query_seq), "".join(template_seq), (query_start, template_start),
                               (query_end, template_end), probability, evalue, score,
                               aligned_cols, identity, similarity, sum_probs)
        results.append(result)

    return results


def read_result(input_file):
    with open(input_file) as fh:
        lines = fh.readlines()
        return parse_result(lines)


def WriteOneAlignment(query, template, query_ali, template_ali, queryNum, templateFirst=True, savefolder=os.getcwd()):
	savefilename = query +'-HHP' + '-m' + str(queryNum).zfill(2) + '-' + template + '.fasta'
	savefile = os.path.join(savefolder, savefilename)

	if templateFirst:
		content = [ '>' + template, template_ali, '>' + query, query_ali ]
	else:
		content = ['>' + query, query_ali, '>' + template, template_ali ]
	
	with open(savefile, 'w') as fh:
		fh.writelines( '\n'.join(content) )

	return savefilename

def WriteOneAlignmentInResult(result, queryName, querySeq, queryNum, savefolder=os.getcwd() ):
	query_ali = ''.join(result.query_ali)
	template_ali = ''.join(result.template_ali)
	query_name = result.query_id
	assert query_name == queryName

	template_name = result.template_id

	## modify alignment to output the full query sequence			
	queryN = querySeq[ : result.start[0]-1 ]
	queryC = querySeq[result.end[0] : ]
	templateN = '-' * len(queryN)
	templateC = '-' * len(queryC)
	query_ali = queryN + query_ali + queryC
	template_ali = templateN + template_ali + templateC
	return WriteOneAlignment(query_name, template_name, query_ali, template_ali, queryNum, savefolder=savefolder)

def ParseAlignments(hhrfile, queryHHMFile, E=0.001, topK=20, savefolder=os.getcwd()): 
	results = read_result(hhrfile)
	if not bool(results):
		return None

	simpleHHM = ExtractSeqFromHHM(queryHHMFile)
	querySeq = simpleHHM['sequence']
	queryName = simpleHHM['name']

	firstLogE = math.log(results[0].evalue)
	resultSummary = []
	for result, i in zip(results, range(len(results)) ):
		#print result.evalue, result.identity
		if result.evalue>E or i>=topK:
			continue

		## skip if E-value is much bigger than the smallest one
		logE = math.log(result.evalue)
		if (logE-10) > firstLogE:
			continue

		savefilename = WriteOneAlignmentInResult(result, queryName, querySeq, i, savefolder)
		resultSummary.append( (savefilename, int(round(math.log(result.evalue) ) ), int(round(100*result.identity)), int(round(100*result.identity*result.aligned_cols/len(querySeq))), int(round(100*result.query_length*1./len(querySeq))) ) )

	return resultSummary


def Usage():
	print 'python ParseAlignmentFromHHR.py [-e Evalue | -k topK | -d savefolder ] hhrfile query.hhm'
	print '	This script extracts pairwise alignments from a .hhr file for template-based modeling. The .hhr file is generated by HHsearch'
	print '	-e: the max E-value for an alignment to be extracted, default 0.001'
	print '	-k: the number of top alignments to be generated, default 20'
	print '	-d: the result saving folder, default current work directory'
	print '	Each pairwise alignment is written to an individual file in FASTA format with template being placed before query'
	print '	The file is named after queryName-HHP-mXX-templateName.fasta where XX is the ranking order starting from 0'
	print '	In the alignment file, the query sequence is complete, but the template sequence may be just a substring of the original template sequence'

if __name__ == "__main__":

	if len(sys.argv) < 2:
		Usage()
		exit(1)
	try:
		opts, args = getopt.getopt(sys.argv[1:], "e:k:d:", ["Evalue=", "topK=", "savefolder=" ])
		#print opts
		#print args
	except getopt.GetoptError as err:
		Usage()
		exit(1)

	E = 0.001
	topK = 20
	savefolder = os.getcwd()

	for opt, arg in opts:
		if opt in ("-e", "--Evalue"):
			E = np.float32(arg)

		elif opt in ("-k", "--topK"):
			topK= np.int32(arg)

		elif opt in ("-d", "--savefolder"):
			savefolder = arg
			if not os.path.isdir(savefolder):
				os.mkdir(savefolder)

		else:
			Usage()
			exit(1)

	if len(args) < 2:
		Usage()
		exit(1)

	hhrfile = args[0]
	queryHHM = args[1]

	if not os.path.isfile(hhrfile) or not os.path.isfile(queryHHM):
		print >> sys.stderr, "ERROR: cannot find one of the files: ", hhrfile, queryHHM
		exit(1)

	resultSummary = ParseAlignments(hhrfile, queryHHM, E=E, topK=topK, savefolder=savefolder)

	resultSummary_sorted = sorted(resultSummary, key=lambda x:(x[3], -x[1], x[2]), reverse=True)

	summaryStrs = []
	for summary in resultSummary_sorted:
		s = '\t'.join([ str(e) for e in summary ])
		summaryStrs.append(s)

	targetName = os.path.basename(queryHHM).split('.')[0]
	summaryFile = targetName + '.summary.txt'
	savefile = os.path.join(savefolder, summaryFile)
	with open(savefile, 'w') as fh:
		fh.write('\n'.join(summaryStrs) )
