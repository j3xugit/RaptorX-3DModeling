#!/usr/bin/env python

"""
extract alignments from the hhr result files created with hhblits|hhsearch|hhalign -o <hhr_file>
Input: .hhr file, query.hhm E-value cutoff, seqID cutoff
this script is adapted from one file in the HHsuite package
"""

import os
import sys
import numpy as np
import getopt
import random

from collections import namedtuple

from ExtractSeqFromHHM import ExtractSeqFromHHM


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


def WriteOneAlignment(query, template, query_ali, template_ali, savefolder, templateFirst=True):
	if templateFirst:
		savefile = template + '-' + query + '.fasta'
		content = [ '>' + template, template_ali, '>' + query, query_ali ]
	else:
		savefile = query + '-' + template + '.fasta'
		content = ['>' + query, query_ali, '>' + template, template_ali ]

	savefile = os.path.join(savefolder, savefile)	
	with open(savefile, 'w') as fh:
		fh.writelines( '\n'.join(content) )

def WriteOneAlignmentInResult(result, queryName, querySeq, savefolder):
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
	WriteOneAlignment(query_name, template_name, query_ali, template_ali, savefolder)

def ReadClusterFile(PDBClusterFile):
	with open(PDBClusterFile, 'r') as fh:
		content = [ line.strip() for line in list(fh) ]

	clusters = dict()
	for row, i  in zip(content, xrange(len(content)) ):
		names = row.split()
		for name in names:
			clusters[name] = i

	return clusters
	

def ReadTemplateListFile(templateFile):
	with open(templateFile, 'r') as fh:
		templates = [ line.strip() for line in list(fh) ]
	return templates	

def TemplateAllowed(templateName, allowedTemplates):
	if allowedTemplates is None:
		return True
	return templateName in allowedTemplates

def SelectTemplate(templateName, selectedTemplates, selectedClusters, clusters, allowedTemplates):
	if not TemplateAllowed(templateName, allowedTemplates):
		return False

	if clusters is None:
		return True

	if clusters.has_key(templateName):
		if clusters[templateName] not in selectedClusters:
			return True
		else:
			return False
	else:
		if templateName not in selectedTemplates:
			return True
		else:
			return False


## filter alignments, only keep those alignments with E<cutoff, seqID<MaxSeqID
## no two chosen templates are in the same cluster 
def FilterAlignments(hhrfile, query, savefolder, PDBClusterFile=None, Evalues=[0.00001, 0.001, 0.1], MaxSeqID=0.8, templateFile=None, bWriteAlignment=True, maxNumTemplates=15): 
	results = read_result(hhrfile)
	if not bool(results):
		return None

	simpleHHM = ExtractSeqFromHHM(query)
	querySeq = simpleHHM['sequence']
	queryName = simpleHHM['name']

	selectedTemplates = set()
	selectedClusters= set()

	if PDBClusterFile is not None:
		clusters = ReadClusterFile(PDBClusterFile)
	else:
		clusters = None

	if templateFile is not None:
		allowedTemplates = set(ReadTemplateListFile(templateFile) )
	else:
		allowedTemplates = None

	resultPool = []
	for E, eindex in zip(Evalues, range(len(Evalues)) ):
		for result, i in zip(results, range(len(results)) ):
			#print result.evalue, result.identity
			if result.evalue > E or ( result.evalue < 0.00001 and result.identity > MaxSeqID):
				continue

			selected=SelectTemplate(result.template_id, selectedTemplates, selectedClusters, clusters, allowedTemplates)
			if not selected:
				continue

			selectedTemplates.add(result.template_id)
			if clusters is not None and clusters.has_key(result.template_id):
				selectedClusters.add( clusters[result.template_id] )

			resultPool.append(result)

			if len(resultPool) > maxNumTemplates:
				break

			## we do not want too many templates with a large E values
			if eindex > 0 and len(resultPool) > 2:
				break

		## if got one, then done; otherwise use a larger E-value as threshold
		if len(resultPool) > 0:
			break

	## if no templates are selected, then return None
	if len(resultPool)<1 or len(selectedTemplates)<1:
		return None

	if bWriteAlignment:
		for result in resultPool:
			WriteOneAlignmentInResult(result, queryName, querySeq, savefolder)

	## write one line for the group file
	content = [ queryName ]
	for template in selectedTemplates:
		content.append( queryName + '-' + template )

	return ' '.join(content)


def Usage():
	print 'python ExtractAlignmentFromHHR.py hhrfile query.hhm [-E Evalue | -I MaxSeqID | -c PDBClusterFile | -t templateListFile | -d savefolder | -n maxNumTemplatesPerQuery ]'
	print ' This script extracts pairwise alignments from a .hhr file, mainly for training'
	print '	-E: a list of Evalues (from small to large, separated by ,) for an alignment to be extracted, default 0.00001,0.001,0.1'
	print '	-I: the max SeqID (0-1) between template and query allowed for an alignment to be extracted, default 1., i.e., output alignments of high seq id'
	print '	-c: this file is used to filter templates so that no two are in the same cluster.'
	print '	     Usually bc-30.out or bc-40.out downloaded from PDB is used. When not provided, no templates will be filtered'
	print '	-t: if provided, only select templates in this file'
	print '	-d: the folder for saving all pairwise alignments, default current work directory. Recommended since lots of alignments will be generated'
	print '	-n: the maximum number of selected templates, default 15'
	print '	'
	print '	Each pairwise alignment is written to an individual file in FASTA format with template being placed before query'
	print '	     In an alignment file, the query sequence shall be complete, but the template sequence may be just a substring of the full template sequence'

if __name__ == "__main__":

	if len(sys.argv) < 3:
		Usage()
		exit(1)

	hhrfile = sys.argv[1]
	if not os.path.isfile(hhrfile):
		print 'ERROR: invalid .hhr file:', hhrfile
		exit(1)

	queryHHM = sys.argv[2]
	if not os.path.isfile(queryHHM):
		print 'ERROR: invalid .hhm file:', queryHHM
		exit(1)

	try:
		opts, args = getopt.getopt(sys.argv[3:],"E:I:c:t:d:n:",["Evalue=", "MaxSeqID=", "PDBClusterFile=", "templateList=", "savefolder=", "maxNumTemplates=" ])
                #print opts, args
        except getopt.GetoptError:
                Usage()
                exit(1)

	Evalues = [0.00001, 0.001, 0.1]
	MaxSeqID = 1.0
	clusterFile=None
	templateFile=None
	maxNumTemplates = 15

	savefolder = os.getcwd()

        for opt, arg in opts:
                if opt in ("-E", "--Evalue"):
			fields = arg.split(',')
                        Evalues = [ np.float32(e) for e in fields ]

                elif opt in ("-I", "--MaxSeqID"):
                        MaxSeqID = np.float32(arg)

                elif opt in ("-c", "--PDBClusterFile"):
                        clusterFile = arg
			if not os.path.isfile(clusterFile):
				print 'ERROR: invalid PDB clustering file', clusterFile
				exit(1)

                elif opt in ("-t", "--templateList"):
                        templateFile = arg
			if not os.path.isfile(templateFile):
				print 'ERROR: invalid template file', templateFile
				exit(1)

                elif opt in ("-d", "--savefolder"):
                        savefolder = arg
			if not os.path.isdir(savefolder):
				os.mkdir(savefolder)

		elif opt in ("-n", "--maxNumTemplates"):
			maxNumTemplates = np.int32(arg)

                else:
                        Usage()
                        exit(1)


	res = FilterAlignments(hhrfile, queryHHM, savefolder, PDBClusterFile=clusterFile, Evalues=Evalues, MaxSeqID=MaxSeqID, templateFile=templateFile, maxNumTemplates=maxNumTemplates)
	if res is not None:
		print res

