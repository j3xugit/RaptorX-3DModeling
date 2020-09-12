import os
import sys
import numpy as np

import Bio.PDB
from Bio.PDB.Polypeptide import *
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum80

def Alignment2Mapping(alignment):
	mapping = []

	numMismatches = 0
	coverage = 0
        j = -1
        for a, b in zip(alignment[0], alignment[1]):
        	if a!='-' and b!='-':
                	j += 1
                	mapping.append(j)
			if a!=b:
				numMismatches += 1
			else:
				coverage += 1
                elif a!='-':
			mapping.append(-1)
                elif b!='-':
                        j += 1
                else:
                        print 'ERROR: abnormal in the alignment'
                        exit(1)

	return mapping, coverage, numMismatches

##build 1-to-1 mapping between two sequences
def MapTwoSequences(querySeq, refSeq):

        alignments = pairwise2.align.localds(querySeq, refSeq, blosum80, -5, -0.2)
        ##alignments = pairwise2.align.localxx(pdbseq, sequence)
        ##print '#alignments:', len(alignments)

        ##find the alignment with maximum coverage
        bestCoverage = 0
	bestMapping = None
	bestMisMatches = len(querySeq)
        for alignment in alignments:
                #print alignment[0]
                #print alignment[1]
		mapping, coverage, numMisMatches = Alignment2Mapping(alignment)
		if bestCoverage < coverage:
			bestCoverage = coverage
			bestMapping = mapping
			bestMisMatches = numMisMatches

	return bestMapping, bestMisMatches

                                                                            
def CalcPhiPsi(PDBfile, querySeq=None, chainID=None, modelNo=None):

	sequence = ""
	resIDs = []
	PhiPsi = []

	for model in Bio.PDB.PDBParser().get_structure("XXX", PDBfile):
		if modelNo is not None and model.id != modelNo:
			continue
    		for chain in model:
			if chainID is not None and chinaID != chain.id:
				continue

        		polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
        		for poly_index, poly in enumerate(polypeptides) :
				"""
            			print "Model %s Chain %s" % (str(model.id), str(chain.id)),
            			print "(part %i of %i)" % (poly_index+1, len(polypeptides)),
            			print "length %i" % (len(poly)),
            			print "from %s%i" % (poly[0].resname, poly[0].id[1]),
            			print "to %s%i" % (poly[-1].resname, poly[-1].id[1])
				"""

            			phi_psi = poly.get_phi_psi_list()
            			for res_index, residue in enumerate(poly):
					PhiPsi.append(phi_psi[res_index])
					resIDs.append(residue.id[1])
					sequence += three_to_one(residue.resname)
					"""
                			res_name = "%s%i" % (residue.resname, residue.id[1])
                			print res_name, phi_psi[res_index]
					"""

	if querySeq is None:
		return PhiPsi, sequence, resIDs

	mapping, numMismatches = MapTwoSequences(querySeq, sequence)	
	if numMismatches > 5:
		print 'ERROR: too many mismatches between the query sequence and the PDB file'
		exit(1)

	assert len(mapping) == len(querySeq)

	qPhiPsi = [ None ] * len(querySeq)
	qResIDs = [ None ] * len(querySeq)
	for i, m in zip(range(len(querySeq)), mapping):
		if m < 0:
			continue
		qPhiPsi[i] = PhiPsi[ m ]
		qResIDs[i] = resIDs[ m ]

	return qPhiPsi, querySeq, qResIDs
	


def Usage():
	print 'CalcPhiPsi.py PDBfile [chainID] [modelNumber]'
	
def main(argv):

	if len(argv) < 1:
		Usage()
		exit(1)

	PDBfile = argv[0]
	if not os.path.isfile(PDBfile):
		print 'Please provide a valid PDB file'
		exit(1)

	chainID=None
	if len(argv) >=2 :
		chainID = argv[1]

	modelNo=None
	if len(argv) >=3:
		modelNo = np.int32(argv[2])
		if modelNo < 0:
			print 'ERROR: please provide a valid model number'
			exit(1)

	PhiPsi, sequence, resIDs = CalcPhiPsi(PDBfile, chainID=chainID, modelNo=modelNo)

	## print
	for angle, res, id in zip(PhiPsi, sequence, resIDs):
		phi_degree = None
		if angle[0] is not None:
			phi_degree = angle[0]/np.pi * 180
	
		psi_degree = None
		if angle[1] is not None:	
			psi_degree = angle[1]/np.pi * 180

		print res, id, angle[0], angle[1], phi_degree, psi_degree


if __name__ == "__main__":
        main(sys.argv[1:])

