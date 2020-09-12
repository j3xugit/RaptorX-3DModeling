from Bio.PDB import *
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import three_to_one

## select CG atoms for an amino acid
def SelectCG(aa, bUseAlternativeCG=True):
	assert len(aa) == 1
        AA = aa.upper()
        a2 = 'CG'
        if AA == 'V' or AA == 'I':
                a2 = 'CG1'
        elif AA == 'T':
                a2 = 'CG2'
        elif AA == 'S':
                a2 = 'OG'
        elif AA == 'C':
                a2 = 'SG'
        elif AA == 'A':
		if bUseAlternativeCG:
                	a2 = 'CB'
		else:
			a2 = ''
        elif AA == 'G':
		if bUseAlternativeCG:
                	a2 = 'CA'
		else:
			a2 = ''
        return a2

def SelectCB(AA, bUseAlternativeCB=True):
	assert len(AA) == 1
	if AA.upper() == 'G':
		if bUseAlternativeCB:
			return 'CA'
		else:
			return ''
	else:
		return 'CB'

## select atom pairs for two residues i and j in a sequence
def SelectAtomPair(sequence, i, j, atomPairType):

        if atomPairType == 'CaCa':
                return 'CA', 'CA'

        if atomPairType == 'NO':
                return 'N', 'O'

        if atomPairType == 'CbCb':
                a1 = SelectCB(sequence[i])
                a2 = SelectCB(sequence[j])
                return a1, a2

        if atomPairType == 'CaCg':
                a1 = 'CA'
                a2 = SelectCG(sequence[j])
                return a1, a2

        if atomPairType == 'CgCg':
                a1 = SelectCG(sequence[i])
                a2 = SelectCG(sequence[j])
                return a1, a2

	return None

## select Ca, N, O, Cb and Cg atoms for one residue
def SelectAtoms4Distance(res):

	ca = None
        n = None
        o = None

        if res.has_id('CA'):
        	ca = res['CA']
        if res.has_id('N'):
                n = res['N']
        if res.has_id('O'):
                o = res['O']

	if ca is None:
		print 'WARNING: Ca atom does not exist in residue: ', res
	if n is None:
		print 'WARNING: N atom does not exist in residue: ', res
	if o is None:
		print 'WARNING: O atom does not exist in residue: ', res

        if res.get_resname().upper() == 'GLY':
                cb = ca
        else:
                if res.has_id('CB'):
                	cb = res['CB']
                else:
			print 'WARNING: Cb atom does not exist in residue: ', res
                        cb = None

        cg = SelectCG( three_to_one(res.get_resname()) ).upper()
        if res.has_id(cg):
        	cg = res[cg]
        else:
		print 'WARNING: CG atom does not exist in residue: ', res
                cg = None

        return ca, n, o, cb, cg

def SelectAtoms4Orientation(sequence, i, j, labelName):
	seqLen = len(sequence)
	resName1 = sequence[i]
	resName2 = sequence[j]

	if labelName == 'Ca1Cb1Cb2Ca2':
		if resName1.upper() == 'G' or resName2.upper() == 'G':
			return None
		return 'Ca', 'Cb', 'Cb', 'Ca'

	if labelName == 'N1Ca1Cb1Cb2':
		if resName1.upper() == 'G':
			return None
		return 'N', 'Ca', 'Cb', SelectCB(resName2)

	if labelName == 'Ca1Cb1Cb2':
		if resName1.upper() == 'G':
			return None
		return 'Ca', 'Cb', SelectCB(resName2)

	if labelName == 'Ca1Ca2Ca3Ca4' or labelName == 'Ca1Ca2Ca4Ca3':
		if i == seqLen-1 or j == seqLen-1:
			return None
		return 'Ca', 'Ca', 'Ca', 'Ca'

	if labelName == 'Ca1Ca2Ca3' or labelName == 'Ca1Ca2Ca4':
		if i == seqLen-1:
			return None
		return 'Ca', 'Ca', 'Ca'

	return None
