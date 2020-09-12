#!/usr/bin/env python
import numpy as np
import pickle
import sys
import os


#------- usage -------#
def Usage():
    print('python DumpTPLPKL.py <TPL_pkl> <TPL_dump> <map_dump> ')
    print('TPL_pkl  : the input file of a pkl formated TPL')
    print('TPL_dump : the dumped PDB content in plain-text')
    print('map_dump : the dumped mapping between SEQRES and DSSP')


#------- main --------#
def main(argv):
    if len(argv) != 3:
        Usage()
        sys.exit(1)

    #----- input arguments ----#
    TPL_pkl = argv[0]
    TPL_dump = argv[1]
    map_dump = argv[2]

    if not os.path.isfile(TPL_pkl):
        print('please provide a valid TPL_pkl file')
        sys.exit(1)

    #----- load TPL_pkl -----#
    with open(TPL_pkl, 'rb') as fh:
    	a = pickle.load(fh)

    #----- dump PDB content ---#
    savefh = open(TPL_dump, 'wb')
    savefh.write('  Num Res  Missing   SSE    CLE   ACC   pACC  CNa CNb   Xca       Yca       Zca       Xcb       Ycb       Zcb\n')
    length=a['length']

    validSeq = ''
    count = 0
    for i in range(length):
	#print "count=", count, "i=", i
        if a['missing'][i] == 0:
	    if a['Ca'][i] is None:
		CaX, CaY, CaZ = -999, -999, -999
	    else:
		CaX, CaY, CaZ = a['Ca'][i][0],a['Ca'][i][1],a['Ca'][i][2]
	    if a['Cb'][i] is None:
		CbX, CbY, CbZ = -999, -999, -999
	    else:
		CbX, CbY, CbZ = a['Cb'][i][0],a['Cb'][i][1],a['Cb'][i][2]

            savefh.write('{} {} {} {} {} {} {} {} {} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n'.format(count+1,a['pdbseq'][count],0,a['SS8_str'][i],a['SS3_str'][i],a['ACC'][i],a['pACC'][i],a['CNa'][i],a['CNb'][i],CaX, CaY, CaZ, CbX, CbY, CbZ))
	    validSeq += a['pdbseq'][count]
            count = count + 1
	else:
	    validSeq += '-'
	    #print 'missing: ', i
    savefh.close()

    #print validSeq, '*'

    #----- dump map string ---#
    savefh = open(map_dump, 'wb')
    savefh.write('{}\n'.format(a['sequence']))
    #savefh.write('{}\n'.format(a['DSSPsequence']))
    savefh.write('{}\n'.format(validSeq))
    savefh.close()


#-------------- python main ---------------#
if __name__ == "__main__":
    main(sys.argv[1:])


