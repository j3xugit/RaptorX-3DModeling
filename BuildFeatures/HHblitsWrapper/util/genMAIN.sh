#!/bin/bash
if [ $# -ne 2 ]
then
        echo "Usage: ./genMAIN.sh <input_a3m> <out_dir> "
        exit
fi

#----- input files ------#
INPUTA3M=$1
DESTDIR=$2
fulnam=`basename $1`
bname=${fulnam%.*}
rootname=R$bname


#------ directory ------#
RaptorX_HOME=~/GitBucket/casp13_pipeline/TGT_Package
A3M_To_PSI=$RaptorX_HOME/util/A3M_To_PSI
MSA_To_PSSM=$RaptorX_HOME/util/MSA_To_PSSM
HHMAKE=$RaptorX_HOME/util/hhmake

#----- generate PSP and MTX -----#
$A3M_To_PSI $INPUTA3M $rootname.psi_tmp
grep -v "ss_pred\|ss_conf\|ss_dssp" $rootname.psi_tmp > $rootname.psi
$MSA_To_PSSM -i $rootname.psi -o $rootname.psp -m $rootname.mtx -c 20
$HHMAKE -i $INPUTA3M -o $rootname.hhm
rm -f $rootname.psi_tmp

#----- move generated files to output ---#
mv $rootname.psi $DESTDIR/$bname.psi
mv $rootname.psp $DESTDIR/$bname.psp
mv $rootname.mtx $DESTDIR/$bname.mtx
mv $rootname.hhm $DESTDIR/$bname.hhm

