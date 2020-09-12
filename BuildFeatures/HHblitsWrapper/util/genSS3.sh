#!/bin/bash
if [ $# -ne 2 ]
then
        echo "Usage: ./PSIPRED <input_mtx> <out_dir> "
        exit
fi

# $1 the sequence name with suffix .seq
# $2 out directory

RaptorX_HOME=~/GitBucket/casp13_pipeline/TGT_Package
PSIPREDDIR=$RaptorX_HOME/util/PSIPRED
INPUTMTX=$1
DESTDIR=$2

###### BY TINA, May 6, 2003
# first make PSP directory if it doesn't exist
if [ ! -d $DESTDIR ] ; then
    mkdir $DESTDIR
fi
###### BY TINA, May 6, 2003

fulnam=`basename $1`
bname=${fulnam%.*}
rootname=R$bname
cp $INPUTMTX $rootname.mtx

echo Pass1 ....
$PSIPREDDIR/bin/psipred $rootname.mtx $PSIPREDDIR/data/weights.dat $PSIPREDDIR/data/weights.dat2 $PSIPREDDIR/data/weights.dat3 > $rootname.ss

echo Pass2
$PSIPREDDIR/bin/psipass2 $PSIPREDDIR/data/weights_p2.dat 1 0.98 1.09 $rootname.ss2 $rootname.ss > $rootname.horiz


echo "Final output files:" $rootname.ss2 $rootname.horiz $rootname.ss
mv $rootname.ss $DESTDIR/$bname.ss
mv $rootname.ss2 $DESTDIR/$bname.ss2
mv $rootname.horiz $DESTDIR/$bname.horiz

#remove temporary files
echo Cleaning up ....
#-rm -f $rootname.pn $rootname.sn $rootname.aux error.log $rootname.mtx
#-rm -f $rootname.chk 
rm -f $rootname.*
