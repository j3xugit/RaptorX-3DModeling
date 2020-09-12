#!/bin/bash
if [ $# -ne 2 ]
then
        echo "Usage: ./DISOPRED <input_mtx> <out_dir> "
        exit
fi

# $1 the sequence name with suffix .seq
# $2 out directory

RaptorX_HOME=~/GitBucket/casp13_pipeline/TGT_Package
PSIPREDDIR=$RaptorX_HOME/util/DISOPRED
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

echo Pass1 ...
echo Pass2 ...
$PSIPREDDIR/bin/disopred $rootname $rootname.mtx $PSIPREDDIR/data/


echo "Final output files:" $rootname.diso $rootname.horiz_d
mv $rootname.diso $DESTDIR/$bname.diso
mv $rootname.horiz_d $DESTDIR/$bname.horiz_d

#remove temporary files
echo Cleaning up ....
#-rm -f $rootname.pn $rootname.sn $rootname.aux error.log $rootname.mtx
#-rm -f $rootname.chk 
rm -f $rootname.*
