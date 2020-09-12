#!/bin/sh

if [ $# -lt 3 ]; then
	echo "Usage: $0 listFile ModelDir NativeDir [method]"
	exit 1
fi

listFile=$1
ModelDir=$2
NativeDir=$3

if [ ! -d $ModelDir ]; then
        echo "$ModelDir does not exist"
        exit 1
fi

if [ ! -d $NativeDir ]; then
        echo "$NativeDir does not exist"
        exit 1
fi


method=CbCb_t3_c15_A5.0
if [ $# -ge 4 ]; then
	method=$4
fi

pid=$$

for i in `cat $listFile`
do
	scoreFile=$i.deepscore.$pid
	$DistanceFoldingHome/Helpers/CalcModelQuality.sh $i $ModelDir $NativeDir $method > $scoreFile &
done

wait

for i in `cat $listFile`
do
	cat $i.deepscore.$pid
	rm -f $i.deepscore.$pid
done
