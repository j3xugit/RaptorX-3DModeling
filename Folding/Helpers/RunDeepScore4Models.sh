#!/bin/bash

if [ $# -lt 2 ]; then
	echo "Usage: $0 modelListFile nativePDB [ID (default RX)], where ID can be used to identify one method"
	exit 1
fi

native=$2
modelList=$1
ID=RX

if [ $# -ge 3 ]; then
	ID=$3
fi

target=`basename $native .pdb `

result=$target.${ID}.quality.txt
cp /dev/null $result

for m in `cat $modelList`
do
        scores=`DeepScore $m $native | python $DistanceFoldingHome/Helpers/CollectModelQuality.py `
        echo $scores $m >> $result
done

sort -k5,5 -rn $result  > $result.sorted
mv $result.sorted $result

