#!/bin/sh

if [ $# -lt 2 ]; then
	echo "Usage: $0 targetName NativeDir "
	exit 1
fi

target=$1
NativeDir=$2

result=${target}.quality
cp /dev/null $result

for i in Run*
do
	for model in $i/Out/S*pdb
	do
		native=${NativeDir}/${target}.pdb
		scores=`DeepScore $model $native | python $DistanceFoldingHome/Helpers/CollectModelQuality.py`
		echo $scores $model >> $result
	done
done

sort -k5,5 -rn $result > $result.sorted

rm -f $result

