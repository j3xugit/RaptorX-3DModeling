#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 "proteinListFile DecoyDir [ AtomPairType [ boundType [ cutoff [ angleConstant ] ] ] ]"
	echo "  default values are All, 3, 15 and 5.0 "
	echo "Example: " $0 "ProteinLists/CASP12DMTPL-Top2.list ModelingResult All 3 15"
	exit
fi

proteins=$1

DataName=`basename $proteins .list`
realDataName=CASP12DM

DecoyDir=$2
DecoyDirBase=`basename $DecoyDir`

atomPairType=All
if [ $# -ge 3 ]; then
	atomPairType=$3
fi

boundType=t3
if [ $# -ge 4 ]; then
	boundType=t$4
fi

cutoff=c15
if [ $# -ge 5 ]; then
	cutoff=c$5
fi

angleFlag=A5.0
if [ $# -ge 6 ]; then
	angleFlag=A$6
fi


method=${atomPairType}_${boundType}_${cutoff}_${angleFlag}

resultfile=${DataName}-${DecoyDirBase}-${method}.modelQuality.txt
cp /dev/null $resultfile

NativeDir=${realDataName}-PDB/

for target in `cat $proteins`
do
	ID=${target}.${DecoyDirBase}-${method}
	scorefile=${ID}.deepscore
	cp /dev/null $scorefile

	part1=`echo $target | cut -f1 -d'-'`
	part2=`echo $target | cut -f2 -d'-'`
	realTargetName=${part1}-${part2}

	native=${NativeDir}/${realTargetName}.pdb
	if [ ! -f $native ]; then
		echo "the native file does not exist: " $native
		exit -1
	fi

	for i in 1 2 3 4 5
	do
		model=${DecoyDir}/${target}_${method}/${target}_model${i}.pdb 
		if [ ! -f $model ]; then
			echo "the model file does not exist: " $model
			#exit -1
		else
			scores=`DeepScore $model $native | python scripts/CollectModelQuality.py`
			echo $scores $model >> $scorefile
		fi

	done

	if [ -s $scorefile ]; then
		python scripts/AnalyzeDeepScore.py $scorefile >> $resultfile
	fi

	rm -f $scorefile
done


top1TM=`  awk '{ total += $3 }  END {print total/NR }' $resultfile `
top1GDT=` awk '{ total += $5 }  END {print total/NR }' $resultfile `
top1GHA=` awk '{ total += $6 }  END {print total/NR }' $resultfile `
top5TM=`  awk '{ total += $8 }  END {print total/NR }' $resultfile `
top5GDT=` awk '{ total += $10 } END {print total/NR }' $resultfile `
top5GHA=` awk '{ total += $11 } END {print total/NR }' $resultfile `

echo $top1TM $top1GDT $top1GHA $top5TM $top5GDT $top5GHA

