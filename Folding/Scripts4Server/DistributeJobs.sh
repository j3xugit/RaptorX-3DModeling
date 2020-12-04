#!/bin/bash

if [ $# -lt 3 ]; then
	echo $0 ' numAllowedJobs listFile command [arguments for command except the protein name]'
	echo '  When you run a protein with the command, the protein name shall always be the first argument of the command'
	echo '  Example: ' $0 ' 10 CASP12FM.list scripts/FoldingWrapper10.sh SeqFolder DistancePredictionFolder  PropertyPredictionFolder CbCb'
	exit
fi

myself=`basename $0 | cut -f1 -d'.'`

numAllowedJobs=$1
listFile=$2
program=$3

## use the first 10 letters as the base name
keywords=`basename $program | cut -f1 -d'.' | cut -c1-10 `

arguments=${@:4}

echo $listFile
echo $program
echo $keywords
echo $arguments


for target in `cat $listFile`
do

	echo 'now scheduling' $target...

	while true
	do
		a=`expr $RANDOM % 10 `
		sleep $a

		##check the number of running jobs
		#echo ${keywords}
		numRunningJobs=`ps -ef | grep ${keywords} | grep -v ${myself} | wc -l  `
		#echo $numRunningJobs

		if [ $numRunningJobs -lt `expr $numAllowedJobs + 1 ` ]; then
			#echo 'Running' $target...
			#echo $target > $target.list
			cmdStr="${program} $target ${arguments}"

			echo Executing $cmdStr ...

			${cmdStr} &

			break

		fi

		sleep 5

	done
done

