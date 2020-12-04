#!/bin/bash

if [ $# -lt 2 ]; then
	echo "$0 targetListFile MSADir [ gpu ]"
	echo "gpu: -1 (default), 0, 1, 2 and 3. If -1, then choose one gpu automatically"
	exit 1
fi

targets=$1
MSADir=$2
GPU=-1
if [ $# -ge 3 ]; then
	GPU=$3
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
parentDir=`dirname $cmdDir`

bMSADir=`basename $MSADir`
DestDir=`echo $bMSADir | sed -r "s/MSA/Features4Train/"`
mkdir -p $DestDir

program=$parentDir/GenDistFeaturesFromMSA.sh 
numAllowedJobs=12
keywords=`basename $program`
myself=`basename $0 `

for i in `cat $targets`
do
	while true
	do
		 ## check the number of running jobs
                numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l`
                if [ $numRunningJobs -lt $numAllowedJobs ]; then
			if [ -f $MSADir/${i}.a3m ]; then
				$program -o $DestDir/feat_${i}_contact -g $GPU $MSADir/${i}.a3m &
			fi
			break
                else
                        a=`expr $RANDOM % 3 `
                        sleep $a
                fi

        done

done

wait
