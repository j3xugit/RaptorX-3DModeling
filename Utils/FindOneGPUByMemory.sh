#!/bin/bash

if [ $# -lt 1 ]; then
	echo $0 "RAM_needed_bytes [mins_to_wait]"
	echo "	This script checks the GPU usage every minute to find one GPU with sufficient memory "
	echo "	RAM_needed: the mount of GPU memory (in bytes) needed by a job"
	echo "	mins_to_wait: the number of minutes to be waited, default 90"
	exit 1
fi

neededRAM=$1

## in total wait for 90 mins
patience=90
if [ $# -ge 2 ]; then
	patience=$2
fi


OneM=1048576
## maxGPURAM unit is Mbytes
maxGPURAM=`nvidia-smi --query-gpu=memory.total --format=csv | tail -n +2 | cut -f1 -d' ' | sort -rn | head -1`
maxGPURAM=`expr $maxGPURAM - 500 `
maxGPURAM=`expr $maxGPURAM \* $OneM `
#echo maxGPURAM=$maxGPURAM

gpu=-1
## if neededRAM is too big, then just exit
if [ $neededRAM -gt $maxGPURAM ]; then
	echo $gpu
	exit 
fi

## t is the time spent so far
t=0
while [ $t -lt $patience ]
do

	## got the GPU with the maximum amount of free memory
	bestGPU=`nvidia-smi --query-gpu=index,memory.free --format=csv | tail -n +2 | cut -f1,2 -d' ' | sort -k2,2 -rn | head -1 `
	# echo $bestGPU

	index=`echo $bestGPU | cut -f1 -d',' `
	freeM=`echo $bestGPU | cut -f2 -d',' `
	fmBytes=`expr $freeM  \* $OneM `
	#echo $index $fmBytes

	if [ $neededRAM -lt $fmBytes ]; then
		gpu=$index
		break
	fi

	## wait for 1 minute to try again
	sleep 1m
	t=`expr $t + 1 `
done

echo $gpu
