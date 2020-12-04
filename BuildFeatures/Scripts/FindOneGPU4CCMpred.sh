#!/bin/bash

if [ $# -lt 1 ]; then
	echo $0 "a2m_file [mins_to_wait] "
	echo "	this script finds one GPU with sufficient memory "
	echo "	a3m_file: the a3m file to be run by CCMpred"
	echo "	mins_to_wait (default 180): #minutes to wait for an appropirate GPU (i.e., a GPU with sufficient memory for your job)"
	exit 1
fi

a2mfile=$1
if [ ! -f $a2mfile ]; then
	echo "ERROR: FindOneGPU.sh cannot find the file: $a2mfile"
	exit 1
fi

## in total wait for 120 mins
patience=120
if [ $# -ge 2 ]; then
	patience=$2
fi

numSeqs=`wc -l $a2mfile | cut -f1 -d' '`
seqLen=`head -1 $a2mfile | wc -c | cut -f1 -d' '`

## we add 500M as the slack 
neededRAM=`expr 320 \* $seqLen + 8500 \* $seqLen \* $seqLen + $numSeqs \* 94 \* $seqLen +  $numSeqs \* 4 + 500000000 `

#echo "neededRAM=$neededRAM"

gpu=-1

OneM=1048576
## maxGPURAM unit is Mbytes
maxGPURAM=`nvidia-smi --query-gpu=memory.total --format=csv | tail -n +2 | cut -f1 -d' ' | sort -rn | head -1`
maxGPURAM=`expr $maxGPURAM - 500 `
maxGPURAM=`expr $maxGPURAM \* $OneM `
#echo maxGPURAM=$maxGPURAM

## if neededRAM is too big, then just exit
if [ $neededRAM -gt $maxGPURAM ]; then
        echo $gpu
        exit
fi


## t is the time spent so far
t=0
while [ $t -le $patience ]
do

	## got the GPU with the maximum amount of free memory
	bestGPU=`nvidia-smi --query-gpu=index,memory.free --format=csv | tail -n +2 | cut -f1,2 -d' ' | sort -k2,2 -rn | head -1 `
	# echo $bestGPU

	index=`echo $bestGPU | cut -f1 -d',' `
	freeM=`echo $bestGPU | cut -f2 -d',' `
	fmBytes=`expr $freeM  \* 1048576 `

	#echo $index $fmBytes

	if [ $neededRAM -lt $fmBytes ]; then
		gpu=$index
		break
	fi

	## wait for 1 minute to try again
	sleep 1m
	t=`expr $t + 2 `

done

#echo "DetectedGPU=$gpu"
echo $gpu
