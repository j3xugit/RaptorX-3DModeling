#!/bin/bash

if [ $# -lt 1 ]; then
	echo $0 MSAfile
	exit 1
fi

input=$1
if [ ! -f $input ]; then
	echo "ERROR: invalid MSA file $input"
	exit 1
fi

OneM=1048576
## maxGPURAM unit is Mbytes
maxGPURAM=`nvidia-smi --query-gpu=memory.total --format=csv | tail -n +2 | cut -f1 -d' ' | sort -rn | head -1`
host=`hostname`
#if [[ "$host" == "raptorx7.uchicago.edu" ]]; then
#        maxGPURAM=12180
#fi
maxGPURAM=`expr $maxGPURAM - 800 `
#echo maxGPURAM=$maxGPURAM

## convert unit of maxGPURAM to bytes
maxGPURAM=`expr $maxGPURAM \* $OneM `

## estimate the allowed number of sequences based upon sequence length and GPU RAM limit
seqLen=`head -2 $input | tail -1 | wc -c | cut -f1 -d' '`
allowedRAM=`expr $maxGPURAM - 320 \* $seqLen - 10756 \* $seqLen \* $seqLen `
tmpvalue=`expr 94 \* $seqLen + 4 `

#echo $allowedRAM, $tmpvalue
numAllowedSeqs=`expr $allowedRAM / $tmpvalue `

## looks like that CCMpred cannot handle more than 250k sequences even if GPU has sufficient memory
## here we allow at most 200k sequences
if [ $numAllowedSeqs -gt 200000 ]; then
        numAllowedSeqs=200000
fi

## we do not want to remove too many sequences from the MSA
if [ $numAllowedSeqs -lt 10000 ]; then
        numAllowedSeqs=10000
fi

numLines=`wc -l $input | cut -f1 -d' '`
numLines=`expr $numLines / 2 `

if [ $numLines -lt $numAllowedSeqs ]; then
        numAllowedSeqs=$numLines
fi

echo $numAllowedSeqs
