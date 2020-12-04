#!/bin/bash

if [ $# -lt 1 ]; then
	echo $0 MSAfile
	echo "MSAfile: MSA file in a3m format"
	exit 1
fi

a3mfile=$1
if [ ! -f $a3mfile ]; then
	echo invalid MSA file $a3mfile
	exit 1
fi

numSeqs=`wc -l $a3mfile | cut -f1 -d' '`
numSeqs=` expr $numSeqs / 2 `
seqLen=`head -2 $a3mfile | tail -1 | wc -c | cut -f1 -d' '`

## we add 500M as the slack
neededRAM=`expr 320 \* $seqLen + 8500 \* $seqLen \* $seqLen + $numSeqs \* 94 \* $seqLen +  $numSeqs \* 4 + 500000000 `

echo $neededRAM
