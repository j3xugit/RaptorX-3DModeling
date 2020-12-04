#!/bin/bash

if [ $# -lt 1 ]; then
	echo $0 rawScoreFile
	exit 1
fi

resultfile=$1

TM=` awk '{ total += $4 }  END {print total/NR }' $resultfile `
GDT=` awk '{ total += $5 }  END {print total/NR }' $resultfile `
TMGDT=` awk '{ total += $6 } END {print total/NR }' $resultfile `
GHA=` awk '{ total += $7 }  END {print total/NR }' $resultfile `
topTM=` awk '{ total += $9 }  END {print total/NR }' $resultfile `
topGDT=` awk '{ total += $10 }  END {print total/NR }' $resultfile `
topTMGDT=` awk '{ total += $11 } END {print total/NR }' $resultfile `
topGHA=` awk '{ total += $12 } END {print total/NR }' $resultfile `

numLines=`wc -l $resultfile | cut -f1 -d' ' `
echo $numLines $TM $GDT $TMGDT $GHA $topTM $topGDT $topTMGDT $topGHA

