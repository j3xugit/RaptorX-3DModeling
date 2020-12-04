#!/bin/bash

resultfile=$1

top1TM=`  awk '{ total += $3 }  END {print total/NR }' $resultfile `
top1GDT=` awk '{ total += $5 }  END {print total/NR }' $resultfile `
top1GHA=` awk '{ total += $6 }  END {print total/NR }' $resultfile `
top5TM=`  awk '{ total += $8 }  END {print total/NR }' $resultfile `
top5GDT=` awk '{ total += $10 } END {print total/NR }' $resultfile `
top5GHA=` awk '{ total += $11 } END {print total/NR }' $resultfile `

echo $top1TM $top1GDT $top1GHA $top5TM $top5GDT $top5GHA

