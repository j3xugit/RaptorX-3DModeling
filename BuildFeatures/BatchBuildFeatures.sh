#!/bin/bash

numAllowedJobs=3
ResDir=`pwd`
gpu=-1

MSAmethod=25

function Usage
{
	echo $0 "[ -o outDir | -g gpu | -m MSAmethod | -n numJobs ] proteinList SeqDir "
	echo "	This script builds contact/distance/orientation/property/threading prediction features for a set of proteins"
	echo "	SeqDir: the folder for sequence files, each in FASTA format"
	echo "	-o: the folder for result saving, default current work directory"
	echo "	-g: -1 (default), 0-3. If -1, automatically select one GPU"
	echo "	-m: the MSA generation method, default $MSAmethod, See BuilFeatures.sh and BuildMSAs.sh for a detailed explanation"
	echo "	-n: the number of sequences to be run simultaneously, default $numAllowedJobs"
	echo "	    Note that all the jobs may compete for the same GPU with a limited amount of memory, so please use a small number"
	echo "	    Set gpu=-1 will allow these jobs to use more than 1 GPUs if available"
}

while getopts ":o:g:n:m:" opt; do
        case $opt in
        o)
                ResDir=$OPTARG
                ;;
        g)
                gpu=$OPTARG
                ;;
	m)	
		MSAmethod=$OPTARG
		;;
       	n)
        	numAllowedJobs=$OPTARG
        	;;
        #-> help
        \?)
                echo "Invalid option: -$OPTARG" >&2
                exit 1
                ;;
        :)
                echo "Option -$OPTARG requires an argument." >&2
                exit 1
                ;;
        esac
done
shift $((OPTIND -1))

if [ $# -ne 2 ]; then
	Usage
	exit 1
fi

proteinList=$1
if [ ! -f $proteinList ]; then
	echo "ERROR: invalid protein list file $proteinList "
	exit 1
fi

SeqDir=$2
if [ ! -d $SeqDir ]; then
	echo "ERROR: invalid sequence folder $SeqDir"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/BuildFeatures.sh
if [ ! -x $program ]; then
	echo "ERROR: $program not executable"
	exit 1
fi

keywords=`basename $program`
myself=`basename $0 `

for target in `cat $proteinList `
do
        while true
        do
                numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l`
                if [ $numRunningJobs -lt $numAllowedJobs ]; then
			$program -o $ResDir -g $gpu -m $MSAmethod $SeqDir/$target.fasta &
			sleep 1
                        break
                else
                        a=`expr $RANDOM % 4 `
                        sleep $a
                fi
        done
done

wait
