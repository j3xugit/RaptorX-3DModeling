#!/bin/bash

if [[ -z "${ModelingHome}" ]]; then
	echo "ERROR: please set the environmental variable ModelingHome to the install folder of RaptorX-3DModeling"
	exit 1
fi

RESDIR=`pwd`

numAllowedJobs=`nproc --all`

function Usage 
{
	echo $0 "[ -n numJobs | -d ResDir ] proteinListFile SeqDir PDBDir"
	echo "	This script calculates the property/distance/orientation ground truth by sequence and structure file"
	echo "	proteinListFile: a file for a list of proteins, each in one row"
	echo "	SeqDir: the folder for protein sequence files in FASTA format"
	echo "	PDBDir: the folder for protein structure file, each ending with .pdb or .cif"
	echo "	numJobs: the number of proteins to be simultanously run, default $numAllowedJobs"
	echo "	ResDir: the folder for result saving, default current work directory"
}

while getopts ":n:d:" opt; do
        case ${opt} in
                n )
                  numAllowedJobs=$OPTARG
                  ;;
                d )
                  RESDIR=$OPTARG
                  ;;
                \? )
                  echo "Invalid Option: -$OPTARG" 1>&2
                  exit 1
                  ;;
                : )
                  echo "Invalid Option: -$OPTARG requires an argument" 1>&2
                  exit 1
                  ;;
        esac
done
shift $((OPTIND -1))

if [ $# -ne 3 ]; then
	Usage
	exit 1
fi

list=$1
if [ ! -f $list ]; then
	echo "ERROR: invalid protein list file $list"
	exit 1
fi

SEQDIR=$2
if [ ! -d $SEQDIR ]; then
	echo "ERROR: invalid protein seq folder $SEQDIR"
	exit 1
fi

PDBDIR=$3
if [ ! -d $PDBDIR ]; then
	echo "ERROR: invalid folder for protein structure files $PDBDIR"
	exit 1
fi

if [ ! -d $RESDIR ]; then
	mkdir -p $RESDIR
fi

program=$ModelingHome/Common/CalcGroundTruthFromSeqPDB.py

keywords=`basename $program`
myself=`basename $0`
	
for i in `cat $list `
do
	while true
	do
                numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l`
                if [ $numRunningJobs -lt $numAllowedJobs ]; then
			pdbfile=${PDBDIR}/${i}.pdb
                        ciffile=${PDBDIR}/${i}.cif
                        if [ -f $pdbfile ]; then
                                structfile=$pdbfile
                        else
                                structfile=$ciffile
                        fi

        		python $program $structfile $SEQDIR/${i}.fasta $RESDIR &
			sleep 1
                        break
                else
                        a=`expr $RANDOM % 4 `
                        sleep $a
                fi
	done
done

wait
