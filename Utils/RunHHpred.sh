#!/bin/bash

numIterations=2
Evalue=0.001
neffmax=6.5

numCPUs=4
ResDir=`pwd`

if [ -z "${HHDIR}" ]; then
        echo "ERROR: please set HHDIR to the install folder of hhblits, e.g., /mnt/data/RaptorXCommon/HHblits/hhsuite-3.2.0-SSE2-Linux"
        exit 1
fi
HHbin=$HHDIR/bin

if [ -z "${PDB70HHM}" ]; then
        echo "ERROR: please set PDB70HHM to the install folder of PDB70 template database provided by hhsuite, e.g., /mnt/data/RaptorXCommon/HHblits/DB/PDB70/pdb70"
	echo "The database can be downloaded from http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/"
        exit 1
fi

DB=$PDB70HHM

function Usage 
{
        echo $0 "[-i nIterations | -e evalue | -m maxneff | -n numCPUs | -d DB | -s savefolder ] inFile"
	echo "	This script builds an .hhm file for query sequence and then searches it against an HHM database (PDB70 by default) downloaded from http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/"
	echo "	inFile: a protein seq file in FASTA format, ending with .seq or .fasta, or an MSA file ending with .a3m or an HHM file ending with .hhm"
	echo "	-d: the template database to be searched by HHblits, default $DB"
	echo "	-i: the number of iterations searched by hhblits for MSA generation, default $numIterations"
	echo "	-m: maximum neff value used by hhblits for MSA generation, default $neffmax"
        echo "	-n: the number of CPUs to be used by hhblits for MSA generation, default $numCPUs"
	echo "	-e: E-value used by hhblits for HHM database search, default $Evalue"
	echo "	-s: the folder for result saving, default current work directroy"
}

while getopts ":i:n:e:m:s:d:" opt; do
        case ${opt} in
                i )
                  numIterations=$OPTARG
                  ;;
                n )
                  numCPUs=$OPTARG
                  ;;
		e )
		  Evalue=$OPTARG
		  ;;
		m )
		  neffmax=$OPTARG
		  ;;
		d )
		  DB=$OPTARG
		  ;;
                s )
                  ResDir=$OPTARG
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

if [ $# -ne 1 ]; then
        Usage
        exit 1
fi

seqFile=$1
if [ ! -f $seqFile ]; then
	echo "ERROR: invalid sequence file $seqFile"
	exit 1
fi

if [ ! -d $ResDir ]; then
	mkdir -p $ResDir
fi

hhmNeeded=0
if [[ $seqFile == *.fasta ]]; then
        seqName=`basename $seqFile .fasta `
	hhmNeeded=1

elif [[ $seqFile == *.a3m ]]; then
	seqName=`basename $seqFile .a3m`
	hhmake -i $seqFile -o $ResDir/$seqName.hhm

elif [[ $seqFile == *.hhm ]]; then
	seqName=`basename $seqFile .hhm`
	if [ ! -f $ResDir/$seqName.hhm ]; then
		cp $seqFile $ResDir/$seqName.hhm
	fi
else
        seqName=`basename $seqFile .seq `
	hhmNeeded=1
fi


cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

if [ ! -f ${DB}_hhm.ffindex ]; then
	echo "ERROR: invalid or damaged template database to be searched by HHblits: $DB"
	exit 1
fi

if [ $hhmNeeded -eq 1 ]; then
	## generate .hhm file
	$cmdDir/BuildHHM.sh -i $numIterations -n $numCPUs -e $Evalue -m $neffmax -s $ResDir $seqFile
	if [ $? -ne 0 ]; then
		echo "ERROR: failed to generate a .hhm file for $seqFile"
		exit 1
	fi
fi
hhmfile=$ResDir/$seqName.hhm

## search the template database
fname=`basename $DB`
hhrfile=$ResDir/${seqName}.$fname.hhr
$HHbin/hhsearch -i $hhmfile -d $DB -mact 0.1 -o $hhrfile -z 300 -b 300
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run $HHbin/hhsearch -i $hhmfile -d $DB -o $hhrfile -z 300 -b 300"
	exit 1
fi

## generate individual alignments
python $ModelingHome/Utils/ParseAlignmentFromHHR.py -d $ResDir $hhrfile $ResDir/${seqName}.hhm
