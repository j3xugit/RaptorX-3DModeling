#!/bin/bash

if [[ -z "${DistFeatureHome}" ]]; then
	echo "ERROR: Please set the environmental variable DistFeatureHome to the installation directory of BuildFeatures, e.g., $HOME/RaptorX-3DModeling/BuildFeatures/"
	exit 1
fi

if [ -z "${HHDIR}" ]; then
        echo "ERROR: please set environmental variable HHDIR to the install folder of HHblits"
        exit 1
fi

if [ ! -d $HHDIR ]; then
        echo "ERROR: invalid folder $HHDIR "
        exit 1
fi

DB=""

out_root=`pwd`   #-> output to current directory
coverage=-2     #-> automatic determine the coverage on basis of input sequence length 
cpu_num=4       #-> use 4 CPUs 

input_fasta=""
iteration=3
e_value=0.001

function usage
{
	echo $0 " [-n iteration] [-e evalue] [-C coverage] [-o out_root] [-c numCPUs ] [-d DB] [-h hhsuite] seqFile"
	echo "	This script run HHblits to build an MSA for a protein sequence"
	echo "	-n: the number of iterations, default $iteration"
	echo "	-e: Evalue, default ${e_value}"
	echo "	-C: -2, -1 or a positive value, default $coverage"
	echo "		if -2, do not use the -cov option in HHblits"
	echo "		if -1, automatically determine coverage"
	echo "		if a positive value, use this value for the -cov option of HHblits"
	echo "	-o: output directory, default current work directory"
	echo "	-c: the number of CPUs to be used, default ${cpu_num}"
	echo "	-d: the hhm database to be searched by HHblits, default $DB"
	echo "	-h: the install folder of hhsuite, default $HHDIR"
	exit 1
}

while getopts ":n:e:C:o:c:d:h:" opt;
do
	case $opt in
	n)
		iteration=$OPTARG
		;;
	e)
		e_value=$OPTARG
		;;
	C)
		coverage=$OPTARG
		;;
	o)
		out_root=$OPTARG
		;;
	c)
		cpu_num=$OPTARG
		;;
	d)
		DB=$OPTARG
		;;
	h)
		HHDIR=`readlink -f $OPTARG`
		;;
	#-> others
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

if [ $# -ne 1 ]; then
        usage
fi

if [ -z "$DB" ]; then
        if [ -z "${HHDB}" ]; then
                echo "ERROR: You did not provide a seq database for search by hhblits, so the default database will be used"
                echo "ERROR: However, the environemental variable HHDB for the default database is not set"
                echo "ERROR: please set HHDB to a seq database searchable by hhblits, e.g., /mnt/data/RaptorXCommon/HHblits/DB/uniref30_2020_02/uniref30_2020_02"
                exit 1
        else
                DB=$HHDB
        fi
fi

DBIndexFile=${DB}_hhm.ffindex
if [ ! -f $DBIndexFile ]; then
        echo "ERROR: invalid or damanged HMM database to be searched by hhblits: $DB"
        exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

input_fasta=$1
if [ ! -f "$input_fasta" ]; then
	echo "ERROR: invalid input sequence file: $input_fasta "
	exit 1
fi

fulnam=`basename $input_fasta`
relnam=${fulnam%.*}

mkdir -p $out_root
currDir=`pwd`

tmpWorkDir=`mktemp -d ${relnam}.tmpWork4HHBlits.XXXXXX`
cp $input_fasta $tmpWorkDir/$relnam.fasta
cd $tmpWorkDir

seq_file=$relnam.seq
$cmdDir/util/Verify_FASTA $relnam.fasta $seq_file
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run $cmdDir/util/Verify_FASTA $relnam.fasta $seq_file"
	exit 1
fi

# ----- calculate minimum alignment coverage  ----
if [ $coverage -eq -1 ]; then
	## the default minimum coverage percentage is set to 60%
	a=60

	## b is the minimum percentage of coverage such that at least 80 residues of the query sequence are covered by homologs
	## $3 in the below sentence is the query sequence length
	b=`tail -n1 $seq_file | wc | awk '{print int(7000/($3-1))}'`
	if [ $a -gt $b ]; then
		coverage=$b
	else
		coverage=$a
	fi
fi

a3m_file=$relnam.a3m
echo "Searching database $DB with evalue $e_value and $iteration iterations"

if [ $coverage -eq -2 ]; then
	echo "Running HHblits without specifying coverage ..."
	$HHDIR/bin/hhblits -i $seq_file -cpu $cpu_num -d $DB -o $relnam.hhr -oa3m $relnam.a3m -n $iteration -e $e_value
	ret=$?
else
	echo "Running HHblits with -maxfilt 500000 -diff inf -id 99 -cov $coverage..."
	$HHDIR/bin/hhblits -i $seq_file -cpu $cpu_num -d $DB -o $relnam.hhr -oa3m $relnam.a3m -n $iteration -e $e_value -maxfilt 500000 -diff inf -id 99 -cov $coverage
	ret=$?
fi
if [ $ret -ne 0 ]; then
	echo "ERROR: failed to run  $HHDIR/bin/hhblits -i $seq_file -cpu $cpu_num -d $DB -o $relnam.hhr -oa3m $relnam.a3m -n $iteration"
	rm -f $relnam.*
	exit 1
fi

cd $currDir

cp $input_fasta $out_root/$relnam.fasta_raw
mv $tmpWorkDir/$seq_file $out_root
mv $tmpWorkDir/$a3m_file $out_root
rm -rf $tmpWorkDir
