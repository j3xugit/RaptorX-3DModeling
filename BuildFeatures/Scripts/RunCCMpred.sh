#!/bin/bash

if ! type "nvidia-smi" > /dev/null; then
	echo "ERROR: looks like GPU not installed on this machine"
	exit 1
fi

if [[ -z "${ModelingHome}" ]]; then
	echo "ERROR: please set environmental variable ModelingHome to the installation folder of RaptorX-3DModeling"
	exit 1
fi

if [[ -z "${DistFeatureHome}" ]]; then
        echo "ERROR: please set environmental variable DistFeatureHome to the installation folder of BuildFeatures"
        exit 1
fi

outDir=`pwd`
gpu=-1

function Usage
{
	echo $0 "[ -d outDir | -g gpu ] a3mfile"
	echo "	This script runs CCMpred on a local GPU"
	echo "	a3mfile: a multiple sequence alignment file in a3m format"
	echo "	outDir: the folder for result saving, default current work directory"
	echo "	gpu: -1 (default and recommended), 0-3. If -1, select a GPU automatically"
}

while getopts ":d:g:" opt; do
        case ${opt} in
                d )
                  outDir=$OPTARG
                  ;;
                g )
                  gpu=$OPTARG
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

input=`readlink -f $1`
if [ ! -f $input ]; then
	echo "ERROR: invalid input a3m file $input"
	exit 1
fi
fulnam=`basename $input`
target=${fulnam%.*}

if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi

output=$outDir/$target.ccmpred
mpkout=$output.mpk
a2mfile=$outDir/${target}.a2m
tmpa2mfile=$a2mfile.tmp

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

numAllowedSeqs=`$cmdDir/EstimateNumAllowedSeqs.sh $input`
numLines=`wc -l $input | cut -f1 -d' '`
numLines=`expr $numLines / 2 `

if [ $numLines -gt $numAllowedSeqs ]; then
	echo "WARNING: #sequences in a3m: " $numLines ", exceeds #allowed sequences: " $numAllowedSeqs
        $DistFeatureHome/bin/meff_filter -i $input -o $outDir/$target.a3m_filter -n $numAllowedSeqs -c 150000
	#python $DistFeatureHome/Helpers/SampleA3MByNumber.py $input $numAllowedSeqs $outDir/${target}.a3m_filter
        $DistFeatureHome/util/A3M_To_PSI $outDir/$target.a3m_filter $tmpa2mfile
	ecode=$?
else
        $DistFeatureHome/util/A3M_To_PSI $input $tmpa2mfile
	ecode=$?
fi

if [ $ecode -ne 0 ]; then
	echo "ERROR: failed to run util/A3M_To_PSI $input $tmpa2mfile"
        exit 1
fi
grep -v "ss_pred\|ss_conf" $tmpa2mfile | awk '{print substr($0,34,length($0)-32) }' > $a2mfile
rm -f $tmpa2mfile

## find one GPU
if [ $gpu -lt 0 ]; then
	numSeqs=`wc -l $a2mfile | cut -f1 -d' '`
        seqLen=`head -1 $a2mfile | wc -c | cut -f1 -d' '`
	neededRAM=`expr 320 \* $seqLen + 8500 \* $seqLen \* $seqLen + $numSeqs \* 94 \* $seqLen +  $numSeqs \* 4 + 500000000 `
	#echo "neededRAM=$neededRAM"

	OneM=1048576
	## maxGPURAM unit is Mbytes
	maxGPURAM=`nvidia-smi --query-gpu=memory.total --format=csv | tail -n +2 | cut -f1 -d' ' | sort -rn | head -1`
	orgMaxGPURAM=`expr $maxGPURAM \* $OneM`

	if [ $neededRAM -gt $orgMaxGPURAM ]; then
		echo "ERROR: insufficient GPU memory for running CCMpred on $input "
		exit 1
	fi

	gpu=`$ModelingHome/Utils/FindOneGPUByMemory.sh $neededRAM 30`
        if [ $? -ne 0 ]; then
                echo "ERROR: failed to run $ModelingHome/Common/FindOneGPUByMemory.sh $neededRAM 30"
                exit 1
        fi

	if [ $gpu -lt 0 ]; then
		echo "ERROR: cannot find an appropriate GPU for CCMpred on $input"
		exit 1
	fi
fi

echo "Running CCMpred on gpu $gpu of `hostname`..."
$DistFeatureHome/CCMpred/bin/ccmpred -d $gpu -b $mpkout -R $a2mfile $output
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run CCMpred on $outDir/$a2m_file at GPU $gpu of `hostname` "
        exit 1
fi
