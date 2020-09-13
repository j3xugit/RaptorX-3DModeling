#!/bin/bash

if [[ -z "${DistFeatureHome}" ]]; then
	echo "Please set the environmental variable DistFeatureHome to its install folder, e.g, $HOME/RaptorX-3DModeling/BuildFeautres "
	exit 1
fi

uniref90=${JackDB}
out_root=`pwd`
cpu_num=5 
input_fasta=""
iteration="3"
#bitscore="1"
#threshold="0.5"
bitscore="0"
threshold="0.001"

function usage
{
	echo $0 "[-n iter | -b bitsco | -t thres | -o out_root | -c numCPUs | -d uniref90] input_fasta"
	echo "	-n: the number of iterations (default $iteration)"
	echo "  -c: the number of CPUs to be used for homology search (default $cpu_num)"
	echo "  -b: using BitScore or E-value (default) as threshold, 1 for Bitscore and 0 (default) for E-value"
	echo "  -t: threshold score for BitScore or E-value (default 0.001). When BitScore is used, 0.5 may be a good value"
	echo "  -o: output folder, default current work directory "
	echo "  -d: The file path for the sequence database to be searched, default $uniref90"
}


JackHome=$DistFeatureHome/EVAlign

while getopts ":i:n:b:t:o:c:d:" opt;
do
        case $opt in
        n)
                iteration=$OPTARG
                ;;
	b)
		bitscore=$OPTARG
		;;
	t)
		threshold=$OPTARG
		;;
        o)
                out_root=$OPTARG
                ;;
        c)
                cpu_num=$OPTARG
                ;;
        d)
                uniref90=$OPTARG
		#uniref90=`basename $uniref90 .fasta`
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
	exit 1
fi

if [ ! -f $uniref90 ]; then
	echo "ERROR: invalid sequence database file $uniref90"
	exit 1
fi

input_fasta=$1
if [ ! -f "$input_fasta" ]; then
	echo "ERROR: invalid sequence file $input_fasta" 
	exit 1
fi

fulnam=`basename $input_fasta`
relnam=${fulnam%.*}

if [ ! -d $out_root ]; then
	mkdir -p $out_root
fi

currDir=`pwd`

## create a temporary work directory. Do all the work in this temp folder
tmpWorkDir=`mktemp -d ${relnam}.tmpWork4Jackhmmer.XXXXXX`

cp $input_fasta $tmpWorkDir/$relnam.fasta
cd $tmpWorkDir

seq_file=$relnam.seq
$JackHome/util/Verify_FASTA $relnam.fasta $seq_file
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run $JackHome/util/Verify_FASTA $relnam.fasta $seq_file"
	exit 1
fi

tmp_root=$relnam"_"$$
mkdir -p $tmp_root

#-> generate config file for EVAlign
$JackHome/util/Gen_EVfold_config2 $tmp_root/$relnam $seq_file $cpu_num $bitscore $threshold $iteration $uniref90 $JackHome > $tmp_root/$relnam.config_file 

#-> run evcouplings_runcfg (align stage)
echo "evcouplings_runcfg start with parameter cpu $cpu_num, bitsco $bitscore, thres $threshold, iter $iteration, uniref $uniref90"
source activate evfold
evcouplings_runcfg $tmp_root/$relnam.config_file
conda deactivate
echo "evcouplings_runcfg done"

#-> generate a3m file
a3m_file=$relnam.a3m

if [ ! -f $tmp_root/$relnam/align/$relnam.a2m ]; then
	echo "ERROR: cannot find $tmp_root/$relnam/align/$relnam.a2m"
	exit 1
fi

$JackHome/util/EVfold_MSA_Trans $tmp_root/$relnam/align/$relnam.a2m > $a3m_file

cd $currDir

cp $input_fasta $out_root/$relnam.fasta_raw
mv $tmpWorkDir/$seq_file $out_root/
mv $tmpWorkDir/$a3m_file $out_root/

rm -rf $tmpWorkDir
