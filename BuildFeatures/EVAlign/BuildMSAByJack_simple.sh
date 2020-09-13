#!/bin/bash

# ----- usage ------ #
usage
{
	echo "Version 3.00 [2018-11-20] "
	echo $0 " <-i input_fasta> [-n iter] [-b bitsco] [-t thres] [-o out_root] [-c CPU_num] [-d uniref90]"
	echo "[note1]: default JackHmmer iteration number is 3, CPU_num is 5 "
	echo "[note2]: thres (x) should between 0 and 1, and default is 0.001 for E-value"
	echo "         When BitScore is used, then set BitScore threshold in JackHmmer as seqlen*x ."
	echo "[note3]: bitsco is set to 0 for E-value by default. "
	echo "         If set to 1, then use BitScore and thres should be set to at least 0.5. "
	echo "[note4]: Default value of job_id is process ID, out_root is current directory. "
	echo "[note5]: Database to be searched (default $DistFeatureHome/EVAlign/databases/uniref90.fasta) "
	exit 1
}

if [ $# -lt 1 ];
then
        usage
fi


if [[ -z "${DistFeatureHome}" ]]; then
	echo "Please set the environmental variable DistFeatureHome, e.g, $HOME/3DModeling/BuildFeautres "
	exit 1
fi

JackHome=$DistFeatureHome/EVAlign

#-> optional arguments
out_root="./"   #-> output to current directory
cpu_num=5       #-> use 1 CPUs

#-> required arguments
input_fasta=""
iteration="3"
bitscore="0"
threshold="0.001"
uniref90=$JackHome/EVAlign/databases/uniref90.fasta

#-> parse arguments
while getopts ":i:n:b:t:o:c:d:" opt;
do
        case $opt in
        #-> required arguments
        i)
                input_fasta=$OPTARG
                ;;
        #-> optional arguments
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

##uniref90=${uniref90}.fasta


if [ ! -f "$input_fasta" ]
then
	echo "input sequence file $input_fasta not found !" >&2
	exit 1
fi

uniref90=`readlink -f $uniref90`

if [ ! -f $uniref90 ]; then
	echo "the sequence database does not exist: $uniref90"
	exit 1
fi

fulnam=`basename $input_fasta`
relnam=${fulnam%.*}

if [ ! -d $out_root ]; then
	mkdir -p $out_root
fi

currDir=`pwd `

## create a temporary work directory. Do all the work in this temp folder
tmpWorkDir=`mktemp -d ${relnam}.tmpWork4Jackhmmer.XXXXXX`

cp $input_fasta $tmpWorkDir/$relnam.fasta
cd $tmpWorkDir

# ---- verify FASTA file -------- #
seq_file=$relnam.seq
$JackHome/util/Verify_FASTA $relnam.fasta $seq_file
if [ $? -ne 0 ]
then
	echo "failed in running $JackHome/util/Verify_FASTA $relnam.fasta $seq_file"
	exit 1
fi


if [ "$bitscore" == "1" ]; then
	echo "Running jackhmmer with bitscore as cutoff..."
	$JackHome/util/jackhmmer -o /dev/null -A $relnam.sto -N $iteration --noali --cpu $cpu_num -T $threshold --domT $threshold --incT $threshold --incdomT $threshold ${seq_file} $uniref90
	ret=$?
else
	echo "Running jackhmmer with E-value as cutoff..."
	$JackHome/util/jackhmmer -o /dev/null -A $relnam.sto -N $iteration --noali --cpu $cpu_num -E $threshold --domE $threshold --incE $threshold --incdomE $threshold ${seq_file} $uniref90
	ret=$?
fi

if [ $ret -ne 0 ]
then
        echo "ERROR: Failed to run $JackHome/util/jackhmmer for ${seq_file} on database $uniref90"
	rm -rf $tmpWorkDir
        exit 1
fi


#-> generate a3m file
a3m_file=$relnam.a3m
$JackHome/util/reformat.pl $relnam.sto $relnam.a3m -M first
ret1=$?

rm -f $relnam.sto

seqLen=`tail -1 ${seq_file} | tr -d '\n' | wc -c`
echo seqLen=$seqLen
$JackHome/util/EVfold_MSA_Trans  $relnam.a3m | cut -c1-$seqLen >  $relnam.a3m.2
ret2=$?

#$JackHome/util/hhfilter -i $relnam.a3m.2 -o $relnam.a3m -id 95 -cov 50 -M first
$DistFeatureHome/MSA_filter/self_filter -i $relnam.a3m.2 -o $relnam.a3m
ret3=$?

rm -f $relnam.a3m.2

if [ $ret1 -ne 0 -o $ret2 -ne 0 -o $ret3 -ne 0 ]; then
	echo "Failed to reformat $relnam.sto to a3m format!"
	rm -f $relnam.*
	exit 1
fi

cd $currDir

cp $input_fasta $out_root/$relnam.fasta_raw
mv $tmpWorkDir/$seq_file $out_root/
mv $tmpWorkDir/$a3m_file $out_root/
rm -rf $tmpWorkDir

