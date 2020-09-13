#!/bin/bash

if [[ -z "${ModelingHome}" ]]; then
        echo "ERROR: Please set environmental variable ModelingHome to the installation folder of the RaptorX-3DModeling package"
        exit 1
fi
if [[ -z "${DistFeatureHome}" ]]; then
        echo "ERROR: Please set environmental variable DistFeatureHome to the installation folder of BuildFeatures"
        exit 1
fi
if [[ -z "${DL4DistancePredHome}" ]]; then
        echo "ERROR: Please set environmental variable DL4DistancePredHome to the installation folder of DL4DistancePrediction4"
        exit 1
fi

#GPUMachineFile=$DistFeatureHome/params/GPUMachines.txt
GPUMachineFile=$ModelingHome/params/GPUMachines.txt
GPUmode=4

out_root=`pwd`
gpu=-1

function Usage
{
	echo $0 "[-o out_root | -g gpu | -r machineMode | -h MachineFile ] input_A3M "
	echo "	This script generates input features for distance/orientation prediction from an MSA file in a3m format"
	echo "	input_A3M: the MSA file in a3m format, in which the first sequence shall be query and not contain any gaps."
	echo "	out_root: the output folder for result saving, default current work directory"	    
	echo "	gpu: -1 (default), 0, 1, 2, 3. If -1, automatically choose one GPU"
 	echo "	-r: specifiy what kind of CPUs and GPUs to use, default $GPUmode"
	echo "	    1: use local GPUs if available"
        echo "	    2: use local GPUs and CPUs"
        echo "	    3: use local GPUs and GPUs of machines defined by the -h option"
        echo "	    4: use local GPUs/CPUs and GPUs of machines defined by the -h option"
        echo "	-h: a file specifying remote machines with GPUs, default $GPUMachineFile"
	echo "		Remote GPUs will not be used if this file does not exist"
}

#-> parse arguments
while getopts ":o:g:r:h:" opt; do
	case $opt in
	o)
		out_root=$OPTARG
		;;
	g)
		gpu=$OPTARG
		;;
	h )
                GPUMachineFile=$OPTARG
                ;;
        r )
                GPUmode=$OPTARG
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

if [ $# -ne 1 ]; then
	Usage
	exit 1
fi

input_file=`readlink -f $1`
if [ ! -f "$input_file" ]; then
	echo "ERROR: invalid input MSA file $input_file " >&2
	exit 1
fi

## get target name
fulnam=`basename $input_file`
relnam=${fulnam%.*}

#-> create out_root
mkdir -p $out_root
out_root=`readlink -f $out_root`

cmd=`readlink -f $0 `
cmdDir=`dirname $cmd`
CCM=${cmdDir}/Scripts/RunCCMpredWrapper.sh

if [ ! -f ${out_root}/${relnam}.ccmpred -o ! -f ${out_root}/${relnam}.extraCCM.pkl ]; then
	$CCM -d $out_root -g $gpu -r $GPUmode -h $GPUMachineFile $input_file
	if [ $? -ne 0 ]; then
		echo "ERROR: failed to $CCM -d $out_root -g $gpu -r $GPUmode -h $GPUMachineFile $input_file"
        	exit 1
	fi
fi

#-> bin folder
BinDir=$DistFeatureHome/bin/

currDir=`pwd`
pid=$$
workDir=$(mktemp -d -p `pwd` -t ${relnam}-GenDistFeature-${pid}-XXXXXX)
mkdir -p $workDir
cd $workDir

# ---- copy A3M file to output folder ----------#
a3m_file=$relnam.a3m
if [ ! -f "$out_root/$a3m_file" ]; then
	cp $input_file $out_root/$a3m_file
	#echo "a3m file copied"
fi

# ---- generate A2M file from A3M file ---------#
a2m_file=$relnam.a2m
$DistFeatureHome/util/A3M_To_PSI $out_root/$a3m_file $a2m_file.tmp 
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run $DistFeatureHome/util/A3M_To_PSI $out_root/$a3m_file $a2m_file.tmp"
	exit 1
fi
grep -v "ss_pred\|ss_conf" $a2m_file.tmp | awk '{print substr($0,34,length($0)-32) }' > $out_root/$a2m_file
rm -f $a2m_file.tmp

# ---- generate FASTA file from A2M file -------#
fasta_file=$relnam.seq
echo ">$relnam" > $out_root/$fasta_file
head -n1 $out_root/$a2m_file >> $out_root/$fasta_file

# ---- from A2M file, generate potential file ------- #
pot_file=$relnam.pot
if [ ! -f "$out_root/$pot_file" ]; then
	#echo "alnstats_omp starts"
	numLines=`wc -l $out_root/$a2m_file | cut -f1 -d' '`
	if [ $numLines -gt 50000 ]; then
		echo "WARNING: calculating $out_root/$pot_file by sampling 50000 seqs from $out_root/$a2m_file"
		python $DistFeatureHome/Helpers/SampleA2MByNumber.py $out_root/$a2m_file 50000 $out_root/${a2m_file}.sampled
		cp $out_root/$a2m_file $out_root/${a2m_file}.original
		cp $out_root/${a2m_file}.sampled $out_root/$a2m_file
		$BinDir/alnstats_omp $out_root/${a2m_file} $relnam.ws1 $out_root/$pot_file
		ret=$?
		mv $out_root/${a2m_file}.original $out_root/$a2m_file
	else
		$BinDir/alnstats_omp $out_root/$a2m_file $relnam.ws1 $out_root/$pot_file
		ret=$?
	fi

	if [ $ret -ne 0 ]; then
		echo "ERROR: failed to run $BinDir/alnstats_omp $out_root/$a2m_file $relnam.ws1 $out_root/$pot_file"
		exit 1
	fi
	rm -f $relnam.ws1
	#echo "alnstats_omp done"
fi

# ---- generate TGT file from A3M file ---------#
tgt_file=$relnam.tgt
if [ ! -f "$out_root/$tgt_file" ]; then
	#echo "a3m_to_tgt start"
	tmp=$relnam"_contact/"
	mkdir -p $tmp
	numLines=`wc -l $out_root/$a3m_file | cut -f1 -d' '`
	if [ $numLines -gt 100000 ]; then
		echo "WARNING: converting A3M to TGT by sampling 50000 seqs from $out_root/$a3m_file"
		python $DistFeatureHome/Helpers/SampleA3MByNumber.py $out_root/$a3m_file 50000 $out_root/${a3m_file}.sampled
		cp $out_root/$a3m_file $out_root/${a3m_file}.original
		cp $out_root/${a3m_file}.sampled $out_root/$a3m_file
		$DistFeatureHome/util/A3M_To_TGT -i $out_root/$fasta_file -I $out_root/${a3m_file} -o $out_root/$relnam.tgt -t $tmp 1> $relnam.ws1 2> $relnam.ws2
		ret=$?
		mv $out_root/${a3m_file}.original $out_root/$a3m_file
	else
		$DistFeatureHome/util/A3M_To_TGT -i $out_root/$fasta_file -I $out_root/$a3m_file -o $out_root/$relnam.tgt -t $tmp 1> $relnam.ws1 2> $relnam.ws2
		ret=$?
	fi

	if [ $ret -ne 0 ]; then
		echo "ERROR: failed to run util/A3M_To_TGT -i $out_root/$fasta_file -I $out_root/$a3m_file -o $out_root/$relnam.tgt -t $tmp"
		exit 1
	fi
	rm -f $tmp/$relnam.*
	rm -f $relnam.ws1 
	rm -f $relnam.ws2
	rm -rf $tmp
	#echo "a3m_to_tgt done"
fi

# ---- generate SS3 and ACC---- #
ss3_file=$relnam.ss3
if [ ! -f "$out_root/$ss3_file" ]; then
        $BinDir/DeepCNF_SS_Con -t $out_root/$tgt_file -s 1 > $out_root/$ss3_file
        if [ $? -ne 0 ]; then
                echo "ERROR: failed to run $BinDir/DeepCNF_SS_Con -t $out_root/$tgt_file -s 1 > $out_root/$ss3_file"
                exit 1
        fi
fi

acc_file=$relnam.acc
if [ ! -f "$out_root/$acc_file" ]; then
        $BinDir/AcconPred $out_root/$tgt_file 1 > $out_root/$acc_file
        if [ $? -ne 0 ]; then
                echo "ERROR: failed to run $BinDir/AcconPred $out_root/$tgt_file 1 > $out_root/$acc_file"
                exit 1
        fi
fi

# ---- make a pseudo DISO file -----#
diso_file=$relnam.diso
if [ ! -f "$out_root/$relnam.diso" ]; then
        tail -n+2 $out_root/$acc_file | awk '{if(NF==8){print $1" "$2" . 0"}else{print $0}}'  > $out_root/$diso_file
fi

python $DL4DistancePredHome/ReadSingleInputFeature.py $relnam $out_root $out_root
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run python $DL4DistancePredHome/ReadSingleInputFeature.py $relnam $out_root $out_root"
	exit 1
fi

cd $currDir
rm -rf $workDir

