#!/bin/bash

if [[ -z "${ModelingHome}" ]]; then
	echo "Please set environmental variable ModelingHome to the installation folder of the RaptorX-3DModeling package"
	exit 1
fi

if [[ -z "${DistFeatureHome}" ]]; then
	echo "Please set environmental variable DistFeatureHome to the installation folder of BuildFeatures"
	exit 1
fi

if [[ -z "${DL4DistancePredHome}" ]]; then
	echo "Please set environmental variable DL4DistancePredHome to the installation folder of DL4DistancePrediction4"
	exit 1
fi

#GPUMachineFile=$DistFeatureHome/params/GPUMachines.txt
GPUMachineFile=$ModelingHome/params/GPUMachines.txt
GPUmode=4

out_root=`pwd`
gpu=-1

function Usage
{
	echo $0 "[-d out_root | -g gpu | -r machineMode | -h machinefile ] input_A3M"
	echo "	iput_a3m: an MSA file in a3m format, in which the first sequence shall be the query and not contain any gaps"
	echo "	out_root: the folder for result saving, default current work directory"	    
	echo "	gpu: -1 (default), 0, 1, 2, 3. If -1, select GPU automatically"
	echo "	-r: specifiy what kind of CPUs and GPUs to use, default $GPUmode"
	echo "	     1: use local GPUs if available"
	echo "	     2: use local GPUs and CPUs"
	echo "	     3: use local GPUs and GPUs of machines defined by the -h option"
	echo "	     4: use local GPUs/CPUs and GPUs of machines defined by the -h option"
	echo "	-h: a file specifying remote machines with GPUs, default $GPUMachineFile"
	echo "		Remote GPUs will not be used if this file does not exist"
}

while getopts ":d:g:h:r:" opt; do
        case ${opt} in
                d )
                  out_root=$OPTARG
                  ;;
                g )
                  gpu=$OPTARG
                  ;;
		h )
		  GPUMachineFile=$OPTARG
		  ;;
		r )
		  GPUmode=$OPTARG
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

cmd=`readlink -f $0 `
cmdDir=`dirname $cmd`

input_file=$1
if [ ! -f $input_file ]; then
	echo "ERROR: invalid input a3m file $input_file " 
	exit 1
fi
input_file=`readlink -f $input_file`
fulnam=`basename $input_file`
relnam=${fulnam%.*}

if [ ! -d $out_root ]; then
	mkdir -p $out_root
fi
out_root=`readlink -f $out_root`

a3m_file=$relnam.a3m
if [ ! -f "$out_root/$a3m_file" ]; then
	cp $input_file $out_root/$a3m_file
fi

txtoutfile=${out_root}/$relnam.ccmpred
mpkoutfile=${out_root}/$relnam.ccmpred.mpk

if [ ! -f $txtoutfile -o ! -f $mpkoutfile ]; then

	if [ ! -f $GPUMachineFile -o $GPUmode -le 2 ]; then
		UseRemoteMachine=false
	else
		machines=""
		neededRAM=`$cmdDir/EstimateGPURAM4CCMpred.sh ${input_file}`
		FifteenG=16106127360
		if [ $neededRAM -gt $FifteenG ]; then
			machines=`grep -w on $GPUMachineFile | grep LargeRAM | cut -f1 -d' ' `
		fi

		## if cannot find any machine with large RAM, then use machine with smaller RAM
		if [[ -z "${machines}" ]]; then
			machines=`grep -w on $GPUMachineFile | grep RAM | cut -f1 -d' ' `
		fi

		## check to see if local machine is in machines or not
		localMachine=`hostname `
		UseRemoteMachine=true
		for machine in $machines
		do
			if [[ $localMachine == "$machine" ]]; then
                       		UseRemoteMachine=false
                       		break;
               		fi
		done
		if [[ -z "${machines}" ]]; then
			UseRemoteMachine=false
		fi

		## randomly pick up one remote machine
		if $UseRemoteMachine; then
			remoteAccount=`echo $machines | sed "s/ /\n/g"  | shuf | head -1`
		fi
	fi

	if $UseRemoteMachine; then
		echo "WARNING: try to run CCMpred at $remoteAccount ..."
		$cmdDir/RunCCMpredRemote.sh -d $out_root -g $gpu -r $remoteAccount $input_file
		if [ $? -ne 0 ]; then
			echo "ERROR: failed to run $cmdDir/RunCCMpredRemote.sh -d $out_root -g $gpu -r $remoteAccount $input_file"
			#exit 1
		fi
	else
        	$cmdDir/RunCCMpred.sh -d $out_root -g $gpu $input_file
		if [ $? -ne 0 ]; then
			echo "ERROR: failed to run $cmdDir/RunCCMpred.sh -d $out_root -g $gpu $input_file on `hostname` "
			#exit 1
		fi
	fi
	
	NotUseCPU=false
	if [ $GPUmode -ne 2 -a $GPUmode -ne 4 ]; then
		NotUseCPU=true
	fi

	if [ ! -f $txtoutfile -o ! -f $mpkoutfile ]; then
		if $NotUseCPU; then
			echo "ERROR: failed to run CCMpred on $input_file on GPUs. You may use CPUs to run CCMpred."
			exit 1
		fi

		a2mfile=${out_root}/$relnam.a2m
		tmpa2mfile=$a2mfile.tmp
		$DistFeatureHome/util/A3M_To_PSI $input_file $tmpa2mfile
		grep -v "ss_pred\|ss_conf" $tmpa2mfile | awk '{print substr($0,34,length($0)-32) }' > $a2mfile
		rm -f $tmpa2mfile

		numCPUs=10
		host=`hostname`
		if [ "$host" == "raptorx10.uchicago.edu" ]; then
			numCPUs=20
		fi
		echo "Running CCMpred on $numCPUs CPUs..."
		$DistFeatureHome/CCMpred/bin/ccmpred -R -t $numCPUs -b $mpkoutfile $a2mfile $txtoutfile 
		if [ $? -ne 0 ]; then
			echo "ERROR: failed to run $DistFeatureHome/CCMpred/bin/ccmpred -R -t $numCPUs -b $mpkoutfile $a2mfile $txtoutfile "
			exit 1
		fi
	fi

	if [ ! -f $txtoutfile -o ! -f $mpkoutfile ]; then
		echo "ERROR: failed to run CCMpred on $input_file"
		exit 1
	fi
fi

#extraCCMfile=${out_root}/relnam.extraCCM.pkl
zscorefile=${out_root}/$relnam.ccmpred_zscore
## convert mpk file to pkl file
python $DL4DistancePredHome/CCMpredUtils.py $mpkoutfile ${out_root}/
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run python $DL4DistancePredHome/CCMpredUtils.py $mpkoutfile ${out_root}/"
        exit 1
fi
rm -f $mpkoutfile

python $DistFeatureHome/bin/normalize_ccmpred_sep.py $txtoutfile > $zscorefile
if [ $? -ne 0 ]; then
        echo "ERROR: failed to run python $bin/normalize_ccmpred_sep.py on $txtoutfile"
        exit 1
fi

#cd $currDir
#rm -rf $workDir
