#!/bin/bash

if [[ -z "${ModelingHome}" ]]; then
        echo "ERROR: Please set environmental variable ModelingHome to the installation folder of the RaptorX-3DModeling package"
        exit 1
fi

#GPUMachineFile=$DistFeatureHome/params/GPUMachines.txt
GPUMachineFile=$ModelingHome/params/GPUMachines.txt
GPUmode=4

numAllowedJobs=3
ResDir=`pwd`
gpu=-1

function Usage 
{
	echo $0 "[ -o ResDir | -g gpu | -n numJobs | -r machineMode | -h MachineFile ] proteinListFile MSADir"
	echo "	This script generates input features for distance/orientation prediction from a set of MSAs, each for one protein"
	echo "	proteinListFile: a file for a list of protein names, each in one row"
	echo "	MSADir: a folder for MSAs, each in a3m format"
	echo " "
	echo "	-o: a folder for result saving, default current work directory"
	echo "	    for each protein, three feature files XXX.inputFeatures.pkl, XXX.extraCCM.pkl and XXX.a2m will be saved to a subfolder ResDir/feat_proteinName_contact/"
	echo "	    these feature files will also be linked to ResDir/ where XXX is protein name"
	echo "	-g: -1 (default), 0, 1, 2 and 3. If -1, choose one gpu automatically"
	echo "	-n: the number of proteins to be run simultaneously, default $numAllowedJobs"
	echo "	    Note that all these jobs may compete for the same GPU with a limited amount of memory, so please use a small number "
	echo "	    Set gpu=-1 will allow these jobs to use more than 1 GPUs if available"
	echo "	-r: specifiy what kind of CPUs and GPUs to use, default $GPUmode"
        echo "	     1: use local GPUs if available"
        echo "	     2: use local GPUs and CPUs"
        echo "	     3: use local GPUs and GPUs of machines defined by the -h option"
        echo "	     4: use local GPUs/CPUs and GPUs of machines defined by the -h option"
        echo "	-h: a file specifying remote machines with GPUs, default $GPUMachineFile"
	echo "          Remote GPUs will not be used if this file does not exist"
}

while getopts ":o:g:n:r:h:" opt; do
        case $opt in
        o)
                ResDir=$OPTARG
                ;;
        g)
                gpu=$OPTARG
                ;;
	n)
		numAllowedJobs=$OPTARG
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

if [ $# -ne 2 ]; then
	Usage
	exit 1
fi

proteinListFile=$1
if [ ! -f $proteinListFile ]; then
	echo "ERROR: invalid protein list file $proteinListFile"
	exit 1
fi
targets=`cat $proteinListFile`

MSADir=$2
if [ ! -d $MSADir ]; then
	echo "ERROR: invalid folder for MSAs: $MSADir"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/GenDistFeaturesFromMSA.sh 
if [ ! -x $program ]; then
	echo "ERROR: $program is not executable"
	exit 1
fi

keywords=`basename $program`
myself=`basename $0 `

if [ ! -d $ResDir ]; then
	mkdir -p $ResDir
fi

for i in $targets
do
	while true
	do
                numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l`
                #numRunningJobs=`pgrep ${keywords} | wc -l`
                if [ $numRunningJobs -le $numAllowedJobs  ]; then
			if [ -f $MSADir/${i}.a3m ]; then
				$program -o $ResDir/feat_${i}_contact -g $gpu -r $GPUmode -h $GPUMachineFile $MSADir/${i}.a3m &
				sleep 2
			fi
			break
                else
                        a=`expr $RANDOM % 4`
                        sleep $a
                fi
        done
done

wait

## link files
currDir=`pwd`

cd $ResDir
for i in $targets
do
	ln -s feat_${i}_contact/${i}.inputFeatures.pkl
	ln -s feat_${i}_contact/${i}.extraCCM.pkl
	ln -s feat_${i}_contact/${i}.a2m
done

cd $currDir
