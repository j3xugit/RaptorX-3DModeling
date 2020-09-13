#!/bin/sh

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

#********************************you may have to edit the following code to set ModelingHome to the install folder of the RaptorX-3DModeling package******************
export ModelingHome=`dirname $cmdDir`
#****************************************************************************************************************************************************************************

. $ModelingHome/raptorx-external.sh

export DistFeatureHome=$ModelingHome/BuildFeatures/
export DL4DistancePredHome=$ModelingHome/DL4DistancePrediction4/
export DL4PropertyPredHome=$ModelingHome/DL4PropertyPrediction/
export DistanceFoldingHome=$ModelingHome/Folding/

export PYTHONPATH=$ModelingHome:$PYTHONPATH
export PATH=$ModelingHome/bin:$PATH

outDir=`pwd`
GPU=-1
MSAmethod=9
runningMode=0
numDecoys=120

## information for remote account, if empty, then run folding on local machine
RemoteAccountInfo=""
machineType=0

## the maximum length of a protein to be folded. If you want to fold a larger protein, you may set a larger value for this parameter
maxLen2BeFolded=1050

function Usage
{
	echo $0 "[ -o outDir | -g gpu | -m MSAmethod | -n numDecoys | -r runningMode | -R remoteAccountInfo | -t machineType | -l maxLen2BeFolded ] inputFile"
	echo "	This script predicts angle/contact/distance/orientation of a protein and optionally folds it"
	echo "		ModelingHome=$ModelingHome"
	echo "		Please make sure that ModelingHome is correctly set to the install folder of the RaptorX-3DModeling package"
	echo "		Before running this script, some external software packages shall be installed in your own account. See README.md for details"
	echo " "
	echo "	inputFile: a protein primary sequence file in FASTA format (ending with .fasta or .seq) or a multiple sequence alignment file in a3m format (ending with .a3m)"
	echo "	-o: the folder for results, default current work directory, in which a subfolder target_OUT will be created where target is the protein name"
	echo "	-g: 0-3, and -1(default). If -1, automatically select GPUs with the maximum amount of free memory"
	echo "		A very large protein may need GPUs with more than 12G memory"
	echo "	-l: the maximum length of a protein for which a 3D model will be built, default $maxLen2BeFolded"
	echo "		If you want to fold a larger protein, you may set a larger limitation through this option"
	echo "		However, bear in mind that when relaxation is applied, it may take ~10 hours and ~10G CPU memory to build a single 3D model for such a large protein"
	echo " "
	echo "	-m: an integer indicating MSA generation methods defined by combining 1, 2, 4, 8 and 16, default $MSAmethod"
        echo "	    	1: run HHblits to generate MSA for protein local structure property (e.g., Phi/Psi) prediction"
	echo "		2: run HHblits 2.0 to generate MSA for contact/distance/orientation prediction (obsolete)"
        echo "	    	4: run Jackhmmer to generate MSA for contact/distance/orientation prediction (slow)"
	echo "		8: run HHblits 3.0 to generate MSA for contact/distance/orientation prediction (fast)"
        echo "	    	16: search MetaGenome database to enhance MSAs for contact/distance/orientation prediction (slow)"
	echo " "
	echo "	    	NOTE that when m=0, inputFile shall be an MSA in a3m format, i.e., no new MSA will be generated"
	echo " "
	echo "	-n: the number of decoys to be generated, default $numDecoys"
	echo "		When <=0, folding will not be done, only contact/distance/orientation will be predicted"
	echo "	-r: 0 for fold only and 1 for fold+relax, default $runningMode"
	echo "		The relaxation step is very time-consuming. When the protein is large (>600AAs), it may take a few hours to relax one initial 3D model on a single CPU."
	echo " "
	echo "	-R: run the folding module at a remote account specified by this option, e.g., raptorx@raptorx3.uchicago.edu:Work4Server/. By default, folding is done locally."
	echo "		RaptorX (at least the folding module) shall be installed and some environmental variables shall be set on the remote machine"
	echo "		Please make sure that you are able to scp/ssh/rsync to this remote account without password."
	echo "		This option allows you to use several machines to run a single job without manually copying data among machines."
	echo "	-t: the type of machine for folding jobs: 0 (default) for self-determination (not recommended since it may not work for your machines), 1 for a mulit-CPU Linux computer with GNU parallel installed,"
	echo "		2 for a slurm cluster with homogenous nodes, 3 for a slurm cluster with hetergenous nodes and 4 for a multi-CPU Linux computer without GNU parallel installed"
}

while getopts ":o:g:m:n:r:R:t:l:" opt; do
        case ${opt} in
                o )
                  outDir=$OPTARG
                  ;;
                g )
                  GPU=$OPTARG
                  ;;
                m )
                  MSAmethod=$OPTARG
                  ;;
		n )
		  numDecoys=$OPTARG
		  ;;
		r )
		  runningMode=$OPTARG
		  ;;
		R )
		  RemoteAccountInfo=$OPTARG
		  ;;
		t )
		  machineType=$OPTARG
		  ;;
		l )
		  maxLen2BeFolded=$OPTARG
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

inFile=$1
if [ ! -f $inFile ]; then
	echo "ERROR: invalid input file: $inFile "
	exit 1
fi

fullname=`basename $inFile`
target=${fullname%.*}

if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi

## buildFeatures always try to use a GPU with the maximum amount of free memory
$DistFeatureHome/BuildFeatures.sh -o $outDir -g $GPU -m $MSAmethod -r 4 $inFile
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run $DistFeatureHome/BuildFeatures.sh -o $outDir -g $GPU -m $MSAmethod -r 4 $inFile "
	exit 1
fi
seqFile=$outDir/${target}_OUT/${target}.seq
if [ ! -f $seqFile ]; then
	echo "ERROR: cannot find sequence file $seqFile"
	exit 1
fi

$DL4PropertyPredHome/Scripts/PredictProperty4Server.sh -g $GPU $target $outDir/${target}_OUT/ 
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run $DL4PropertyPredHome/Scripts/PredictProperty4Server.sh -g $GPU $target $outDir/${target}_OUT/"
	exit 1
fi
## the predicted property file is saved in $outDir/${target}_OUT/PropertyPred/
predPropertyFile=$outDir/${target}_OUT/PropertyPred/${target}.predictedProperties.pkl
if [ ! -f $predPropertyFile ]; then
	echo "ERROR: failed to predict structure properties such as Phi/Psi angles for $seqFile"
	exit 1
fi

$DL4DistancePredHome/Scripts/PredictPairRelation4Server.sh -g $GPU $target $outDir/${target}_OUT/
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run $DL4DistancePredHome/Scripts/PredictPairRelation4Server.sh -g $GPU $target $outDir/${target}_OUT/"
	exit 1
fi
## the predicted distance/orientation file is saved in $outDir/${target}_OUT/DistancePred/
predMatrixFile=$outDir/${target}_OUT/DistancePred/${target}.predictedDistMatrix.pkl
if [ ! -f $predMatrixFile ]; then
	echo "ERROR: failed to predict distance/orientation matrix for $seqFile"
	exit 1
fi

echo "Finished predicting angle/contact/distance/orientation for $target. All results are in $outDir/${target}_OUT/"
if [ $numDecoys -le 0 ]; then
	exit 0
fi

seqLen=`tail -n +1 $seqFile | wc -c`
if [ $seqLen -gt $maxLen2BeFolded ]; then
	echo "WARNING: the protein sequence in $seqFile is too long (more than $maxLen2BeFolded residues) and no 3D models will be built"
	echo "	If you want to build 3D models for such a large protein, please increase the limit through the -l option"
	exit 0
fi

## generate decoys
decoyFolder=$outDir/${target}_OUT/${target}-RelaxResults
mkdir -p $decoyFolder

if [[ -z "${RemoteAccountInfo}" ]]; then
	command=$DistanceFoldingHome/LocalFoldNRelaxOneTarget.sh
else
	command="$DistanceFoldingHome/RemoteFoldNRelaxOneTarget.sh -R $RemoteAccountInfo "
fi
$command -t $machineType -d $decoyFolder -n $numDecoys -r $runningMode $seqFile $predMatrixFile $predPropertyFile
if [ $? -ne 0 ]; then
        echo "ERROR: failed to run $command -d $decoyFolder -n $numDecoys -r $runningMode $seqFile $predMatrixFile $predPropertyFile"
        exit 1
fi

## Cluster all decoys
$DistanceFoldingHome/Scripts4SPICKER/SpickerOneTarget.sh -d $outDir/${target}_OUT/${target}-SpickerResults $seqFile $decoyFolder
if [ $? -ne 0 ]; then
        echo "ERROR: failed to run $DistanceFoldingHome/Scripts4SPICKER/SpickerOneTarget.sh -d $outDir/${target}_OUT/${target}-SpickerResults $seqFile $decoyFolder"
        exit 1
fi

echo "Finished folding $target. Please check out results in $outDir/${target}_OUT/"
