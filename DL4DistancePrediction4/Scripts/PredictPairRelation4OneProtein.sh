#!/bin/bash

if [[ -z "${DL4DistancePredHome}" ]]; then
        echo "ERROR: please set environmental variable DL4DistancePredHome to the install folder of DL4DistancePrediction4"
        exit 1
fi
if [[ -z "${CUDA_ROOT}" ]]; then
        echo "ERROR: please set environmental variable CUDA_ROOT"
        exit 1
fi

#DeepModelFile=/mnt/data/RaptorXCommon/TrainTestData/Distance_Contact_TrainData/Jinbo_Folder/result4HBBeta/DistanceV3ModelFiles.txt
DeepModelFile=$DL4DistancePredHome/params/ModelFile4PairwisePred.txt
#ModelName=EC47C31C16CL99LargerS35V2020CbCbTwoRModels
#ModelName=EC47C37C19CL99S35V2020MidModels
DefaultModel4FM=EC47C37C19CL99S35V2020MidModels
DefaultModel4HHP=HHEC47C37C19CL99S35PDB70Models
DefaultModel4NDT=NDTEC47C37C19CL99S35BC40Models

ModelName=""

GPU=-1
ResultDir=`pwd`
UseMetaGenomeData=true
MSAmethod=4

tplStr=""
aliStr=""

alignmentType=0

Usage () 
{
        echo $0 "[ -f DeepModelFile | -m ModelName | -d ResultDir | -g gpu | -s MSA_mode | -M ] proteinName inputFolder"
        echo Or $0 "[ -f DeepModelFile | -m ModelName | -d ResultDir | -g gpu | -s MSA_mode | -M | -T alignmentType ] proteinName inputFolder aliFile/aliFolders tplFile/tplFolder"
        echo "	This script predicts distance/orientation for one protein using information in inputFolder"
	echo "	inputFolder: a folder generated by BuildFeatures.sh (e.g., T0955_OUT/T0955_contact/) containing subfolders feat_proteinName_XXX where XXX is an MSA generation method (e.g., uce3, ure3_meta)"
	echo "	-s: indicates which MSAs to be used, 0 for MSAs generated by hhblits, 1 for jackhmmer, 2 for both, 3 for user-provided MSA and 4 for all three, default 4"
	echo "	-M: if specified, do not use MSAs derived from metagenome data, otherwise use it when available"
	echo " "
	echo "	aliFile/aliFolders: optional, specify a query-template alignment file or one/multiple folders for many query-template alignment files"
        echo "		an alignment file shall be in FASTA format and has name proteinName-*.fasta"
	echo "		When mutliple folders are used, they shall be separated by semicolon without whitespace, e.g., Folder1;Folder2;Folder3"
	echo "		Two different alignment files shall have differnt names even if they are in different folders"
        echo "	tplFile/tplFolder: optional, specify a template file or a folder containing template files. A template file shall end with .tpl.pkl and may be generated by Common/MSA2TPL.sh"
        echo "  -T: indicate how query-template alignments are generated: 1 for alignments generated by HHpred and 2 for alignments generated by RaptorX threading"
        echo "          This option will be used only if both aliFile/aliFolders and tplFile/tplFolder are present"
	echo " "
	echo "  -f: a file containing path info for deep learning models, default $DeepModelFile"
        echo "  -m: a model name defined in DeepModelFile representing a set of deep learning models. Below is the default setting:"
	echo "          When aliFile/aliFolders are not used, $DefaultModel4FM will be used. Otherwise, when alignmentType is not set, $DefaultModel4HHP will be used"
        echo "          When aliFile/aliFolders are used, if alignmentType=2, $DefaultModel4NDT will be used; otherwise $DefaultModel4HHP will be used"
        echo "	-d: the folder for result saving, default current work directory"
        echo "	-g: -1 (default), 0-3. If -1, select a GPU automatically"
}

while getopts ":f:m:d:g:s:T:M" opt; do
        case ${opt} in
                f )
                  DeepModelFile=$OPTARG
                  ;;
                m )
                  ModelName=$OPTARG
                  ;;
                d )
                  ResultDir=$OPTARG
                  ;;
                g )
                  GPU=$OPTARG
                  ;;
		T )
		  alignmentType=$OPTARG
		  ;;
		M )
		  UseMetaGenomeData=false
		  ;;
		s )
		  MSAmethod=$OPTARG
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

if [ $# -ne 2 -a $# -ne 4 ]; then
        Usage
        exit 1
fi

proteinName=$1

metaFolder=$2
if [ ! -d $metaFolder ]; then
	echo "ERROR: invalid input folder for distance/orientation prediction: $metaFolder"
	exit 1
fi

if [ $# -eq 4 ]; then
	aliStr=$3
	tplStr=$4
fi

if [ ! -f $DeepModelFile ]; then
        echo "ERROR: invalid file for deep model path information: $DeepModelFile"
        exit 1
fi

if [ $MSAmethod -eq 4 ]; then
	methods="uce3 uce5 ure3 ure5 user"
elif [ $MSAmethod -eq 2 ]; then
	methods="uce3 uce5 ure3 ure5"
elif [ $MSAmethod -eq 0 ]; then
	methods="uce3 uce5"
elif [ $MSAmethod -eq 1 ]; then
	methods="ure3 ure5"
elif [ $MSAmethod -eq 3 ]; then
	methods="user"
else
	echo "ERROR: please specify which MSAs shall be used, hhblits MSA or jackhmmer MSA ?"
	exit 1
fi

for method in `echo $methods`
do
	#echo "MSA method=$method"
	if $UseMetaGenomeData; then
		featureFolder=$metaFolder/feat_${proteinName}_${method}_meta
		if [ ! -d $featureFolder ]; then
			featureFolder=$metaFolder/feat_${proteinName}_${method}
		fi
	else
		featureFolder=$metaFolder/feat_${proteinName}_${method}
	fi

	#echo $featureFolder
	if [ ! -d $featureFolder ]; then
		echo "WARNING: feature folder $featureFolder does not exist"
		continue
	fi
	if [ -z $inputFolders ]; then
		inputFolders=$featureFolder
	else
		inputFolders=$inputFolders';'$featureFolder
	fi
done

if [ -z $inputFolders ]; then
	echo "ERROR: there are no feature folders for $proteinName in $metaFolder"
	exit 1
fi
#echo $inputFolders

## load model file names
. $DeepModelFile

if [ $# -eq 2 ]; then
	if [ -z "$ModelName" ]; then
        	ModelName=$DefaultModel4FM
	fi
elif [ $# -eq 4 ]; then
	if [ -z "$ModelName" ]; then
        	if [ $alignmentType -eq 2 ]; then
                	ModelName=$DefaultModel4NDT
		else
                	ModelName=$DefaultModel4HHP
        	fi
	fi
fi

ModelFiles=`eval echo '$'${ModelName}`
#echo ModelFiles=$ModelFiles
if [ $ModelFiles == "" ]; then
	echo "ERROR: ModelFiles for $ModelName is empty"
	exit 1
fi

program=$DL4DistancePredHome/RunPairwisePredictor.py
if [ ! -f $program ]; then
	echo ERROR: invalid program $program
	exit 1
fi

if [ ! -d $ResultDir ]; then
	mkdir -p $ResultDir
fi

command=" python $program -m $ModelFiles -p $proteinName -i $inputFolders -d $ResultDir "
if [ ! -z "$aliStr" -a ! -z "$tplStr" ]; then
	command=$command" -a $aliStr -t $tplStr "
fi

if [ $GPU == "-1" ]; then
        ## here we assume 6G is sufficient, although sometimes not
        neededRAM=6403741824
        GPU=`$ModelingHome/Utils/FindOneGPUByMemory.sh $neededRAM 40`
fi

if [ $GPU == "-1" ]; then
        echo "WARNING: cannot find an appropriate GPU to run distance/orientation prediction!"
        exit 1
else
        GPU=cuda$GPU
fi

THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.include_path=${CUDA_ROOT}/include,dnn.library_path=${CUDA_ROOT}/lib64 $command
if [ $? -ne 0 ]; then
        echo "ERROR: failed to run $command"
        exit 1
fi
