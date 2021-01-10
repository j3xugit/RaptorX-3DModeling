#!/bin/sh

savefolder=''
selectModel=false
QA=false

function Usage 
{
        echo $0 "[ -s | -a | -d savefolder ] SeqFile ModelFolder1 ModelFolder2 ModelFolder3 ..."
	echo "	-s: when specified, select models by energy, default no"
	echo "	-a: when specified, assign model quality, default no"
	echo "	-d: the folder for result saving, default bname-SpickerResults/ in current work directory where bname is the base name of the first modelFolder"
	echo "	SeqFile: the primary seq file in FASTA format"
	echo "	ModelFolder: the folder for decoys to be clustered"
}

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
#cmdDir=$DistanceFoldingHome/Scripts4SPICKER

if [ "$selectModel" = true ]; then
	if [[ -z "$DistanceFoldingHome" ]]; then
		echo "ERROR: please set environmental variable DistanceFoldingHome to the install folder of Folding/"
		exit 1
	fi

fi

while getopts ":sad:" opt; do
       case ${opt} in
        s )
          selectModel=true
          ;;
	a )
	  QA=true
	  ;;
	d )
	  savefolder=$OPTARG
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

if [ $# -lt 2 ]; then
	Usage
	exit 1
fi

seqFile=$1
if [ ! -f $seqFile ]; then
	echo "ERROR: invalid sequence file $seqFile"
	exit 1
fi

program=$cmdDir/GenInputInfo4SPICKER.py
if [ ! -f $program ]; then
	echo "ERROR: invalid program $program"
	exit 1
fi

modelFolder=$2
if [ ! -d $modelFolder ]; then
	echo "ERROR: invalid folder for decoys to be clustered $modelFolder"
	exit 1
fi

if [ -z "$savefolder" ]; then
	bname=`basename $modelFolder -RelaxResults`
	savefolder=${bname}-SpickerResults
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

modelListFile=$savefolder/models4clustering.list
cp /dev/null $modelListFile

modelFolders=${@:2}
for modelFolder in $modelFolders
do
	if [ "$selectModel" = true ]; then
		bname=`basename $modelFolder`
		scorefile=$savefolder/${bname}.score
		cp /dev/null $scorefile
		$DistanceFoldingHome/Scripts4Rosetta/ExtractScoreFromRelaxedModels.sh $modelFolder $scorefile
		if [ $? -ne 0 ]; then
			echo "ERROR: failed to run $DistanceFoldingHome/Scripts4Rosetta/ExtractScoreFromRelaxedModels.sh $modelFolder $scorefile"
			exit 1
		fi
		python $DistanceFoldingHome/Scripts4Rosetta/SelectModels4Clustering.py $scorefile >> $modelListFile
		if [ $? -ne 0 ]; then
			echo "ERROR: failed to run python $DistanceFoldingHome/Scripts4Rosetta/SelectModels4Clustering.py $scorefile >> $modelListFile"
			exit 1
		fi
	else
		ls $modelFolder/*.pdb >> $modelListFile
	fi
done

python $program -l -s $savefolder $seqFile $modelListFile
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run python $program -l -s $savefolder $seqFile $modelListFile"
	exit 1
fi

currDir=`pwd`
cd $savefolder
spicker
cd $currDir

## parse clustering result
name=`basename $savefolder -SpickerResults`
python $cmdDir/ParseSpickerResult.py $name $savefolder
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run python $cmdDir/ParseSpickerResult.py $name $savefolder"
	exit 1
fi

if [ "$QA" == false ]; then	
	echo "Do not assign model quality to the selected 3D models in $savefolder"
	exit 0
fi

echo "Assigning model quality..."
if [[ -z "$ModelingHome" ]]; then
	echo "ERROR: please set environmental variable ModelingHome to the install folder of RaptorX-3DModeling"
	exit 1
fi

QualityAssessProgram=${ModelingHome}/Utils/AssessModelByRef.sh
for i in $savefolder/*_center?.pdb
do
	$QualityAssessProgram $i $modelFolders &
done

wait
