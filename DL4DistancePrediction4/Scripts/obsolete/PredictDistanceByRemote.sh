#!/bin/sh

if [ $# -lt 2 ]; then
	echo $0 "seqFile inputFeature_PKL [remoteAccount] [gpu]"
	echo "   this script predicts distance using a GPU in a remote machine"
	echo "   seqFile and inputFeature: the input sequence file in FASTA format and the input feature file in PKL format"
	echo "   remoteAccount: an account in the remote machine, e.g., RaptorX@raptorx5.uchicago.edu (default)"
	echo "   gpu: -1 (default), 0, 1, 2, 3. -1 indicates automatically select a GPU."
	exit 1
fi

seqFile=$1
if [ ! -f $1 ]; then
	echo "the input sequence file does not exist: $1"
	exit 1
fi
target=`basename $seqFile | cut -f1 -d'.'`

inputFeature=$2
if [ ! -f $2 ]; then
	echo "the input feature file does not exist: $2"
	exit 1
fi

RemoteAccount="RaptorX@raptorx5.uchicago.edu"
if [ $# -ge 3 ]; then
	RemoteAccount=$3
fi

gpu=-1
if [ $# -ge 4 ]; then
	gpu=$4
fi

pid=$$
localMachine=`hostname | cut -f1 -d'.' `

RemoteWorkDir="tmpWorkDir4RemoteDistancePrediction-${target}-$localMachine-$pid"

ssh -o StrictHostKeyChecking=no $RemoteAccount "mkdir -p $RemoteWorkDir"
if [ $? -ne 0 ]; then
	echo "Failed to create $RemoteWorkDir in a remote account $RemoteAccount!"
	exit 1
fi

scp $inputFeature $RemoteAccount:$RemoteWorkDir/
if [ $? -ne 0 ]; then
	echo "Failed to scp input feature file to a remote account $RemoteAccount!"
	exit 1
fi

bname=`basename $inputFeature`
scp $seqFile $RemoteAccount:$RemoteWorkDir/${target}.seq
if [ $? -ne 0 ]; then
	echo "Failed to scp input seq file to a remote account $RemoteAccount!"
	exit 1
fi

ssh -o StrictHostKeyChecking=no $RemoteAccount "\$DL4DistancePredHome/Utils/PredictDistance4PKLInput.sh $RemoteWorkDir/${target}.seq $RemoteWorkDir/$bname $gpu "
if [ $? -ne 0 ]; then
	echo "Failed to remotely predict distance using a remote account $RemoteAccount!"
	exit 1
fi

dname=`dirname $inputFeature`
scp -r $RemoteAccount:$RemoteWorkDir/*Server* $dname/
if [ $? -ne 0 ]; then
	echo "Failed to copy prediction results from the remote account $RemoteAccount!"
	exit 1
fi

ssh -o StrictHostKeyChecking=no $RemoteAccount "rm -rf $RemoteWorkDir"
