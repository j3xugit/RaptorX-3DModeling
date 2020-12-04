#!/bin/bash

outDir=`pwd`
gpu=-1
RemoteAccount="raptorx7.uchicago.edu"

function Usage
{
	echo $0 "[ -d outDir | -g gpu | -r remoteAccount ] a3mFile"
	echo "	This script runs CCMpred on a GPU of a remote machine"
	echo "	a3mFile: a multiple sequence alignment file in a3m format"
	echo "	outDir: the folder for result saving, default current work directory"
	echo "	remoteAccount: an account in the remote machine, e.g., youraccount@raptorx7.uchicago.edu (default)"
	echo "	gpu: -1 (default), 0, 1, 2, 3. if -1, select a GPU automatically"
}

while getopts ":d:g:r:" opt; do
        case ${opt} in
                d )
                  outDir=$OPTARG
                  ;;
                g )
                  gpu=$OPTARG
                  ;;
                r )
                  RemoteAccount=$OPTARG
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

if [ $# -lt 1 ]; then
	Usage
	exit 1
fi

a3mfile=$1
if [ ! -f $a3mfile ]; then
	echo "ERROR: invalid input a3m file $a3mfile "
	exit 1
fi
bname=`basename $a3mfile `
target=${bname%.*}

localMachine=`hostname | cut -f1 -d'.' `

pid=$$
RemoteWorkDir="tmpWorkDir4RemoteCCMpred-${target}-$localMachine-$pid"
ssh -o StrictHostKeyChecking=no $RemoteAccount "mkdir -p $RemoteWorkDir"
if [ $? -ne 0 ]; then
	echo "ERROR: failed to create $RemoteWorkDir in a remote account $RemoteAccount!"
	exit 1
fi

scp $a3mfile $RemoteAccount:$RemoteWorkDir/
if [ $? -ne 0 ]; then
	echo "ERROR: failed to scp $a3mfile to a remote account $RemoteAccount!"
	exit 1
fi

ssh -o StrictHostKeyChecking=no $RemoteAccount "\$DistFeatureHome/Scripts/RunCCMpred.sh -d $RemoteWorkDir -g $gpu $RemoteWorkDir/$bname"
if [ $? -ne 0 ]; then
	echo "ERROR: failed to remotely run CCMPred on $a3mfile at $RemoteAccount!"
	ssh -o StrictHostKeyChecking=no $RemoteAccount "rm -rf $RemoteWorkDir"
	exit 1
fi

if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi

scp -r $RemoteAccount:$RemoteWorkDir/${target}.ccmpred $outDir/
ret1=$?
scp -r $RemoteAccount:$RemoteWorkDir/${target}.ccmpred.mpk $outDir/
ret2=$?
scp -r $RemoteAccount:$RemoteWorkDir/${target}.a2m $outDir/
ret3=$?

if [ $ret1 -ne 0 -o $ret2 -ne 0 -o $ret3 -ne 0 ]; then
	echo "ERROR: failed to scp CCMpred results from $RemoteAccount to local $outDir"
	ssh -o StrictHostKeyChecking=no $RemoteAccount "rm -rf $RemoteWorkDir"
	exit 1
fi

ssh -o StrictHostKeyChecking=no $RemoteAccount "rm -rf $RemoteWorkDir"
