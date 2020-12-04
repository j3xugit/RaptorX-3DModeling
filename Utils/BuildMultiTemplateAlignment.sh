#!/bin/bash

## this is a wrapper for building sequence to multiple-template alignment 

RXHOME=/mnt/data/RaptorXCommon/RaptorX-Threading/

DeepAlignDir=$RXHOME/DeepAlign_Package
MRFAlignDir=$RXHOME/DeepThreader_Package
MultiTempDir=$RXHOME/MultiTemp_Package

PDBDIR=$RXHOME/databases/pdb_BC100
TPLDIR=$RXHOME/databases/TPL_BC100

function Usage {
	echo "$0 [ -p PDBDir | -t TPLDir ] tgtFile tplListFile aliFolder"
	echo "	This script builds a multi-template alignment for a query sequence"
	echo "	tgtFile: a .tgt file for the query sequence"
	echo "	tplListFile: a list of templates to be used, each in one orw"
	echo "	aliFolder: the folder for all pairwise query-to-template alignments"
	echo "	-p: the folder for structure files, only PDB format supported, default $PDBDIR"
	echo "	-t: the folder for template files, currently only .tpl file supported, default $TPLDIR"
	echo "	The resultant files will be saved to a folder with name XXX_MULTITEMP where XXX is the query protein name"
}
while getopts ":p:t:" opt; do
        case ${opt} in
                p )
                  PDBDIR=$OPTARG
                  ;;
                t )
                  TPLDIR=$OPTARG
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

if [ $# -lt 3 ]; then
	Usage
	exit 1
fi

tgtFile=$1
tempList=$2
aliFolder=$3

$MultiTempDir/MultiTemp_Proc.sh -q $tgtFile -i $tempList -f $aliFolder -p $PDBDIR -l $TPLDIR -D $DeepAlignDir -M $MRFAlignDir -H $MultiTempDir
