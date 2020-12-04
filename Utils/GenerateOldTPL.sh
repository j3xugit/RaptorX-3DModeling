#!/bin/bash

## this script generates an old tpl file from an .a3m  or .tgt file and a pdb file
## the tpl file is named after proteinName.tpl

if [[ -z "$ModelingHome" ]]; then
        echo "ERROR: please set environmental variable ModelingHome to the instllation directory of RaptorX-3DModeling"
        exit 1
fi

if [[ -z "$DistFeatureHome" ]]; then
        echo "ERROR: please set environmental variable DistFeatureHome to the instllation directory of BuildFeatures"
        exit 1
fi

if [ $# -lt 2 ]; then
	echo $0 "MSAfile/TGTfile structFile [ResultDir]"
	echo "	This script generates an old TPL file from an MSA (or TGT) file and a structure file"
	echo "	MSAfile/TGTfile: a .a3m  or .tgt file"
	echo "	structFile: a protein structure file in PDB or mmCif format"
	echo "	ResultDir: the folder for result saving, default current work directory"
	echo "	The result file is named after proteinName.tpl"
	exit 1
fi

MSAfile=$1
if [ ! -f $MSAfile ]; then
	echo 'ERROR: the input MSA or TGT file does not exist: ' $MSAfile
	exit 1
fi

structFile=$2
if [ ! -f $structFile ]; then
	echo "ERROR: invalid protein structure file $structFile"
	exit 1
fi

ResultDir=`pwd`
if [ $# -ge 3 ]; then
	ResultDir=$3
	mkdir -p $ResultDir
fi

fulnam=`basename $MSAfile`
if [[ "$fulnam" == *.a3m ]]; then
	target=${fulnam%.a3m}
	$DistFeatureHome/GenTGTFromA3M.sh $MSAfile $ResultDir
	tgtFile=$ResultDir/${target}.tgt

elif [[ "$fulnam" == *.tgt ]]; then
	tgtFile=$MSAfile
	target=${fulnam%.tgt}
else
	echo "ERROR: unsupproted file suffix for input file: " $MSAfile
	exit 1
fi

if [[ "$target" == *_* ]]; then
	fullChainID=`echo $target | cut -f2 -d'_'`
else
	fullChainID=`echo $target | cut -c5-`
fi

if [[ -z "$fullChainID" ]]; then
	echo "ERROR: empty chain ID in $target"
	exit 1
fi

chainID=`echo $fullChainID | cut -c1 `

fulnam=`basename $structFile`
removePDBFile=0
if [[ "$fulnam" == *.cif ]]; then
	bname=${fulnam%.cif}
	pdbFile=$ResultDir/${bname}.pdb
	python $ModelingHome/Common/CIF2PDB.py $structFile $pdbFile
	removePDBFile=1

elif [[ "$fulnam" == *.pdb ]]; then
	pdbFile=$structFile
else
	echo "ERROR: unsupproted file suffix for protein structure file: " $structFile
	exit 1
fi	

cmd=`readlink -f $0 `
cmdDir=`dirname $cmd`

tplFile=$ResultDir/$target.tpl
$cmdDir/TPLGenerator/TGT_To_TPL2 -i $pdbFile -I $tgtFile -h $chainID -H $cmdDir/TPLGenerator/ -o $tplFile

if [ ! -f $tplFile ]; then
	echo ERROR: failed to generate $tplFile
	exit 1
fi

if [[ "$fullChainID" != "$chainID" ]]; then
	## replace chainID by fullChainID
	oldStr="Chain ID = $chainID"
	newStr="Chain ID = $fullChainID"
	sed -i "s/$oldStr/$newStr/g" $tplFile
fi

## remove temporary files
if [ -f ${target}.pdb_post ]; then
	rm -f ${target}.pdb_post
fi
if [ $removePDBFile -eq 1 ]; then
	rm -f $ResultDir/${target}.pdb
fi
