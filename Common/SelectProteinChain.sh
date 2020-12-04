#!/bin/bash

if [ $# -lt 2 ]; then
	echo "$0 proteinChain PDBDir [ ResDir ]"
	echo "proteinChain: it has form XXXX_Y or XXXX_YY"
	exit 1
fi

proteinChain=$1

PDBDIR=$2
if [ ! -d $PDBDIR ]; then
	echo "ERROR: invalid folder for raw PDB files $PDBDIR"
	exit 1
fi

RESDIR=`pwd`
if [ $# -ge 3 ]; then
	RESDIR=$3
fi

if [ ! -d $RESDIR ]; then
	mkdir -p $RESDIR
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/SelectProteinChains.py

pdbid=`echo $proteinChain | cut -f1 -d'_' | tr [:upper:] [:lower:] `
chainName=`echo $proteinChain | cut -f2 -d'_' `
			
pdbfile=${PDBDIR}/${pdbid}.pdb
ciffile=${PDBDIR}/${pdbid}.cif
if [ -f $pdbfile ]; then
  	python $program ${pdbfile} $chainName $RESDIR &
elif [ -f $ciffile ]; then 
	python $program ${ciffile} $chainName $RESDIR &
else
	echo "ERROR: no structure file for chain $proteinChain"
fi

