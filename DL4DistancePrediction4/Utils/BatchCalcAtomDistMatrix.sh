#!/bin/bash

list=$1
TPLDIR=PDB25-7952-TPL
PDBDIR=PDB25-7952-PDB

for i in `cat $list `
do
	python ../Utils/CalcAtomDistMatrixFromTPLPDB.py $TPLDIR/$i.tpl $PDBDIR/$i.pdb
	mv $i.atomDistMatrix.pkl PDB25-7952-atomDistMatrix/
done

