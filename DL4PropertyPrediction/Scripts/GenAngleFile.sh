#!/bin/bash

if [ $# -lt 1 ]; then
        echo $0 "chainName, e.g. 5xyig "
        exit -1
fi

protein=$1
PDBTool=../../Common/PDB_Tool_angles/PDB_Tool
PDBDir=pdb_BC100/
TPLDir=TPL_BC100/
PDBTool2=../Utils/PropertyGenerator/GenerateProperty4OneProtein.py

if [ ! -f $TPLDir/${protein}.tpl ]; then
	echo "failed to generate the angle file due to lack of tpl file"
	exit -1
fi

${PDBTool} -i ${PDBDir}/${protein}.pdb -F 8 -o ${protein}.pdb
mv ${protein}.ang ${protein}.raw_ang
python $PDBTool2 ${protein} $TPLDir/${protein}.tpl ${protein}.raw_ang  ${protein}.ang
chmod a+r ${protein}.ang

HOST=`hostname`
rsync ${protein}.ang RaptorX@${HOST}:/mnt/data/RaptorXCommon/TrainTestData/ProteinProperty_Project/Ang_BC100/${protein}.ang
rm -f ${protein}.raw_ang ${protein}.pdb
