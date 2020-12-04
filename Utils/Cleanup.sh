#!/bin/bash

WorkDir1=Work4Server-raptorx/
WorkDir2=Work4Server-raptorx2/
WorkDir3=Work4Server-raptorx6/

#for WorkDir in $WorkDir1 $WorkDir2 $WorkDir3
for WorkDir in $WorkDir1 $WorkDir3
do

	if [ -d $WorkDir ]; then

		currDir=`pwd`
		cd $WorkDir
		seqIDs=`ls *.* | cut -f1 -d'.' | sort | sort -u -m `
		for id in $seqIDs
		do
			numLines=`ps -ef | grep RaptorX | grep $id | wc -l`
			if [ $numLines -le 1 ]; then
				echo "files for $id shall be removed"
				rm -rf ${id}.* ${id}_*
			fi
		done
		cd $currDir
	fi

done

WorkDir=Server/RaptorX_Server_2.11_result/
if [ -d $WorkDir ]; then
	find $WorkDir -type d -ctime +5 | xargs rm -rf
	find $WorkDir -type f -ctime +5 | xargs rm -rf
fi
