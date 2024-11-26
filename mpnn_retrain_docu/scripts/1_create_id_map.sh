#!/bin/bash

input=../input_pdb/
mkdir -p ../data_files/
touch ../data_files/id_map.txt

COUNTER=0
base36="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
for item in ${input}*.pdb;
do
	COUNT=$COUNTER
	result=""
	while true; do
		result=${base36:((COUNT%36)):1}${result}
		if [ $((COUNT=${COUNT}/36)) -eq 0 ]; then
			break
		fi
	done
	printf -v padded "%04s" "$result"
	result="${padded// /0}"
	
	#python3 parse_pdb_noX.py $item 
	outname="${item/.pdb/}"
	outname="${outname##*/}"
	echo "${result},${outname}" >> ../data_files/id_map.txt
	
	COUNTER=$((COUNTER+1))
done
