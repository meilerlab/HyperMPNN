#!/bin/bash

input=$1
id=$(echo "${input}" | cut -d ',' -f1)
outname=$(echo "${input}" | cut -d ',' -f2)

mkdir -p ../pdb_pt

input="../input_pdb/${outname}.pdb"
echo "./parse_pdb_noX.py ${input} ../pdb_pt/${outname} ${id}"

python3 -W ignore ./parse_pdb_noX.py $input ../pdb_pt/${outname} "${id}"
