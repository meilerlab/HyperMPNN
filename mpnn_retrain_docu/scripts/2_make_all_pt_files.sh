#!/bin/bash

n_splits=2000
n_process=20

mkdir -p ../chunk
split -l $n_splits ../data_files/id_map.txt ../chunk/chunk_
# ONLY RUN ONES !!!
ls ../chunk/chunk_* | xargs -n 1 -P $n_process bash -c 'while read file; do bash helper_create_pt_files.sh $file; done < $1' _
