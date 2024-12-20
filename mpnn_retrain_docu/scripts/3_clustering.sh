#!/bin/bash

mkdir -p ../cluster_data

# create fasta file
python3 -W ignore ../scripts/get_multifasta_from_pdb_path.py \
	-p ../input_pdb/ -c A -o ../cluster_data/input_pdb_seqs.fasta

# Clustering the fasta sequences
echo "Start clustering"
mmseqs easy-cluster ../cluster_data/input_pdb_seqs.fasta ../cluster_data/input_pdb_seqs ./tmp \
 	--min-seq-id 0.3 -c 0.8 --cov-mode 0
