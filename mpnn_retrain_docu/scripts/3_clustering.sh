#!/bin/bash

# create fasta file
python3 -W ignore ../scripts/get_multifasta_from_pdb_path.py \
	-p ../pdb_trimmed/ -c A -o ../cluster_data/pdb_trimmed_seqs.fasta

# Clustering the fasta sequences
echo "Start clustering"
mmseqs easy-cluster ../cluster_data/pdb_trimmed_seqs.fasta ../cluster_data/pdb_trimmed_seqs ./tmp \
 	--min-seq-id 0.3 -c 0.8 --cov-mode 0
