#!/bin/bash

python3 consensus_seq_pmpnn.py \
	--sort_by global_score \
	--output design_output/default_consensus.txt \
	design_output/default/seqs/6p6f.pentamer.renumbered.fa

python3 consensus_seq_pmpnn.py \
	--sort_by global_score \
	--output design_output/hyper_consensus.txt \
	design_output/hyper/seqs/6p6f.pentamer.renumbered.fa
