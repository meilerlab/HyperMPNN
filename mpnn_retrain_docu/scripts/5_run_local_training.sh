#!/bin/bash

python3 ../training.py \
    --path_for_training_data ../pdb_pt \
    --path_for_data_files ../data_files \
    --path_for_outputs ../train_output \
    --num_neighbors 48 \
    --backbone_noise 0.20 \
    --num_epochs 300 \
    --batch_size 3000 \
    --max_protein_length 5000 \
    --num_examples_per_epoch 1000 \
    --rerun
