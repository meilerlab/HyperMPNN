#!/bin/bash
export MPNN=/media/data/software/ProteinMPNN-main/

folder_with_pdbs=$1

output_dir="design_output/hyper"
mkdir -p $output_dir


path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
chains_to_design="B"
design_only_positions="1 2 3 4 5 6 7 8 9 10 11 12 13 14 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 50 53 57 61 62 64 66 67 68 69 70 71 72 73 74 75 76 78 79 83 101 102 103 108 114 115 116 117 118 119 120 121 122 123 126 127 129 130 131 132 133 134 136 137 138 140 141 144 145 148"

python $MPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains
python $MPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
python $MPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed

python $MPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 200 \
        --sampling_temp "0.001" \
        --save_probs 1 \
        --batch_size 1 \
        --seed 1 \
        --path_to_model_weights ../retrained_models \
        --model_name v48_020_epoch300_hyper
