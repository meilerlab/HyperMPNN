#!/bin/bash

IDS_FILE=$1
OUTPUT_DIR=$2

mkdir -p "$OUTPUT_DIR"

while IFS= read -r uniprot_id
do
    PDB_FILE="${OUTPUT_DIR}/AF-${uniprot_id}-F1-model_v4.pdb"

    if [[ -f "$PDB_FILE" ]]; then
        echo "File for UniProt ID $uniprot_id already exists, skipping download."
    else
        echo "Downloading structure for UniProt ID: $uniprot_id"
        wget -q -P "$OUTPUT_DIR" "https://alphafold.ebi.ac.uk/files/AF-${uniprot_id}-F1-model_v4.pdb"
    fi
done < "$IDS_FILE"

echo "Download completed. Files saved in $OUTPUT_DIR."

