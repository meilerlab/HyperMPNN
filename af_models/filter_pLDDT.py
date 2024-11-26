import argparse
import os
import numpy as np
import pandas as pd
import pyrosetta
import pickle
from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor
from typing import List, Tuple
from tqdm import tqdm
from math import pi

def remove_low_pLDDT_stretches(pose: pyrosetta.Pose, 
                               per_res_pLDDT: List[float], 
                               low_pLDDT_threshold: float) -> Tuple[pyrosetta.Pose, bool]:
    
    start, end = 1, len(per_res_pLDDT)
    modified = False

    # Calculate length of start and end stretches
    start_stretch_length = sum(1 for plddt in per_res_pLDDT if plddt < low_pLDDT_threshold)
    end_stretch_length = sum(1 for plddt in reversed(per_res_pLDDT) if plddt < low_pLDDT_threshold)

    # Skip deletion if stretches are more than half of the pose
    if start_stretch_length + end_stretch_length > len(per_res_pLDDT) / 2:
        return pose, False

    # Identifying and removing start stretch
    while start <= len(per_res_pLDDT) and per_res_pLDDT[start - 1] < low_pLDDT_threshold:
        start += 1
    if start - 1 >= 10:
        pose.delete_residue_range_slow(1, start - 1)
        modified = True
        # Adjust per_res_pLDDT list after removing start stretch
        per_res_pLDDT = per_res_pLDDT[start - 1:]

    # Recalculate the end index based on the updated per_res_pLDDT list
    end = len(per_res_pLDDT)

    # Identifying and removing end stretch
    while end > 0 and per_res_pLDDT[end - 1] < low_pLDDT_threshold:
        end -= 1
    if len(per_res_pLDDT) - end >= 10:
        pose.delete_residue_range_slow(end + 1, len(per_res_pLDDT))
        modified = True

    return pose, modified

def process_pdb_file(pdb_file: str, output_dir: str, pLDDT_threshold: float = 70) -> str:

    pose = pyrosetta.pose_from_pdb(pdb_file)
    num_residues = pose.total_residue()

    per_res_pLDDT = [pose.pdb_info().bfactor(i, 1) for i in range(1, num_residues + 1) if pose.residue(i).is_protein()]
    modified_pose, modified = remove_low_pLDDT_stretches(pose, per_res_pLDDT, pLDDT_threshold)
    
    total_pLDDT = sum([modified_pose.pdb_info().bfactor(i, 1) for i in range(1, modified_pose.total_residue() + 1)])
    new_average_pLDDT = total_pLDDT / modified_pose.total_residue() if modified_pose.total_residue() > 0 else 0
    
    if new_average_pLDDT >= pLDDT_threshold:
        pyrosetta.dump_pdb(modified_pose, os.path.join(output_dir, os.path.basename(pdb_file)))
        return pdb_file
    else:
        return None
        
        
def process_file_wrapper(pdb_file, output_dir, pLDDT_threshold):
    pyrosetta.init('-mute all')
    return process_pdb_file(pdb_file, output_dir, pLDDT_threshold)


        
def filter_based_on_pLDDT(pdb_files: List[str], output_dir: str, n_cores: int, pLDDT_threshold: float = 70) -> List[str]:
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        # Map the process_file_wrapper function to each pdb file
        results = list(tqdm(executor.map(process_file_wrapper, pdb_files, [output_dir]*len(pdb_files), [pLDDT_threshold]*len(pdb_files)), total=len(pdb_files)))

    filtered_pdb_files = [result for result in results if result is not None]
    return filtered_pdb_files

    
def main():
    parser = argparse.ArgumentParser(description="Filter out PDBs with low pLDDT (after removing N/C termini stretches of low pLDDT) and save the others to a specified folder.")
    parser.add_argument("input_folder", help="Folder with PDBs to filter")
    parser.add_argument("output_folder", help="Output folder to save resulting PDBs to")
    parser.add_argument("n_cores", help="The number of cores to use")
   
    args = parser.parse_args()

    pdbs = os.listdir(args.input_folder)
    pdbs = [args.input_folder+'/' + file for file in pdbs]
    filtered_pdbs = filter_based_on_pLDDT(pdbs, args.output_folder+'/', int(args.n_cores), pLDDT_threshold=70)


if __name__ == "__main__":
    main()
    
    

