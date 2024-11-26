#!/usr/bin/env python3
import torch
import pandas as pd
from tqdm import tqdm
from sklearn.model_selection import train_test_split
import os, sys


def get_files_from_path(path:str, filesuffix:str):
    '''Return all files (with path) in a given path containing the filesuffix.'''
    return [os.path.join(path, file) for file in os.listdir(path) if file.endswith(filesuffix)]


def _get_pt_dict(input_path:str):
    # first get all ids from the generated pt files with theirs sequences
    # we only have chain A, therefore this should be simple
    # example pair:
    #   AF-A0A0E3GSV3-F1-model_v4.pt (ids, chains, sequences)
    #   AF-A0A0E3GSV3-F1-model_v4_A.pt (dont need that, filter them out)
    pt_files = get_files_from_path(input_path, filesuffix='.pt')
    pt_files = [file for file in pt_files if 'xyz' not in torch.load(file).keys()]

    pt_dict = {}
    print('[INFO] Read pt files.')
    for pt_file in tqdm(pt_files):
        basename = pt_file.split('/')[-1].replace('.pt', '')
        torch_file = torch.load(pt_file)

        # get sequence and id from it, just 6 letters
        # only the first sequence
        seqs = [torch_file['seq'][i][0] for i in range(len(torch_file['seq']))]
        seq_hash = str(torch_file['id']) + '_id'
        
        chain_list = [f'{basename}_{letter}' for letter in torch_file['chains']]
        fasta_list = [f'{basename}_chain_{letter}' for letter in torch_file['chains']]
        
        # other properties
        date = torch_file['date']
        resolution = torch_file['resolution']
        
        for file, chain, seq in zip(fasta_list, chain_list, seqs):
            pt_dict.update({file : {'model' : chain, 'hash': seq_hash, 'seq': seq, 
                            'date': date, 'resolution': resolution}})

    return pt_dict


def _create_clusters(cluster_file:str, pt_dict:dict):
    cluster_df = pd.read_csv(filepath_or_buffer=cluster_file, sep='\t', names=['represent', 'member'])
    # display(cluster_df)

    cluster_names = cluster_df['represent'].unique()
    n_clusters = len(cluster_names)
    print(f'N clusters: {n_clusters}')

    # iterate over all representatives and get the cluster members
    counter_cluster = 1
    member_dict = {'CHAINID': [], 'DEPOSITION': [], 'RESOLUTION': [], 'HASH': [], 'CLUSTER': [], 'SEQUENCE': []}
    print('[INFO] Create sequence list.')
    for id_cluster in cluster_names:
        # store current representative and get all cluster members
        represent = id_cluster
        members = cluster_df.loc[cluster_df['represent'] == id_cluster, : ]['member']
        # for every member create new entry in list
        for member in members:
            member_dict['CHAINID'].append(pt_dict[member]['model'])
            member_dict['DEPOSITION'].append(pt_dict[member]['date'])
            member_dict['RESOLUTION'].append(pt_dict[member]['resolution'])
            member_dict['HASH'].append(pt_dict[member]['hash'])
            member_dict['CLUSTER'].append(counter_cluster)
            member_dict['SEQUENCE'].append(pt_dict[member]['seq'])
        
        # add counter + 1 for new cluster
        counter_cluster += 1
        
    df_member = pd.DataFrame(data=member_dict)
    if not os.path.exists('../data_files'):
        os.makedirs('../data_files/')
    df_member.to_csv(path_or_buf='../data_files/list_cluster.csv', index=False)
    return df_member



if __name__ == '__main__':
    inputpt = sys.argv[1]
    clusterfile = sys.argv[2]
    
    pt_dict = _get_pt_dict(input_path=inputpt)
    df_member = _create_clusters(cluster_file=clusterfile, pt_dict=pt_dict)
    
    X_train, X_test, = train_test_split(df_member['CLUSTER'], test_size=0.15, random_state=1)
    X_test, X_val, = train_test_split(X_test, test_size=0.6, random_state=1)

    print(f'Train: {X_train.shape[0]}, Valid: {X_val.shape[0]}, Test: {X_test.shape[0]}')
    X_train.to_csv('../data_files/train_cluster.txt', index=False, header=False)
    X_val.to_csv('../data_files/valid_cluster.txt', index=False, header=False)
    X_test.to_csv('../data_files/test_cluster.txt', index=False, header=False)
    
    print('[INFO] Done.')
        