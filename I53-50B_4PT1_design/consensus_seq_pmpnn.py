import os
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Get consensus seqeunce from MPNN sequence sampling run. Best sequences ' + 
                                     'are selected by sorting local or global scores. Consensus sequence is calculated ' +
                                     'by finding the most common base at each position in the top N sequences.\n\n' +
                                     'Outfile contains the native and consensus sequence, a match case indicating ' +
                                     'differences between both and the mutations between the native and consensus sequence.')
    parser.add_argument('fasta_file', type=str, help='Path to the fasta file containing the sequences from a ProteinMPNN run.')
    parser.add_argument('--top_p', type=float, help='Top %% of sequences to consider for consensus. Default: 10(%%)', default=10.0)
    parser.add_argument('--sort_by', type=str, help='Sort sequences by "global_score" or "local_score" for consensus. Default: local_score', default='local_score')
    parser.add_argument('--output', type=str, help='Path to the output file to save the consensus sequence. Default: consensus.txt', default='consensus.txt')
    parser.add_argument('--print_consensus', action='store_true', help='Print consensus sequence on the console (maybe good for parsing). Default: False', default=False)
    parser.add_argument('--print_mutations', action='store_true', help='Print mutations between native and consensus sequence on the console (maybe good for parsing). Default: False', default=False)
    parser.add_argument('-v', '--verbose', action='store_true', help='Print informations on the console. Default: False', default=False)
    args = parser.parse_args()
    return args


def parse_mpnn(file_path:str, start:int=0, end:int=None):
    '''
    Parse a ProteinMPNN fasta file and return a dictionary of sequences and their corresponding features.
    Can be parsed into a pandas DataFrame for further analysis. Extracting only sequence fraction is also
    possible.
    Args:
        file_path (str):    Path to the fasta file.
        start (int):        Start position of the sequence. Default is 0.
        end (int):          End position of the sequence. Default is None.
    Returns:
        com (dict):         Dictionary of sequences and their corresponding features.
    '''
    def parse_line(line:str, start:int=0, end:int=None):
        if '>' in line:
            if line.startswith('>T='):
                elements = [float(c.strip().split('=')[-1]) for c in line.split(',')]
                return {'T': elements[0], 'local_score': elements[2], 'global_score': elements[3], 'seq_rec': elements[-1]}
            return {'T': 0, 'local_score': 0, 'global_score': 0, 'seq_rec': 1}
        else:
            return line[start:end].replace('\n', '')
    
    dic = {}
    com = {'T': [], 'local_score': [], 'global_score': [], 'seq_rec': [], 'seq': []}
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            content = parse_line(line, start, end)
            
            if isinstance(content, dict):
                dic = content
            elif isinstance(content, str):
                dic['seq'] = content
            else:
                continue
            
            if 'seq' in dic.keys() and 'T' in dic.keys():
                for key in com.keys():
                    com[key].append(dic[key])
    
    if len(com['T']) == 0:
        raise ValueError(f'No sequences found in the file "{file_path}" for the ProteinMPNN format. Return none.')
    return com


def extract_consensus_top_sequences(fasta_file:str, top:float=10, feature:str='local_score', verbose:bool=False):
    '''
    Extract consensus sequence from the top N sequences in a ProteinMPNN fasta file, by sorting the sequences 
    based on a feature (T, local_score, global_score, seq_rec) and calculating the most common base at each position.
    Args:
        fasta_file (str):   Path to the fasta file.
        top (float):        Top % of sequences to be used for consensus sequence calculation. Default is 10(%).
        feature (str):      Feature to be used for sorting the sequences. Default is 'local_score'.
        verbose (bool):     Print out the extracted sequences and consensus sequence. Default is False.
    Returns:
        native_seq (str):   Native sequence of the ProteinMPNN fasta file.
        consensus (str):    Consensus sequence of the top N sequences.
        match_case (str):   Match case of the native and consensus sequence.
        mutations (list):   List of mutations between the native and consensus sequence.
    '''
    if not os.path.exists(fasta_file):
        print(f'File "{fasta_file}" not found. Return none.')
        return None
    if feature not in ['local_score', 'global_score', 'seq_rec']:
        print(f'Feature "{feature}" not found. Must be one of the following: local_score, global_score, seq_rec. Return none.')
        exit(1)
    
    # Get all sequences as table to analyse it easier with pandas
    # by sorting for a specifc feature
    try:
        df_seq = pd.DataFrame(parse_mpnn(fasta_file))
    except:
        print(f'Error while parsing the fasta file "{fasta_file}". Exit.')
        return None
    native_seq = df_seq['seq'][0]
    df_seq.sort_values(feature, inplace=True, ascending=True)
    
    top_n = int(len(df_seq) * (top/100))
    if top_n == 0:
        top_n = 1
    elif top_n > len(df_seq):
        top_n = len(df_seq)
    elif top_n < 0:
        print(f'Invalid top value specified: "{top}" percent. Return none.')
        return None

    top_seqs = df_seq['seq'][1:top_n+1]
    if verbose and not args.print_consensus: 
        print(f'Extracted top {top_n} sequences from a total of {df_seq.shape[0]-1}.')
    # Calculate consensus sequence from the top N sequences
    # Find most common base at each position
    consensus = ""
    for i in range(min(len(seq) for seq in top_seqs)):
        column = [seq[i] for seq in top_seqs]
        most_common = max(set(column), key=column.count)  # Find most common base
        consensus += most_common
    
    match_case = ''.join(['-' if native_seq[i] == consensus[i] else '*' for i in range(len(native_seq))])
    mutations = [f'{native_seq[i]}{i+1}{consensus[i]}' for i in range(len(native_seq)) if native_seq[i] != consensus[i] ]
    
    if verbose:
        print(f'Native sequence:\t{native_seq}')
        print(f'Consensus sequence:\t{consensus}')
        print(f'Match case:\t\t{match_case}')
        print(f'Mutations:\t\t{", ".join(mutations)}')
    
    if args.print_consensus:
        print(consensus)
    
    if args.print_mutations:
        print(f'{",".join(mutations)}')
    
    return native_seq, consensus, match_case, mutations


if __name__ == '__main__':
    args = parse_args()
    native_seq, consensus, match_case, mutations = extract_consensus_top_sequences(args.fasta_file, args.top_p, args.sort_by, args.verbose)
    if native_seq is None:
        print('Error while extracting the consensus sequence. Something bad happened. Exit.')
        exit()
    
    with open(args.output, 'w') as f:
        f.write(f'# Parameters:\ttop_p={args.top_p}, sort_by={args.sort_by}\n')
        f.write(f'# Mutations:\t{", ".join(mutations)}\n')
        f.write(f'# Match case:\t{match_case}\n')
        f.write(f'>Native sequence\n{native_seq}\n')
        f.write(f'>Consensus sequence\n{consensus}\n')
