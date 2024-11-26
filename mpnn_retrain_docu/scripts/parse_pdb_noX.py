from Bio.PDB import PDBParser
import gzip
import numpy as np
import torch
import os,sys
import glob
import re
from scipy.spatial import KDTree
from itertools import combinations,permutations
import tempfile
import subprocess
import random, string

# --------------------------------
TM_ALIGN_PATH = "/home/iwe34/software/TMalign/TMalign"
# --------------------------------

RES_NAMES = [
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL'
]

RES_NAMES_1 = 'ARNDCQEGHILKMFPSTWYV'

to1letter = {aaa:a for a,aaa in zip(RES_NAMES_1,RES_NAMES)}
to3letter = {a:aaa for a,aaa in zip(RES_NAMES_1,RES_NAMES)}

ATOM_NAMES = [
    ("N", "CA", "C", "O", "CB"), # ala
    ("N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"), # arg
    ("N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"), # asn
    ("N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"), # asp
    ("N", "CA", "C", "O", "CB", "SG"), # cys
    ("N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"), # gln
    ("N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"), # glu
    ("N", "CA", "C", "O"), # gly
    ("N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"), # his
    ("N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"), # ile
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"), # leu
    ("N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"), # lys
    ("N", "CA", "C", "O", "CB", "CG", "SD", "CE"), # met
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"), # phe
    ("N", "CA", "C", "O", "CB", "CG", "CD"), # pro
    ("N", "CA", "C", "O", "CB", "OG"), # ser
    ("N", "CA", "C", "O", "CB", "OG1", "CG2"), # thr
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE2", "CE3", "NE1", "CZ2", "CZ3", "CH2"), # trp
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"), # tyr
    ("N", "CA", "C", "O", "CB", "CG1", "CG2") # val
]
        
idx2ra = {(RES_NAMES_1[i],j):(RES_NAMES[i],a) for i in range(20) for j,a in enumerate(ATOM_NAMES[i])}

aa2idx = {(r,a):i for r,atoms in zip(RES_NAMES,ATOM_NAMES) 
          for i,a in enumerate(atoms)}
aa2idx.update({(r,'OXT'):3 for r in RES_NAMES})


def writepdb(f, xyz, seq, bfac=None):
    # f = open(f, "w")
    f.seek(0)
    
    ctr = 1
    # seq = str(seq)
    L = len(seq)
    
    if bfac is None:
        bfac = np.zeros((L))

    idx = []
    for i in range(L):
        for j,xyz_ij in enumerate(xyz[i]):
            key = (seq[i],j)
            
            if key not in idx2ra.keys():
                continue
            if np.isnan(xyz_ij).sum()>0:
                continue
            r,a = idx2ra[key]
            f.write ("%-6s%5s  %-3s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(
                    "ATOM", ctr, a, r, 
                    "A", i+1, xyz_ij[0], xyz_ij[1], xyz_ij[2],
                    1.0, bfac[i,j] ) )
            if a == 'CA':
                idx.append(i)
            ctr += 1
    
    # f.write('TER\nEND\n')
    # f.close()
    # print(idx)
    f.flush()
    return np.array(idx)


def TMalign(chainA, chainB):
    # print('Run TMalign')
    
    # temp files to save the two input protein chains 
    # and TMalign transformation
    fA = tempfile.NamedTemporaryFile(mode='w+t', dir='/dev/shm')
    fB = tempfile.NamedTemporaryFile(mode='w+t', dir='/dev/shm')
    mtx = tempfile.NamedTemporaryFile(mode='w+t', dir='/dev/shm')
    
    # fA = 'tmp/' + ''.join([random.choice(string.ascii_uppercase) for _ in range(4)]) + '_cA.pdb'
    # fB = 'tmp/' + ''.join([random.choice(string.ascii_uppercase) for _ in range(4)]) + '_cB.pdb'
    # mtx = 'tmp/' + ''.join([random.choice(string.ascii_uppercase) for _ in range(4)]) + '_mtx.txt'

    # create temp PDB files keep track of residue indices which were saved
    idxA = writepdb(fA, chainA['xyz'], chainA['seq'], bfac=chainA['bfac'])
    idxB = writepdb(fB, chainB['xyz'], chainB['seq'], bfac=chainB['bfac'])
    
    # run TMalign
    # print('/home/iwe34/software/TMalign/TMalign %s %s -m %s'%(fA.name, fB.name, mtx.name))
    print('%s %s %s -m %s'%(TM_ALIGN_PATH, fA, fB, mtx))
    tm = subprocess.Popen('%s %s %s -m %s'%(TM_ALIGN_PATH, fA.name, fB.name, mtx.name), 
                          shell=True, 
                          stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE, 
                          encoding='utf-8')
    stdout, stderr = tm.communicate()
    # print(f'\n\n==========================\nSTDOUT\n{stdout}\n\n=========================')
    lines = stdout.split('\n')
    # print(f'\n\n==========================\nLINES: {len(lines)}\n{lines}\n\n=========================')
    # for i, item in enumerate(lines):
    #     print(i, item)
    
    # if TMalign failed
    if len(stderr) > 0:
        return None,None

    # parse transformation
    mtx.seek(0)
    # mtx = open(mtx, 'r')
    tu = np.fromstring(''.join(l[2:] for l in mtx.readlines()[2:5]), 
                       dtype=float, sep=' ').reshape((3,4))
    t = tu[:,0]
    u = tu[:,1:]

    # parse rmsd, sequence identity, and two TM-scores 
    # PROBLEM: The script for TMalign maybe has changed -> means
    # that the line indices in the original script are not true
    # rmsd = float(lines[16].split()[4][:-1])
    # seqid = float(lines[16].split()[-1])
    # tm1 = float(lines[17].split()[1])
    # tm2 = float(lines[18].split()[1])
    # print(f'Lines 16: {lines[16]} (empty!)')
    # print('***************')
    rmsd = float(lines[12].split()[4][:-1])
    seqid = float(lines[12].split()[-1])
    tm1 = float(lines[13].split()[1])
    tm2 = float(lines[13].split()[1])
    # print(rmsd, seqid, tm1, tm2)
    

    # parse alignment
    seq1 = lines[18]
    seq2 = lines[20]
    # print(seq1, '\n', seq2)
    
    ss1 = np.array(list(seq1.strip()))!='-'
    ss2 = np.array(list(seq2.strip()))!='-'
    #print(ss1)
    #print(ss2)
    mask = np.logical_and(ss1, ss2)

    alnAB = np.stack((idxA[(np.cumsum(ss1)-1)[mask]],
                      idxB[(np.cumsum(ss2)-1)[mask]]))

    alnBA = np.stack((alnAB[1],alnAB[0]))

    # clean up
    # fA.close()
    # fB.close()
    mtx.close()
    
    resAB = {'rmsd':rmsd, 'seqid':seqid, 'tm':tm1, 'aln':alnAB, 't':t, 'u':u}
    resBA = {'rmsd':rmsd, 'seqid':seqid, 'tm':tm2, 'aln':alnBA, 't':-u.T@t, 'u':u.T}
    
    return resAB,resBA


def get_tm_pairs(chains):
    """run TM-align for all pairs of chains"""

    tm_pairs = {}
    for A, B in combinations(chains.keys(),r=2):
        print(A, B)
        resAB,resBA = TMalign(chains[A],chains[B])
        #if resAB is None:
        #    continue
        tm_pairs.update({(A,B):resAB})
        tm_pairs.update({(B,A):resBA})
        
    # add self-alignments
    for A in chains.keys():
        L = chains[A]['xyz'].shape[0]
        aln = np.arange(L)[chains[A]['mask'][:,1]]
        aln = np.stack((aln,aln))
        tm_pairs.update({(A,A):{'rmsd':0.0, 'seqid':1.0, 'tm':1.0, 'aln':aln}})
    
    return tm_pairs
        

def parse_pdb(file:str, id:str):
    print(file)
    parser = PDBParser(QUIET=True)
    model = parser.get_structure(file=file, id='random')
    complete_structure = {}
    meta_data = {'seq': []}
    chains = []
    # per chain mpnn parsed two sequences
    # first is the real sequence, the second is the struture
    # sequence containing gaps for not resolved residues
    # did moritz structures contain gaps ???
    seq_chains = [] # per chain: [full, with_gaps]

    for chain in model.get_chains():
        chains.append(chain.get_id())
        
        full_seq, gaps_seq = [], []
        full_atoms_xyz = []
        bfac = []
        occ = []
        mask = []
        for residue in chain.get_residues():
            # adding the sequences for every chain
            # ToDo: add gaps for not resolved residues
            # skip water HOH
            aa_3_name = residue.get_resname()
            if aa_3_name == 'HOH':
                continue
            
            try:
                aa_1_name = to1letter[aa_3_name]
            except KeyError:
                continue
            full_seq.append(aa_1_name)
            gaps_seq.append(aa_1_name)
            
            # range(14) because 14 different atom types, maximum is in
            # Tryptophan (see list)
            atoms_xyz = [[np.nan for _ in range(3)] for _ in range(14)]
            atoms_bfac = [np.nan for _ in range(14)]
            atoms_occ = [0 for _ in range(14)]
            atoms_mask = [False for _ in range(14)]
            for atom in residue.get_atoms():
                # index for final list, same order for every residue
                # from the parse_cif_noX.py script
                try:
                    atom_idx = aa2idx[aa_3_name, atom.get_id()]
                except KeyError:
                    continue
                x, y, z = atom.get_coord()
                
                atoms_xyz[atom_idx] = [np.round(n, 3) for n in [x,y,z]]
                atoms_bfac[atom_idx] = atom.bfactor
                atoms_occ[atom_idx] = atom.occupancy
                atoms_mask[atom_idx] = True
            
            # atoms_xyz = np.array(atoms_xyz)
            full_atoms_xyz.append(atoms_xyz)
            bfac.append(atoms_bfac)
            occ.append(atoms_occ)
            mask.append(atoms_mask)
            # print(residue)
        
        chain_name = chain.get_id()
        tmp_d = {
                'seq': ''.join(gaps_seq),
                'xyz': np.array(full_atoms_xyz),
                'occ': np.array(occ),
                'bfac': np.array(bfac),
                'mask': np.array(mask)
        }
        complete_structure.update({chain_name: tmp_d})
        meta_data['seq'].append([''.join(full_seq), ''.join(gaps_seq)])
        
    # addint pseudo values for the other fields
    meta_data.update({'method' : 'biopython'})
    meta_data.update({'date': '2023-01-01'})
    meta_data.update({'resolution': 0.0})
    meta_data.update({'chains': chains})
    meta_data.update({'id': str(id)})
    meta_data.update({'asmb_chains': [','.join(chains)]})
    meta_data.update({'asmb_details': ['own_parsing']})
    meta_data.update({'asmb_method': ['none']})
    meta_data.update({'asmb_ids': ['1']})
    xform = [
        [[1. ,0. ,0. ,0.],
        [0. ,1. ,0. ,0.],
        [0. ,0. ,1. ,0.],
        [0. ,0. ,0. ,1.]]
    ]
    meta_data.update({'asmb_xform0': np.array(xform)})
    
    return complete_structure, meta_data        

IN = sys.argv[1]
OUT = sys.argv[2]
ID = str(sys.argv[3])

chains,metadata = parse_pdb(IN, ID)
ID = metadata['id']

tm_pairs = get_tm_pairs(chains)
if 'chains' in metadata.keys() and len(metadata['chains'])>0:
    chids = metadata['chains']
    tm = []
    for a in chids:
        tm_a = []
        for b in chids:
            tm_ab = tm_pairs[(a,b)]
            # print('---------------------')
            # print(tm_ab)
            if tm_ab is None:
                tm_a.append([0.0,0.0,999.9])
            else:
                tm_a.append([tm_ab[k] for k in ['tm','seqid','rmsd']])
        tm.append(tm_a)
    metadata.update({'tm':tm})

for k,v in chains.items():
    nres = (v['mask'][:,:3].sum(1)==3).sum()
    print(">%s_%s %s %s %s %d %d\n%s"%(ID,k,metadata['date'],metadata['method'],
                                       metadata['resolution'],len(v['seq']),nres,v['seq']))
    
    torch.save({kc:torch.Tensor(vc) if kc!='seq' else str(vc)
                for kc,vc in v.items()}, f"{OUT}_{k}.pt")

meta_pt = {}
for k,v in metadata.items():
    if "asmb_xform" in k or k=="tm":
        v = torch.Tensor(v)
    meta_pt.update({k:v})
torch.save(meta_pt, f"{OUT}.pt")
