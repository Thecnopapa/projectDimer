import os
import numpy as np
import pandas as pd
import Bio.Align as Align
from Bio.Align import substitution_matrices


d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def create_aligner(matrix = "BLOSUM62"):

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load(matrix)
    return aligner

def get_sequence(chain, id=None, out_fasta=False):
    fasta = ""
    if id is not None and out_fasta:
        fasta = "> " + id + '\n'
    for res in chain.get_residues():
        try:
            fasta += d3to1[res.resname]
        except:
            pass
    if out_fasta:
        fasta += "\n"
    return fasta


def get_score(seq1, seq2, matrix = "BLOSUM62"):
    aligner = create_aligner(matrix=matrix)
    return aligner.score(seq1, seq2)

def get_alignment(seq1, seq2, matrix = "BLOSUM62", first_only=True):
    aligner = create_aligner(matrix=matrix)
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    if first_only:
        return aligner.align(seq1,seq2)[0]
    return aligner.align(seq1,seq2)

def get_alignment_map(keys_seq, target_seq, matrix = "BLOSUM62"):
    al = get_alignment(keys_seq, target_seq, matrix)
    #print(al)
    al_keys = al[0]
    al_target = al[1]
    al_map = {}
    n_keys = 0
    n_target = 0
    for k, target in zip(al_keys, al_target):
        #print(k, target)
        if k == "-":
            continue
        if target == "-":
            al_map[n_keys] = None
        else:
            al_map[n_keys] = n_target
            n_target += 1

        n_keys += 1
    #print(al_map)
    return al_map






