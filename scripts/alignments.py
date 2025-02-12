import os
import numpy as np
import pandas as pd
import Bio.Align as Align
from Bio.Align import substitution_matrices
from Bio.Align.PairwiseAligner import PairwiseAligner
from Bio.Seq import Seq

d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def create_aligner(matrix = 'BLASTP' ):

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
    fasta += "\n"

    return fasta



