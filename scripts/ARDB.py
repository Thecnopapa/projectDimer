import os, sys
import numpy as np
import pandas as pd
from utilities import *
from Globals import root, local, vars


# ARDB: https://androgendb.mcgill.ca/
# Citation to ARDB:
# Gottlieb B, Beitel LK, Nadarajah A, Palioura M, Trifiro M. 2012. The Androgen Receptor Gene Mutations Database(ARDB):2012 Update Human Mutation 33:887-894.
# https://onlinelibrary.wiley.com/doi/full/10.1002/humu.22046



class Mutation:
    def __init__(self, ids, mutation_type, phenotype, position, wt_res, mut_res):
        self.ids = ids
        self.id = self.ids[0]
        self.type = mutation_type
        self.phenotype = phenotype
        self.position = position
        self.wt_res = wt_res
        self.mut_res = mut_res

    def __repr__(self):
        return "Mutation {} ({}): \tr.{} {} --> {} \tleads to: {}".format(self.id, self.type, self.position, self.wt_res, self.mut_res, self.phenotype)


def parse_ardb_sequence():
    with open(os.path.join(root.data, 'ARDB_sequence'), 'r') as f:
        s = []
        for line in f:
            if line.startswith('#'):
                continue
            line = clean_string(line)
            if len(line) == 0:
                continue
            print()
            s.extend([l for l in line])
        return "".join(s)


ardb_sequence = parse_ardb_sequence()

def parse_ardb():
    ardb = pd.read_excel(os.path.join(root.data, 'ARDB.xls'))
    ardb.rename(columns={o: n for o,n in zip(ardb.columns, map(clean_string, ardb.columns))}, inplace=True)
    #[print(e) for e in set(ardb["Phenotype"])]
    #print(ardb.columns)
    ar_mutant_list = []
    for row in ardb.itertuples():
        #print(row.Mutationtype)
        try:

            if "LBD" in row.Exon:
                mut_type = clean_string(row.Mutationtype, allow=[]).lower()
                wt = row.Fromaminoacid.split("(")[1]
                wt_res = wt[:3]
                position = int(wt[3:])
                mut_res = row.Toaminoacid.split(")")[0]
                if len(mut_res) != 3 or mut_res in ["del" ,"dup"]:
                    continue
                phenotype = clean_string(row.Phenotype.replace("\n", " "), allow=[" "]).upper()
                redundant = False
                for mutation in ar_mutant_list:
                    if (mutation.type == mut_type
                            and mutation.phenotype == phenotype
                            and mutation.position == position
                            and mutation.wt_res == wt_res
                            and mutation.mut_res == mut_res):
                        redundant = True
                        mutation.ids.append(row.Accession)
                        break
                if redundant:
                    continue
                ar_mutant_list.append(Mutation(ids = [row.Accession],
                                               mutation_type=mut_type,
                                               phenotype=phenotype,
                                               position=position,
                                               wt_res=wt_res,
                                               mut_res=mut_res,
                                               ))
        except:
            continue

    for mutation in ar_mutant_list:
        print(mutation)
    return ar_mutant_list










