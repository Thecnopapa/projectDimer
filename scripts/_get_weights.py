import os
import sys



from Globals import root, local, vars
from utilities import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


import setup
from Globals import root, local, vars
from imports import *
from pyMol import *


file = "./GR-weights.tsv"

ref = load_list_1by1(identifier="GR", pickle_folder=local.refs).list()[0]
cluster  = load_clusters(identifier="GR-all-all", first_only=False)[0]
print("REF:")
print(ref)
print(ref.__class__)
print(dir(ref))
print(ref.structure)
print()
print("CLUSTER:")
print(cluster)
print(dir(cluster))
print(cluster.__class__)
print(cluster.preference_array)
print()



atoms = list(ref.structure.get_atoms())
weights = cluster.preference_array

assert sum([1 for _ in atoms]) == len(weights)


print("Atoms:")

with open(file, "w") as f:
    print(file)

    for n, (atom, weight) in enumerate(zip(atoms, weights)):
        #print(n, atom, weight)
        if n == 0:
            print(dir(atom))
            print(atom.id, atom.get_full_id())
            #exit()
        resnum = atom.get_full_id()[3][1]
        chain = atom.get_full_id()[2]
        line = f"{n}\t{resnum}\t{chain}\t{weight:5.3f}\n"
        print(line, end="")
        f.write(line)