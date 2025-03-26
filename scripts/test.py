
from utilities import *
import pandas as pd
import os
import setup

# Imports that need globals initialised:
from Globals import root, local, vars


from imports import pickle, export, download_pdbs, load_references, load_single_pdb


# Load/Import molecule references
sprint("Loading References")
vars["references"] = load_references()
print1("References loaded")

# Getting molecule list
sprint("Obtaining molecule list")

molecule_folder = local.many_pdbs

molecule_list = sorted(os.listdir(molecule_folder))
#print(len(vars.do_only), vars.do_only)
#[print(m) for m in sorted(molecule_list)]
print1("Molecule list obtained:", len(molecule_list), "molecules")
for m in molecule_list:
    if "lock" in m:
        sprint(".lock file detected:", m)
        continue
    filename = m.split(".")[0]
    sprint(filename)
    molecules = load_single_pdb(filename, local.molecules)
    for molecule in molecules:
        dimers = molecule.dimers
        for dimer in dimers:
            pass