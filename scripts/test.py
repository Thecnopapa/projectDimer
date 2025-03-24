
from utilities import *
import setup

# Imports that need globals initialised:
from Globals import root, local, vars


from imports import pickle, export, download_pdbs, load_references, load_single_pdb


# Load/Import molecule references
sprint("Loading References")
vars["references"] = load_references()
print1("References loaded")

for reference in vars["references"]:
    print(reference.__dict__)

from dataframes import save_dfs
