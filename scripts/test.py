
import pandas as pd
import os
import setup
import time
import sys
# Imports that need globals initialised:
from Globals import root, local, vars
from utilities import *
from imports import *

sprint("Loading References")
vars["references"] = load_references(force_reload=False)
print1("References loaded")


from faces import *
for reference in vars.references:
    if reference.name != "GR":
        continue
    from faces import get_pca
    pca = get_pca(reference.structure)
    plot_atoms(reference.structure, pca)

monomers = load_single_pdb("1OJ5", pickle_folder=local.monomers)
for monomer in monomers:
    print("Best Fit:", monomer.best_fit)
    pca = get_pca(monomer.replaced)
    plot_atoms(monomer.replaced, pca)

