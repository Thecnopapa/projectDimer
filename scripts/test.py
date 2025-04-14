
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
    from faces import get_pca
    pca = get_pca(reference.structure)
    #plot_atoms(reference.structure, pca)

dimers = load_single_pdb("1A52", pickle_folder=local.dimers)
pcas = []
for dimer in dimers:
    pcas.append(dimer.pca)
plot_pcas(pcas)

