
from utilities import *
import pandas as pd
import os
import setup
import time
import sys
# Imports that need globals initialised:
from Globals import root, local, vars

from imports import *

sprint("Loading References")
vars["references"] = load_references(force_reload=False)
print1("References loaded")


from faces import *
for reference in vars.references:
    if reference.name != "GR":
        continue
    from faces import get_ref_pca
    pca = get_ref_pca(reference)
    plot_atoms(reference.structure, pca)
