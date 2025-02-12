# Project_B by Iain Visa @ IBMB-CSIC / UB

# Essential imports
import os
import pandas as pd
from utilities import *


# Initialise globals
import globals
if os.name == "nt":
    globals.set_root("../")
    globals.set_local("C:/Users/iainv/localdata/_local/projectB")
elif os.name == "posix":
    import platform
    if "aarch" in platform.platform():
        globals.set_local("localdata/projectB")
        globals.set_root("projectB")
    else:
        globals.set_root("../")
        globals.set_local("/localdata/iain/_local/projectB")


# Imports that need globals initialised:
from globals import root, local, vars
print(local.list())
print(root.list())
from dataframes import save_dfs, create_dfs
from imports import pickle, export


# Some essential variables
PROCESS_ALL = False # Ignore saved pickles and generate everything from scratch
LARGE_DATASET = True # Delete all saved data previously to avoid errors
DO_ONLY = "" # Names of pdbs to be processed (CAPS sensitive, separated by space) e.g "5N10 1M2Z"


# Set up
vars["do_only"] = DO_ONLY
print(local.list())


# Download large dataset
tprint("Downloading large data")
from imports import download_pdbs
if "many_pdbs" not in local.list() and LARGE_DATASET:
    download_pdbs(os.path.join(root.pdb_lists,"list_1500"), "many_pdbs", terminal = True)
eprint("Large data downloaded")


# Load/Import molecule files
tprint("Loading files")
from imports import load_from_files
from molecules import Reference, PDB
references = load_from_files(root.references,
                             pickle_extension= ".reference",
                             is_reference=True,
                             ignore_selection = True,
                             force_reload = PROCESS_ALL)
if LARGE_DATASET:
    molecules = load_from_files(local.many_pdbs, force_reload = PROCESS_ALL)
else:
    molecules = load_from_files(force_reload = PROCESS_ALL)
pickle(molecules)
eprint("Files loaded")


# Create dataframes
sprint("Creating dataframes")
create_dfs(references)
for df in os.listdir(root.dataframes):
    print1(df)
print1("Dataframes created")


# Load/Generate monomer files
tprint("Loading monomers")
from imports import load_monomers
monomers = load_monomers(molecules = molecules, force_reload=PROCESS_ALL)
pickle(monomers)
eprint("Monomers loaded")


# Align references to monomers
tprint("Aligning monomers")
progress = ProgressBar(len(monomers))
for monomer in monomers:
    pass
    #monomer.sequence_align(references, force_align = PROCESS_ALL)
    progress.add()
pickle(monomers)
save_dfs()
eprint("Monomers aligned")


# Superpose references to monomers
tprint("Superposing monomers")
progress = ProgressBar(len(monomers))
for monomer in monomers:
    monomer.superpose(references, force_superpose = PROCESS_ALL)
    progress.add(info=monomer.id)
pickle(monomers)
save_dfs()
eprint("Monomers superposed")


# Load/Generate dimer files
tprint("Loading dimers")
from imports import load_dimers
dimers = load_dimers(molecules = molecules, force_reload=PROCESS_ALL)
pickle(dimers)
eprint("Dimers loaded")


# Processing dimers
tprint("Processing dimers")
progress = ProgressBar(len(dimers))
for dimer in dimers:
    #dimer
    progress.add()
pickle(dimers)
save_dfs()
eprint("Dimers processed")







# Save and exit
tprint("Saving data")
save_dfs() # Save dataframes and generate figures

all_files = references + molecules + monomers + dimers # select all loaded files
progress = ProgressBar(len(all_files))
for file in all_files:
    file.pickle() # Save the object
    file.export() # Save the pdb
    progress.add()
eprint("Finished")

# Done xd
