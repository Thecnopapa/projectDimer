import os

from utilities import *
import pandas as pd

import globals
globals.set_root("../")
if os.name == "nt":
    globals.set_local("C:/Users/iainv/localdata/_local/projectB")
elif os.name == "posix":
    globals.set_local("/localdata/iain/_local/projectB")

# Imports that need globals initialised:
from globals import root, local, vars
from dataframes import save_dfs, create_dfs



PROCESS_ALL = True
LARGE_DATASET = True # Delete all saved data previously to avoid errors




from imports import *
if "many_pdbs" not in local.list() and LARGE_DATASET:
    download_pdbs(os.path.join(root.pdb_lists,"list_1500"), "many_pdbs")


tprint("Loading files")
from imports import load_from_files
from molecules import Reference, PDB
references = load_from_files(root.references,
                             pickle_extension= ".reference",
                             get_monomer=True,
                             force_reload = PROCESS_ALL)
if LARGE_DATASET:
    molecules = load_from_files(local.many_pdbs, force_reload = PROCESS_ALL)
else:
    molecules = load_from_files(force_reload = PROCESS_ALL)
pickle(molecules)
eprint("Files loaded")


sprint("Creating dataframes")
create_dfs(references)
print1("Dataframes created")


tprint("Loading monomers")
from imports import load_monomers
monomers = load_monomers(molecules = molecules, force_reload=PROCESS_ALL)
pickle(monomers)
eprint("Monomers loaded")

tprint("Aligning monomers")
progress = ProgressBar(len(monomers))
for monomer in monomers:
    monomer.sequence_align(references, force_align = PROCESS_ALL)
    progress.add()
pickle(monomers)
save_dfs()
eprint("Monomers aligned")



tprint("Superposing monomers")
progress = ProgressBar(len(monomers))
for monomer in monomers:
    monomer.superpose(references, force_superpose = PROCESS_ALL)
    progress.add()
pickle(monomers)
save_dfs()
eprint("Monomers superposed")


tprint("Loading dimers")
from imports import load_dimers
dimers = load_dimers(molecules = molecules, force_reload=True)
eprint("Dimers loaded")






# Save and exit
save_dfs()
from visualisation import *
generate_charts()
tprint("Saving data")
all_files = references + molecules + monomers + dimers
progress = ProgressBar(len(all_files))
for file in all_files:
    file.pickle()
    file.export()
    progress.add()
eprint("Finished")