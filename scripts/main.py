import os
import globals
from utilities import *
import pandas as pd


globals.set_root("../")
if os.name == "nt":
    globals.set_local("C:/Users/iainv/localdata/_local/projectB")
elif os.name == "posix":
    globals.set_local("/localdata/iain/_local/projectB")

# Imports that need globals initialised:
from globals import root, local, vars
from dataframes import save_dfs

PROCESS_ALL = False


from imports import *
if "many_pdbs" not in local.list():
    download_pdbs(os.path.join(root.pdb_lists,"list_1500"), "many_pdbs")


tprint("Loading files")
from imports import load_from_files
from molecules import Reference, PDB
references = load_from_files(root.references,
                           Reference,
                           pickle_extension= ".reference",
                           force_reload = PROCESS_ALL)

#molecules = load_from_files(local.many_pdbs, force_reload = PROCESS_ALL)
molecules = load_from_files(force_reload = PROCESS_ALL)
pickle(molecules)
eprint("Files loaded")


sprint("Creating dataframes")
root["dataframes"] = "dataframes"
columns_raw = ["ID", "name", "chain"]
for reference in references:
    columns_raw.extend(["rmsd_" + reference.name, "align_len_" + reference.name])
vars["raw_monomers_df"] = pd.DataFrame(columns = columns_raw)

colums_filtered = ["ID", "best_fit", "coverage (%)","rmsd", "identity (%)",  "Rx", "Ry", "Rz", "T"]
vars["monomers_df"] = pd.DataFrame(columns = colums_filtered)

vars["failed_df"] = pd.DataFrame(columns= ["ID", "reason", "details"])
print1("Dataframes created")


tprint("Loading monomers")
from imports import load_monomers
monomers = load_monomers(molecules = molecules, force_reload=PROCESS_ALL)
pickle(monomers)
eprint("Monomers loaded")


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