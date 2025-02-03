import os
import globals
from utilities import *
import pandas as pd

globals.set_root("../")
if os.name == "nt":
    globals.set_local("C:/Users/iainv/localdata/_local/projectB")
elif os.name == "posix":
    globals.set_local("/localdata/iain/_local/projectB")
from globals import root, local


PROCESS_ALL = False


from imports import *
#download_pdbs(os.path.join(root.pdb_lists,"list_1500"), "many_pdbs")

tprint("Loading files")
from imports import load_from_files
from molecules import Reference, PDB
references = load_from_files(root.references,
                           Reference,
                           pickle_extension= ".reference",
                           force_reload = PROCESS_ALL)

molecules = load_from_files(local.many_pdbs,
                          force_reload = PROCESS_ALL)
#molecules = load_from_pdb(force_reload = PROCESS_ALL)

for file in references+molecules:
    file.pickle()
eprint("Files loaded")


tprint("Loading monomers")
from imports import load_monomers
monomers = load_monomers(molecules = molecules, force_reload=PROCESS_ALL)
eprint("Monomers loaded")

tprint("Superposing monomers")
columns_raw = ["ID", "name", "chain"]
for reference in references:
    columns_raw.extend(["rmsd_" + reference.name, "align_len_" + reference.name])
super_raw_df = pd.DataFrame(columns = columns_raw)

colums_filtered = ["ID", "best_fit", "coverage (%)","rmsd", "identity (%)",  "Rx", "Ry", "Rz", "T"]
super_filtered_df = pd.DataFrame(columns = colums_filtered)

globals.vars["failed_df"] = pd.DataFrame(columns= ["ID", "reason", "details"])

progress = ProgressBar(len(monomers))
for monomer in monomers:
    monomer.superpose(references, super_raw_df, super_filtered_df,force_superpose = True)
    progress.add()
root["dataframes"] = "dataframes"
super_raw_df.to_csv(os.path.join(root.dataframes,"super_raw.csv"), header = True, index = False)
super_filtered_df.to_csv(os.path.join(root.dataframes,"super_filtered.csv"), header = True, index = False)
eprint("Monomers superposed")


tprint("Loading dimers")
from imports import load_dimers
dimers = load_dimers(molecules = molecules, force_reload=True)
eprint("Dimers loaded")






globals.vars.failed_df.to_csv(os.path.join(root.dataframes,"failed_df.csv"), header = True, index = False)
# Save and exit
tprint("Saving data")
all_files = references + molecules + monomers + dimers
progress = ProgressBar(len(all_files))
for file in all_files:
    file.pickle()
    file.export()
    progress.add()
eprint("Finished")