import os
import globals
from utilities import *
import pandas as pd

globals.set_root("../")
if os.name == "nt":
    globals.set_local("C:/Users/iainv/localdata/_local/projectB")
from globals import root, local



tprint("Loading files")
from imports import load_from_pdb
from molecules import Reference, PDB
references = load_from_pdb(root.references,
                           Reference,
                           ".reference",
                           force_reload = True)

molecules = load_from_pdb(force_reload = True)
eprint("Files loaded")


tprint("Loading monomers")
from imports import load_monomers
monomers = load_monomers(molecules = molecules, force_reload=True)
eprint("Monomers loaded")


columns_raw = ["ID", "name", "chain"]
for reference in references:
    columns_raw.extend(["rmsd_" + reference.name, "align_len_" + reference.name])
super_raw_df = pd.DataFrame(columns = columns_raw)

colums_filtered = ["ID", "best_fit", "coverage","rmsd", "align_len",  "Rx", "Ry", "Rz", "T"]
super_filtered_df = pd.DataFrame(columns = colums_filtered)
progress = ProgressBar(len(monomers))
for monomer in monomers:
    sprint(monomer.id)
    monomer.superpose(references, super_raw_df)
    monomer.choose_superposition(super_filtered_df)
    progress.add()
root["dataframes"] = "dataframes"
super_raw_df.to_csv(os.path.join(root.dataframes,"super_raw.csv"), header = True, index = False)
super_filtered_df.to_csv(os.path.join(root.dataframes,"super_filtered.csv"), header = True, index = False)

tprint("Superposing Monomers")








# Save and exit
tprint("Saving data")
progress = ProgressBar(len(references + molecules + monomers))
for file in references + molecules + monomers:
    file.pickle()
    file.export()
    progress.add()
eprint("Finished")