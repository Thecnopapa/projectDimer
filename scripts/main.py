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

tprint("Superposing Monomers")
columns = ["ID", "name", "chain"]
for reference in references:
    columns.extend(["rmsd_" + reference.name, "align_len_" + reference.name])
super_df = pd.DataFrame(columns = columns)
progress = ProgressBar(len(monomers))
for monomer in monomers:
    monomer.superpose(references, super_df)
    progress.add()
root["dataframes"] = "dataframes"
super_df.to_csv(os.path.join(root.dataframes,"super_df.csv"), header = True, index = False)









# Save and exit
tprint("Saving data")
progress = ProgressBar(len(references + molecules + monomers))
for file in references + molecules + monomers:
    file.pickle()
    file.export()
    progress.add()
eprint("Finished")