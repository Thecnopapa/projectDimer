import os
import globals
from utilities import *


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










# Save and exit
tprint("Saving data")
progress = ProgressBar(len(references + molecules + monomers))
for file in references + molecules + monomers:
    file.pickle()
    file.export()
    progress.add()
eprint("Finished")