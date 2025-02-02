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




for file in references + molecules + monomers:
    file.export()








# Save and exit
tprint("Saving pickles")
for reference in references:
    reference.pickle()
for molecule in molecules:
    molecule.pickle()
for monomer in monomers:
    monomer.pickle()
eprint("Finished")