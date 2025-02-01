import os
import globals



globals.set_root("../")
if os.name == "nt":
    globals.set_local("C:/Users/iainv/localdata/_local/projectB")
from globals import root, local




from imports import load_references, load_experimental
references = load_references(force_reload=True)

molecules = load_experimental(force_reload=True)


for reference in references:
    reference.parse_structure()











# Save and exit
for reference in references:
    reference.pickle()
for molecule in molecules:
    molecule.pickle()