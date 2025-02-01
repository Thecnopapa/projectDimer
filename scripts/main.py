import os
import globals
from utilities import *


globals.set_root("../")
if os.name == "nt":
    globals.set_local("C:/Users/iainv/localdata/_local/projectB")
from globals import root, local



tprint("Loading files")
from imports import load_references, load_experimental
references = load_references(force_reload=True)

molecules = load_experimental(force_reload=True)
eprint("Files loaded")
















# Save and exit
for reference in references:
    reference.pickle()
for molecule in molecules:
    molecule.pickle()