import os
from globals import root, local
from utilities import *

from molecules import Reference

def load_references(force_reload=False):
    sprint("Loading references, force reload: ", force_reload)
    references = []
    if "molecules" in local.list() and not force_reload:
        for file in local.molecules:
            if file.endswith(".reference"):
                references.append(unpickle(file))
    if len(references) == 0:
        print1("No saved references found, generating them now")
        for file in os.listdir(root.references):
            if file.endswith((".pdb", ".cif")):
                references.append(Reference(os.path.join(root.references, file)))
    print1("{} references loaded:".format(len(references)))
    for reference in references:
        print2(reference)
    return references


