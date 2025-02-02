import os
from globals import root, local
from utilities import *

from molecules import Reference, PDB

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
            if file.endswith((".pdb", ".cif", "pdb1")):
                references.append(Reference(os.path.join(root.references, file)))
    print1("{} references loaded:".format(len(references)))
    for reference in references:
        print2(reference)
    return references

def load_experimental(force_reload=False):
    sprint("Loading experimental files, force reload: ", force_reload)
    molecules = []
    if "molecules" in local.list() and not force_reload:
        for file in local.molecules:
            if file.endswith(".molecule"):
                molecules.append(unpickle(file))
    if len(molecules) == 0:
        print1("No saved molecules found, generating them now")
        for file in os.listdir(root.experimental):
            if file.endswith((".pdb", ".cif", "pdb1")):
                molecules.append(PDB(os.path.join(root.experimental, file)))
    print1("{} molecules loaded:".format(len(molecules)))
    for reference in molecules:
        print2(reference)
    return molecules


def load_monomers(force_reload=False):
    sprint("Loading monomer files, force reload: ", force_reload)
    monomers = []

