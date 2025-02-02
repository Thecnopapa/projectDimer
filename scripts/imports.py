import os
from globals import root, local
from utilities import *

from molecules import PDB, Monomer



def load_pickles(folder, extension = (".pickle")):
    print1("Looking for pickles in {}, with extension: {}".format(
        folder, extension))
    pickles = []
    if folder in local.list():
        for file in os.listdir(local[folder]):
            if file.endswith(extension):
                pickles.append(unpickle(file))
    return pickles


def load_from_pdb(pdb_folder = root.experimental, obj = PDB, pickle_folder = "molecules", pickle_extension = ".molecule", pdb_extension = (".pdb", ".pdb1", ".cif"), force_reload=False):
    sprint("Loading pdbs, force reload:", force_reload)
    loaded = []
    if not force_reload:
        loaded = load_pickles(pickle_folder, pickle_extension)
    if len(loaded) == 0:
        print1("No saved pickles found, importing from:", pdb_folder)
        for file in os.listdir(pdb_folder):
            if file.endswith(pdb_extension):
                loaded.append(obj(os.path.join(pdb_folder, file)))
    print1("{} objects loaded:".format(len(loaded)))
    for obj in loaded:
        print2(obj)
    return loaded

def load_monomers(molecules = None, folder = "monomers", extension = ".monomers",force_reload=False):
    sprint("Loading monomers, force reload:", force_reload)
    loaded = []
    if not force_reload:
        loaded = load_pickles(folder, extension)
    if len(loaded) == 0 and molecules is not None:
        print1("No saved monomers found, generating monomers now")
        for molecule in molecules:
            loaded.extend(molecule.get_monomers())
    print1("{} objects loaded:".format(len(loaded)))
    for obj in loaded:
        print2(obj)
    return loaded


