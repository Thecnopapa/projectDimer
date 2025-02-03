import os
from globals import root, local
from utilities import *

from molecules import PDB, Monomer



def load_pickles(folder, extension = (".pickle")):
    print1("Looking for pickles in {}, with extension: {}".format(
        folder, extension))
    pickles = []
    if folder in local.list():
        print("Files found:", len(os.listdir(local[folder])))
        for file in os.listdir(local[folder]):
            if file.endswith(extension):
                pickles.append(unpickle(os.path.join(local[folder],file)))
    return pickles


def load_from_files(pdb_folder = root.experimental, obj = PDB, pickle_folder = "molecules", pickle_extension = ".molecule", pdb_extension = (".pdb", ".pdb1", ".cif"), force_reload=False):
    sprint("Loading pdbs, force reload:", force_reload)
    loaded = []
    if not force_reload:
        loaded = load_pickles(pickle_folder, pickle_extension)
    if len(loaded) == 0:
        print1("No saved pickles found, importing from:", pdb_folder)
        progress = ProgressBar(len(os.listdir(pdb_folder)))
        for file in os.listdir(pdb_folder):
            if file.endswith(pdb_extension):
                loaded.append(obj(os.path.join(pdb_folder, file)))
            progress.add()
    print1("{} objects loaded:".format(len(loaded)))
    for obj in loaded:
        print2(obj)
    return loaded

def load_monomers(molecules = None, folder = "monomers", extension = ".monomer",force_reload=False):
    sprint("Loading monomers, force reload:", force_reload)
    loaded = []
    if not force_reload:
        loaded = load_pickles(folder, extension)
    if len(loaded) == 0 and molecules is not None:
        print1("No saved monomers found, generating monomers now")
        progress = ProgressBar(len(molecules))
        for molecule in molecules:
            loaded.extend(molecule.get_monomers())
            progress.add()
    print1("{} objects loaded:".format(len(loaded)))
    for obj in loaded:
        print2(obj)
    return loaded

def load_dimers(molecules = None, folder = "dimers", extension = ".dimer",force_reload=False):
    sprint("Loading dimers, force reload:", force_reload)
    loaded = []
    if not force_reload:
        loaded = load_pickles(folder, extension)
    if len(loaded) == 0 and molecules is not None:
        print1("No saved dimers found, generating dimers now")
        for molecule in molecules:
            loaded.extend(molecule.get_dimers())
    print1("{} objects loaded:".format(len(loaded)))
    for obj in loaded:
        print2(obj)
    return loaded

def download_pdbs(list_path, save_folder):
    import subprocess
    local[save_folder] = save_folder
    pdb_list = []
    with open(list_path, "r") as f:
        for line in f:
            pdb_list.extend(line.split(","))
    print(pdb_list)
    pdb_links = list_path+"_links.txt"
    open(pdb_links, "w")
    with open(pdb_links, "a") as f:
        for pdb in pdb_list:
            f.write("https://files.rcsb.org/download/{}.pdb\n".format(pdb))
    print(local[save_folder])
    subprocess.run(["wget","-i", pdb_links, "-P", local[save_folder]])



