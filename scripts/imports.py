import os


from globals import root, local, vars
from utilities import *

from molecules import PDB, Monomer,Reference



def load_pickles(folder, extension = (".pickle"), ignore_selection = False):
    print1("Looking for pickles in {}, with extension: {}".format(
        folder, extension))
    if "do_only" in vars:
        print2("Loading only:", vars.do_only)
    pickles = []
    if folder in local.list():
        print("Files found:", len(os.listdir(local[folder])))
        progress = ProgressBar(len(os.listdir(local[folder])))
        for file in os.listdir(local[folder]):
            selection = None
            if "do_only" in vars and not ignore_selection:
                if len(vars.do_only) > 0:
                    selection = vars.do_only.split(" ")
            if file.endswith(extension) and (selection is None or file.split(".")[0] in selection):
                p = unpickle(os.path.join(local[folder],file))
                p.restore_dfs()
                pickles.append(p)
                progress.add()
    pickles.sort(key = lambda p: p.id)
    return pickles


def load_from_files(pdb_folder = root.experimental, obj = PDB, ignore_selection = False, pickle_folder = "molecules",is_reference = False, pickle_extension = ".molecule", pdb_extension = (".pdb", ".pdb1", ".cif"), force_reload=False):
    sprint("Loading pdbs, force reload:", force_reload)
    print1("Loading only:", vars.do_only)
    loaded = []
    if not force_reload:
        loaded = load_pickles(pickle_folder, pickle_extension, ignore_selection=ignore_selection)
    if len(loaded) == 0:
        print1("No saved pickles found, importing from:", pdb_folder)
        progress = ProgressBar(len(os.listdir(pdb_folder)))
        for file in os.listdir(pdb_folder):
            selection = None
            if "do_only" in vars and not ignore_selection:
                if len(vars.do_only) > 0:
                    selection = vars.do_only.split(" ")
            if file.endswith(pdb_extension) and (selection is None or file.split(".")[0] in selection):
                obj = obj(os.path.join(pdb_folder, file))
                if is_reference:
                    obj = obj.get_monomers(as_reference=True)
                loaded.append(obj)
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

def download_pdbs(list_path, save_folder, terminal = False):

    sprint("Downloading PDBs:")

    local[save_folder] = save_folder
    pdb_list = []
    with open(list_path, "r") as f:
        for line in f:
            pdb_list.extend(line.split(","))
    #print(pdb_list)
    print1("{} pdbs for download".format(len(pdb_list)))
    pdb_links = list_path+"_links.txt"
    open(pdb_links, "w")
    with open(pdb_links, "a") as f:
        for pdb in pdb_list:
            f.write("https://files.rcsb.org/download/{}.pdb\n".format(pdb))
    print1("Links saved at {}".format(pdb_links))
    if not terminal:
        import wget
        counter =0
        with open(pdb_links, "r") as f:
            for line in f:
                print2("({}/{})".format(counter,len(pdb_list)),line)
                try:
                    wget.download(line, out=local[save_folder])
                except:
                    print3("Failed to import", line)
                counter+=1
    else:
        import subprocess
        subprocess.run(["wget","-i", pdb_links, "-P", local[save_folder]])


def pickle(list):
    for item in list:
        item.pickle()

def export(list):
    for item in list:
        item.export()



