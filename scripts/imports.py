import os

import numpy as np
from Globals import root, local, vars
from utilities import *

from molecules import PDB, Monomer,Reference



def load_pickles(folder, extension = (".pickle"), ignore_selection = False):
    print1("Looking for pickles in {}, with extension: {}".format(
        folder, extension))
    if "do_only" in vars:
        print2("Loading only:", vars.do_only)
    pickles = []
    if folder in local.list():
        print2("Files found:", len(os.listdir(local[folder])))
        selection = None
        if "do_only" in vars and not ignore_selection:
            if type(vars.do_only) is str:
                selection = vars.do_only.split(" ")
            elif type(vars.do_only) is list:
                selection = vars.do_only
        print(selection)
        if selection is not None:
            print2("Selection:", selection, "\n")

        progress = ProgressBar(len(os.listdir(local[folder])), silent=True, title=False)
        for file in sorted(os.listdir(local[folder])):
            if file.endswith(extension) and (selection is None or any([s in file.split(".")[0] for s in selection])):
                p = unpickle(os.path.join(local[folder],file))
                try:
                    p.restore_dfs()
                except:
                    pass
                pickles.append(p)
            progress.add(info=file)
    pickles.sort(key = lambda p: p.id)
    return pickles




def load_from_files(pdb_folder, load_class = PDB, ignore_selection = False, pickle_folder = "molecules",is_reference = False, pickle_extension = ".molecule", pdb_extension = (".pdb", ".pdb1", ".cif"), force_reload=False):
    sprint("Loading pdbs, force reload:", force_reload)
    try:
        print1("Loading only:", vars.do_only)
    except:
        pass
    loaded = []
    if not force_reload:
        loaded = load_pickles(pickle_folder, pickle_extension, ignore_selection=ignore_selection)
    if len(loaded) == 0:
        print1("No saved pickles found, importing from:", pdb_folder)
        progress = ProgressBar(len(os.listdir(pdb_folder)), silent=True, title=False)
        for file in sorted(os.listdir(pdb_folder)):
            selection = None
            if "do_only" in vars and not ignore_selection:
                if vars.do_only is str:
                    selection = vars.do_only.split(" ")
                elif vars.do_only is list:
                    selection = vars.do_only
            if file.endswith(pdb_extension) and (selection is None or any([s in file.split(".")[0] for s in selection])):
                obj = load_class(os.path.join(pdb_folder, file))
                if is_reference:
                    obj = obj.get_monomers(as_reference=True)
                loaded.append(obj)
            progress.add(info=file)
    print1("{} objects loaded:".format(len(loaded)))
    for obj in loaded:
        print2(obj)
    return loaded

class PickleIterator:
    def __init__(self, id_list, quiet=True, silent=True, **kwargs):
        self.id_list = sorted(id_list, reverse=False)
        self.kwargs = kwargs
        self.kwargs["quiet"] = quiet
        self.kwargs["silent"] = silent

    def __iter__(self):
        self.n = 0
        self.progress = ProgressBar(len(self.id_list), silent=self.kwargs["silent"])
        return self

    def __next__(self):
        if self.n >= len(self.id_list):
            raise StopIteration
        else:
            if not vars.quiet and not self.kwargs["quiet"]:
                sprint(self.id_list[self.n])
                self.progress.add(info=self.id_list[self.n])
            loaded = load_single_pdb(self.id_list[self.n], **self.kwargs)
            assert len(loaded) <= 1, "More than one ({}) object with the same id! \n{}:\n{}".format(len(loaded), self.id_list[self.n], loaded)

            self.n += 1
            if len(loaded) == 0:
                return None
            return loaded[0]

    def list(self):
        return list(iter(self))


def load_clusters(identifier = "all", onebyone=False, **kwargs):
    if onebyone:
        return load_list_1by1(identifier=identifier, pickle_folder=local.cluster_pickles, ignore_do_only=True, quiet = True, **kwargs)
    else:
        return load_single_pdb(identifier, pickle_folder=local.cluster_pickles, quiet=True, **kwargs)


def load_list_1by1(identifier = None, id_list=None, quiet=True, ignore_do_only=False, keep_extension=False, **kwargs):
    if id_list is None:
        assert "pickle_folder" in kwargs, "pickle folder or id list must be provided"
        id_list = os.listdir(kwargs["pickle_folder"])
        if not keep_extension:
            id_list = [file.split(".")[0] for file in id_list if "." in file]
    if "do_only" in vars and not ignore_do_only:
        if len(vars.do_only) > 0:
            id_list = [f for f in id_list if any([s in f for s in vars.do_only])]
    if identifier is not None:
        identifier = identifier.upper()
        id_list = [f for f in id_list if (identifier == "ALL" or identifier in f.upper()) and "lock" not in f and not any([bl in f for bl in vars.blacklist])]
    return iter(PickleIterator(id_list, quiet, **kwargs))

def load_single_pdb(identifier = "all", pickle_folder = None, pdb_folder = None, force_reload=False, object_class = PDB, quiet=False, first_only=False, **kwargs):
    if not quiet:
        print1("Loading single:", identifier, "-Force reload:", force_reload, "-Pickle folder:", pickle_folder)
    objects = []
    identifier = identifier.upper()
    if not force_reload and pickle_folder is not None:
        if not quiet:
            print2("Loading PDB pickle from:", pickle_folder)
        for file in os.listdir(pickle_folder):
            if (identifier == "ALL" or identifier in file.upper()) and "lock" not in file and not any([bl in file for bl in vars.blacklist]):
                p = unpickle(os.path.join(pickle_folder, file))
                try:
                    p.restore_dfs()
                except:
                    pass
                if object_class == Reference:
                    #p.restore_reference_dfs(reset=force_reload)
                    pass
                objects.append(p)
    if len(objects) == 0 and pdb_folder is not None:
        for file in os.listdir(pdb_folder):
            if (identifier == "ALL" or identifier in file.upper()) and "lock" not in file:
                print2("Generating PDB object from:", os.path.join(pdb_folder, file))
                objects.append(object_class(os.path.join(pdb_folder, file)))

    if not quiet:
        if len(objects) == 0:
            print2("No objects loaded")
        else:
            print2("Objects loaded: {}".format(objects))
    if len(objects) == 1 and first_only:
        objects = objects[0]
    return objects

def load_references(force_reload = False, identifier = "all"):
    local["refs"] = "pickles/refs"
    return load_single_pdb(identifier = identifier, pickle_folder = local.refs, pdb_folder=root.references, force_reload = force_reload, object_class = Reference)

def load_monomers(molecules = None, folder = "monomers", extension = ".monomer",force_reload=False):
    sprint("Loading monomers, force reload:", force_reload)
    loaded = []
    if not force_reload:
        loaded = load_pickles(folder, extension)
    if len(loaded) == 0 and molecules is not None:
        print1("No saved monomers found, generating monomers now")
        progress = ProgressBar(len(molecules), silent=True, title=False)
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

def download_pdbs(list_path:str=None , save_folder=None, pdb_list:list = None, terminal = True):
    sprint("Downloading PDBs:")
    local[save_folder] = save_folder
    if list_path is not None:
        pdb_list = []
        with open(list_path, "r") as f:
            for line in f:
                pdb_list.extend(line.split(","))
        pdb_links = list_path + "_links.txt"
    elif pdb_list is not None:
        pdb_links = os.path.join(root.pdb_lists, "temp_list_links.txt")
    else:
        print("Please provide either a list or a path to a PDB list")
        return
    #print(pdb_list)
    print1("{} pdbs for download".format(len(pdb_list)))
    open(pdb_links, "w")
    with open(pdb_links, "a") as f:
        for pdb in pdb_list:
            if not any([bl in pdb for bl in vars.blacklist]):
                f.write("https://files.rcsb.org/download/{}.pdb\n".format(pdb))
    print1("Links saved at {}".format(pdb_links))

    if not terminal:
        import wget
        counter =0
        print(pdb_links)
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
    progress = ProgressBar(len(list), silent=True, title=False)
    sprint("Pickling...")
    for item in list:
        item.pickle()
        progress.add()

def export(list):
    sprint("Exporting...")
    progress = ProgressBar(len(list), silent=True, title=False)
    for item in list:
        item.export()
        progress.add()

def import_references():
    sprint("loading references")
    return load_from_files(root.references,
                                 pickle_extension=".reference",
                                 is_reference=True,
                                 ignore_selection=True,
                                 force_reload=False)

if __name__ == "__main__":
    import setup
    from Globals import *

    imported = {}
    for folder in os.listdir(local.pickles):
        sprint("Extracting {}...".format(folder))
        progress = ProgressBar(len(os.listdir(os.path.join(local.pickles, folder))))
        if os.path.isdir(os.path.join(local.pickles, folder)):
            for file in os.listdir(os.path.join(local.pickles, folder)):
                if os.path.isfile(os.path.join(local.pickles, folder, file)):
                    try:
                        imported[folder].append(unpickle(os.path.join(local.pickles, folder, file)))
                    except:
                        imported[folder] = [unpickle(os.path.join(local.pickles, folder, file))]
                    progress.add(info=file)


