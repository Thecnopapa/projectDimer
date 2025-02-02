import os
from utilities import *
from globals import root, local
from Bio.PDB import PDBParser, MMCIFParser, PDBIO


class BioObject:
    pickle_extension = '.pickle'
    pickle_folder = "other"
    name = "BioObject"
    path = None


    def pickle(self):
        import pickle
        local["pickles"] = "pickles"
        pickle_folder = os.path.join(local.pickles, self.pickle_folder)
        os.makedirs(pickle_folder, exist_ok=True)
        file_name = "{}{}".format(self.name, self.pickle_extension)
        self.pickle_path = os.path.join(pickle_folder, file_name)
        with open(self.pickle_path, 'wb') as f:
            pickle.dump(self, f)

    def __repr__(self):
        return "{} ({})".format(self.name, self.__class__.__name__)

    def parse_structure(self):
        if self.path is None:
            path = self.o_path
        if path.endswith('.cif'):
            self.structure = MMCIFParser(QUIET=True).get_structure(self.name[:4], path)
        else:
            self.structure = PDBParser(QUIET=True).get_structure(self.name[:4], path)


    def export(self, subfolder = None):
        exporting = PDBIO()
        exporting.set_structure(self.structure)
        local["exports"] = "exports"
        if subfolder is not None:
            path = os.path.join(local.exports, subfolder)
        else:
            path = os.path.join(local.exports, self.pickle_folder)
        os.makedirs(path, exist_ok=True)
        path = os.path.join(path, self.id+".pdb")
        exporting.save(path)
        self.path = path

class PDB(BioObject):
    pickle_extension = '.molecule'
    pickle_folder = "molecules"
    def __init__(self, path):
        self.o_path = path
        self.name = self.id = clean_string(os.path.basename(path).split(".")[0], allow = ["_"])
        self.parse_structure()

    def get_monomers(self):
        pass


class Reference(PDB):
    pickle_extension = '.reference'
    pickle_folder = "molecules"




class Monomer(BioObject):
    pickle_extension = '.monomer'
    pickle_folder = "monomers"

    def __init__(self, name, chain, structure):
        self.name = name
        self.chain = chain
        self.structure = structure
        self.id = "{}_{}".format(name, chain)
        self.path = self.export()


class Dimer(BioObject):
    pickle_extension = '.dimer'
    pickle_folder = "dimers"
    pass



