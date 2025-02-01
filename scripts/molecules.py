import os
from utilities import *
from globals import root, local
from Bio.PDB import PDBParser, MMCIFParser


class BioObject:
    pickle_extension = '.pickle'
    pickle_folder = "other"
    name = "BioObject"


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
        if self.path.endswith('.cif'):
            self.structure = MMCIFParser(QUIET=True).get_structure(self.name[:4], self.path)
        else:
            self.structure = PDBParser(QUIET=True).get_structure(self.name[:4], self.path)



class PDB(BioObject):
    pickle_extension = '.molecule'
    pickle_folder = "molecules"
    def __init__(self, path):
        self.path = path
        self.name = clean_string(os.path.basename(path).split(".")[0], allow = ["_"])
        self.parse_structure()


class Reference(PDB):
    pickle_extension = '.reference'
    pickle_folder = "molecules"




class Monomer(BioObject):
    pickle_extension = '.monomer'
    pickle_folder = "monomers"
    pass


class Dimer(BioObject):
    pickle_extension = '.dimer'
    pickle_folder = "dimers"
    pass



