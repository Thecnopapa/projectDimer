import os


from utilities import *
from globals import root, local
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
import numpy as np


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
        file_name = "{}{}".format(self.id, self.pickle_extension)
        self.pickle_path = os.path.join(pickle_folder, file_name)
        with open(self.pickle_path, 'wb') as f:
            pickle.dump(self, f)

    def __repr__(self):
        return "{} ({})".format(self.id, self.__class__.__name__)

    def parse_structure(self):
        if self.path is None:
            path = self.o_path
        if path.endswith('.cif'):
            self.structure = MMCIFParser(QUIET=True).get_structure(self.name[:4], path)
        else:
            self.structure = PDBParser(QUIET=True).get_structure(self.name[:4], path)
        for model in self.structure.get_list()[1:]:
            self.structure.__delitem__(model.id)



    def export(self, subfolder = None):
        exporting = PDBIO()
        exporting.set_structure(self.structure)
        local["exports"] = "exports"
        if subfolder is not None:
            path = os.path.join(local.exports, subfolder)
        else:
            path = os.path.join(local.exports, "pdb_" + self.pickle_folder)
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
        self.export()

    def get_monomers(self):
        monomers = []
        for chain in self.structure.get_chains():
            if len(chain.get_list()) > 50:
                monomers.append(Monomer(self.name, chain.id, chain))
        return monomers





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
        self.export()

        self.rmsds = {}
        self.al_res = {}
        self.rotations = {}

    def superpose(self, references, df, force_superposed=False):
        from superpose import superpose_single
        contents = [self.id, self.name, self.chain]
        self.superpositions = {}
        for reference in references:
            ref_name = reference.id
            if not ref_name in self.rmsds.keys() or force_superposed:
                super_id = self.id + "_x_" + ref_name
                #print(super_id, self.path)
                data = superpose_single(super_id, self.path, reference.path)
                contents.extend([data["rmsd"], data["aligned_residues"]])
                self.superpositions[ref_name] = data
        df.loc[self.id] = contents

    def choose_superposition(self, df):
        finalists = []
        rmsds = []
        print1("Choosing superpositions in", self.id)
        #print(self.superpositions.items())
        for ref_name, data in self.superpositions.items():
            #print(ref_name, data)
            data["coverage"] = data["aligned_residues"] / data["nres"]
            if  data["coverage"] >= 0.8:
                finalists.append((ref_name,data))
                rmsds.append(data["rmsd"])
            else:
                print2("Dropped:", ref_name, round(data["coverage"],2))
        self.super_data = finalists[np.argmin(rmsds)]
        data = self.super_data[1]
        contents = [self.id,self.super_data[0], round(data["coverage"]*100), data["rmsd"], data["aligned_residues"]]
        contents.extend(data["t_matrix"].values())
        print3(contents)
        df.loc[self.id] = contents







class Dimer(BioObject):
    pickle_extension = '.dimer'
    pickle_folder = "dimers"
    pass



