import os


from utilities import *
from globals import root, local, vars
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, StructureBuilder, Structure
import numpy as np
import pandas as pd


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

    def restore_dfs(self):
        for key, value in self.__dict__.items():
            if "_entries" in key:
                df_name = key.split("_entries")[0] + "_df"
                for id, contents in value.items():
                    vars[df_name].loc[id] = contents



    def __repr__(self):
        return "{} ({})".format(self.id, self.__class__.__name__)

    def parse_structure(self, parse_original = False):
        if self.path is None or parse_original:
            self.path = self.o_path
        if self.path.endswith('.cif'):
            self.structure = MMCIFParser(QUIET=True).get_structure(self.name[:4], self.path)
        else:
            self.structure = PDBParser(QUIET=True).get_structure(self.name[:4], self.path)
        for model in self.structure.get_list()[1:]:
            self.structure.__delitem__(model.id)
        for chain in self.structure[0].get_list():
            for residue in chain.get_list():
                for atom in residue.get_list():
                    if atom.name != "CA":
                        residue.__delitem__(atom.id)


    def export(self, subfolder = None, structure = None, extra_id = ""):
        exporting = PDBIO()
        if structure is None:
            structure = self.structure
        exporting.set_structure(structure)
        local["exports"] = "exports"
        os.makedirs(local.exports, exist_ok=True)
        if subfolder is not None:
            path = os.path.join(local.exports, subfolder)
        else:
            path = os.path.join(local.exports, "pdb_" + self.pickle_folder)
        os.makedirs(path, exist_ok=True)
        path = os.path.join(path, self.id+extra_id+".pdb")
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
        self.monomers = []
        self.dimers = []

    def get_monomers(self):
        if len(self.monomers) == 0 or True:
            self.monomers = []
            for chain in self.structure.get_chains():
                if len(chain.get_list()) > 50:
                    self.monomers.append(Monomer(self.name, chain.id, chain))
        return self.monomers

    def get_dimers(self):
        done = []
        if len(self.dimers) == 0 or True:
            self.dimers = []
            if len(self.monomers) >= 2:
                for monomer1 in self.monomers:
                    for monomer2 in self.monomers:
                        if monomer2.id != monomer1.id and monomer2.id not in done:
                            self.dimers.append(Dimer(monomer1, monomer2))
                    done.append(monomer1.id)
        return self.dimers




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
        self.super_path = None

        self.raw_monomers_entries = {}
        self.failed_entries = {}
        self.monomers_entries = {}


    def superpose(self, references, force_superpose=False):
        #sprint("Force superpose:", force_superpose)
        if self.super_path is None or force_superpose:
            #print(self.super_path)
            from superpose import superpose_single
            contents = [self.id, self.name, self.chain]
            self.superpositions = {}
            for reference in references:
                ref_name = reference.id
                if not ref_name in self.rmsds.keys() or force_superpose:
                    super_id = self.id + "_x_" + ref_name
                    #print(super_id, self.path)
                    data = superpose_single(super_id, self.path, reference.path)
                    #print(self.id)
                    #print(data)
                    if "rmsd" in data:
                        contents.extend([data["rmsd"], data["aligned_residues"]])
                        self.superpositions[ref_name] = data
                    else:
                        vars.failed_df.loc[self.id] = [self.id, "gesamt error", "are DISSIMILAR and cannot be reasonably aligned"]
                        contents.extend([99,0])
            vars.raw_monomers_df.loc[self.id] = contents
            self.raw_monomers_entries[self.id] = contents
            self.choose_superposition()

    def choose_superposition(self):
        finalists = []
        criteria = []
        coverages = []
        #print1("Choosing superpositions in", self.id)
        #print(self.superpositions.items())
        for ref_name, data in self.superpositions.items():
            #print(ref_name, data)
            data["coverage"] = data["aligned_residues"] / data["nres"]
            if  data["coverage"] >= 0.8:
                finalists.append((ref_name,data))
                criteria.append(data["identity"])
            else:
                coverages.append(data["coverage"])
                pass
        if len(finalists) >0:
            self.super_data = finalists[np.argmax(criteria)]
            data = self.super_data[1]
            contents = [self.id,self.super_data[0], round(data["coverage"]*100), data["rmsd"], round(data["identity"]*100)]
            contents.extend(data["t_matrix"].values())
            #print3(contents)
            vars.monomers_df.loc[self.id] = contents
            self.monomers_entries[self.id] = contents
            self.super_path = data["out_path"]
        else:
            self.super_path = ""
            vars.failed_df.loc[self.id] = [self.id, "no reference meets coverage (80%)", coverages]
            self.failed_entries[self.id] = [self.id, "no reference meets coverage (80%)", coverages]






class Dimer(BioObject):
    pickle_extension = '.dimer'
    pickle_folder = "dimers"

    def __init__(self, monomer1, monomer2):
        self.monomer1 = monomer1
        self.monomer2 = monomer2
        self.name = monomer1.name
        self.id = "{}_{}{}".format(self.name, monomer1.chain, monomer2.chain)
        self.incomplete = False

        self.failed_entries = {}



        if monomer1.super_path is None or monomer2.super_path is None or monomer1.super_path == "" or monomer2.super_path == "":
            vars.failed_df.loc[self.id] = [self.id, "Missing superposition", "At least one superposition is missing, {}:{}, {}:{}".format(monomer1.chain,monomer1.super_path,monomer2.chain,monomer2.super_path)]
            self.failed_entries[self.id] = [self.id, "Missing superposition", "At least one superposition is missing, {}:{}, {}:{}".format(monomer1.chain,monomer1.super_path,monomer2.chain,monomer2.super_path)]
            self.incomplete = True
        else:
            self.original_structure, self.replaced_structure, self.merged_structure = self.merge_structures(monomer1, monomer2)

    def merge_structures(self, monomer1, monomer2):
        print1("Merging monomers:", monomer1, monomer2)
        structureA = PDBParser(QUIET=True).get_structure(monomer1.id, monomer1.super_path)
        structureB = PDBParser(QUIET=True).get_structure(monomer2.id, monomer2.super_path)
        originalA = structureA[0]
        originalB = structureB[0]
        replacedA = structureA[1]
        replacedB = structureB[1]
        replacedA.id = replacedB.id = 1

        replacedA.get_list()[0].id = originalA.get_list()[0].id
        replacedB.get_list()[0].id = originalB.get_list()[0].id

        originalA.add(originalB.get_list()[0])
        replacedA.add(replacedB.get_list()[0])

        original_structure = Structure.Structure(self.id+"_o")
        original_structure.add(originalA)
        #print1(original_structure, original_structure.get_list())
        replaced_structure = Structure.Structure(self.id+"_r")
        replaced_structure.add(replacedA)
        #print1(replaced_structure, replaced_structure.get_list())
        merged_structure = Structure.Structure(self.id + "_m")
        merged_structure.add(originalA)
        merged_structure.add(replacedA)

        return original_structure, replaced_structure, merged_structure

    def export(self):
        if not self.incomplete:
            super().export(subfolder="dimers_original", structure=self.original_structure)
            super().export(subfolder="dimers_replaced", structure=self.replaced_structure, extra_id="_replaced")
            super().export(subfolder="dimers_merged", structure=self.merged_structure, extra_id="_merged")


