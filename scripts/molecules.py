import os
from utilities import *
from Globals import root, local, vars
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
        for key, entries in self.__dict__.items():
            if "_entries" in key:
                df_name = key.split("_entries")[0] + "_df"
                for contents in entries:
                    if df_name in vars:
                        vars[df_name].loc[len(vars[df_name])] = contents



    def __repr__(self):
        return "{} ({})".format(self.id, self.__class__.__name__)

    def parse_structure(self, parse_original = False):
        if self.path is None or parse_original:
            self.path = self.o_path
        try:
            if self.path.endswith('.cif'):
                self.structure = MMCIFParser(QUIET=True).get_structure(self.name[:4], self.path)
            else:
                self.structure = PDBParser(QUIET=True).get_structure(self.name[:4], self.path)
        except:
            from scripts.imports import download_pdbs
            print1("Could not parse {}".format(self.path))
            import subprocess
            subprocess.run(["rm", self.path])
            print2("Re-downloading {}".format(self.name))
            download_pdbs(save_folder=os.path.dirname(self.path), pdb_list=[self.name])
            try:

                self.structure = PDBParser(QUIET=True).get_structure(self.name[:4], self.path)
            except:
                print2("Could not parse redownloaded {}".format(self.path))
            return False
        if len(self.structure.get_list()) >0: # TODO: Should add to failed df <<<
            for model in self.structure.get_list()[1:]:
                self.structure.__delitem__(model.id)
            for chain in self.structure.get_list()[0].get_list():
                for residue in chain.get_list():
                    for atom in residue.get_list():
                        if atom.name != "CA":
                            residue.__delitem__(atom.id)
        return True


    def export(self, subfolder = None, in_structure = None, extra_id = ""):
        exporting = PDBIO()
        if in_structure is None:
            structure = self.structure
        else: structure = in_structure
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
        if in_structure is None:
            self.path = path
        return path

class PDB(BioObject):
    pickle_extension = '.molecule'
    pickle_folder = "molecules"
    def __init__(self, path):
        self.o_path = path
        self.name = self.id = clean_string(os.path.basename(path).split(".")[0], allow = ["_"])
        self.complete = self.parse_structure()
        self.export()
        self.monomers = []
        self.dimers = []
        self.fractional = None
        self.fractional_path = None
        self.card = None
        self.params = None
        self.space_group = None # [group, key]

    
    def read_card(self):
        from symmetries import get_crystal, calculate_parameters, get_space_group
        try:
           self.card = get_crystal(self.o_path)
           self.params = calculate_parameters(self.card)
           self.space_group = get_space_group(self.card)
        except:
            self.card = None
            self.params = None
            self.space_group = None
        return self.card

    def get_monomers(self, as_reference=False):
        if len(self.monomers) == 0 or True:
            self.monomers = []
            for chain in self.structure.get_chains():
                if len(chain) > 200:
                    if as_reference:
                        return Reference(self.name, chain.id, chain)
                    self.monomers.append(Monomer(self.name, chain.id, chain))
        return self.monomers

    def get_dimers(self):
        done = []
        self.dimers = []
        if len(self.monomers) == 0:
            self.get_monomers()
        if len(self.monomers) >= 2:
            for monomer1 in self.monomers:
                for monomer2 in self.monomers:
                    if monomer2.id != monomer1.id and monomer2.id not in done:
                        self.dimers.append(Dimer(monomer1, monomer2))
                done.append(monomer1.id)
        return self.dimers


    def get_all_dimers(self, force = False):
        self.dimers = []
        self.read_card()
        from symmetries import find_relevant_mates
        self.mates = find_relevant_mates(self.structure, self.params, self.space_group[1])
        for mate in self.mates:
            mate.process(self)
            self.dimers.extend(mate.dimers)
        dimer_paths = []
        for dimer in self.dimers:
            dimer_paths.append(dimer.export())
        self.dimer_paths = dimer_paths
        return self.dimers








    def export_fractional(self):
        print2("Exporting fractional")
        if self.fractional is None:
            return None
        self.fractional_path = self.export(subfolder="fractional", in_structure=self.fractional, extra_id="_fractional")
        return self.fractional_path
    
    def export_neighbour(self):
        print2("Exporting neighbour")
        if self.neighbour is None:
            print3("Neighbour not found")
            return None
        #original_model = self.structure[0].copy()
        #original_model.id = 1
        #self.neighbour.add(original_model)
        self.neighbour_path = self.export(subfolder="neighbour", in_structure=self.neighbour, extra_id="_neighbour")
        return self.neighbour_path

    def generate_fractional(self):
        print1("Generating fractional")
        if self.params is None:
            return None
        if self.structure is None:
            return None
        self.fractional = self.structure.copy()
        from symmetries import convertFromOrthToFrac
        for atom in self.fractional.get_atoms():
            if atom.is_disordered():
                for d_atom in atom:
                    #print("dis", d_atom.coord, end=" -> ")
                    d_atom.coord=convertFromOrthToFrac(d_atom.coord, self.params)
                    #print(d_atom.coord)
            else:
                #print(atom.coord, end=" -> ")
                atom.coord=convertFromOrthToFrac(atom.coord, self.params)
                #print(atom.coord)
        self.export_fractional()
        return self.fractional

    def get_neighbour(self, force = False):
        sprint(self.id)
        print1("Getting neighbour")
        if self.params is None or self.fractional is None:
            print("Params or fractional not found:", self.params,self.fractional)
            return None
        from symmetries import  convertFromFracToOrth, find_nearest_neighbour_by_chain
        neighbour, lines = find_nearest_neighbour_by_chain(self.fractional,self.params, self.space_group[1], self.structure)
        self.lines = lines
        if neighbour is None:
            self.neighbour = None
            return None
        
        self.neighbour = neighbour
        self.export_neighbour()
        return self.neighbour


class Monomer(BioObject):
    pickle_extension = '.monomer'
    pickle_folder = "monomers"

    def __init__(self, name, chain, structure, is_mate = False):
        self.name = name
        self.is_mate = is_mate
        if self.is_mate:
            self.chain = chain.lower()
        else:
            self.chain = chain.upper()

        self.structure = structure
        self.id = "{}_{}".format(name, chain)

        self.export()
        self.rmsds = {}
        self.al_res = {}
        self.rotations = {}
        self.super_path = None
        self.super_data = None
        self.best_fit= None
        self.previews = None

        self.raw_monomers_entries = []
        self.failed_entries = []
        self.monomers_entries = []

        self.sequence = None

        self.scores = None



    def sequence_align(self, references, force_align = False):
        from alignments import create_aligner, get_sequence
        self.sequence = clean_string(get_sequence(self.structure))

        if self.scores is None or force_align:
            if len(self.sequence) > 0 :
                aligner = create_aligner(matrix="BLOSUM62")
                scores = []
                for reference in references:
                    scores.append(aligner.score(self.sequence, reference.sequence))
                vars.alignments_df.loc[self.id] = [self.id] + scores
                self.scores = scores
            else:
                vars.failed_df.loc[len(vars.failed_df)] = [self.id, "sequence", "sequence of length 0", self.sequence]
                self.failed_entries.append([self.id, "sequence", "sequence of length 0", self.sequence])





    def superpose(self, references = None, force_superpose=False):
        #sprint("Force superpose:", force_superpose)
        if self.super_path is None or force_superpose:
            #print(self.super_path)
            from superpose import superpose_single
            contents = [self.id, self.name, self.chain]
            self.superpositions = {}
            if references is None:
                references = self.top_refs_for_super
            for reference in references:
                ref_name = reference.name
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
                        vars.failed_df.loc[len(vars.failed_df)] = [self.id, "gesamt", "gesamt error", "are DISSIMILAR and cannot be reasonably aligned"]
                        contents.extend([99,0])
            vars.raw_monomers_df.loc[self.id] = contents
            self.raw_monomers_entries.append(contents)
            self.choose_superposition()

    def choose_superposition(self):
        finalists = []
        criteria = []
        coverages = []
        #print1("Choosing superpositions in", self.id)
        #print(self.superpositions.items())
        for ref_name, data in self.superpositions.items():
            #print(ref_name, data)
            data["coverage"] = (data["aligned_residues"] / data["nres"]*100, data["aligned_residues"] / data["ref_nres"]*100)
            if  data["coverage"] >= (80, 80):
                finalists.append((ref_name,data))
                criteria.append(data["identity"])
            else:
                coverages.append(data["coverage"])
                pass
        if len(finalists) > 0:
            self.super_data = finalists[np.argmax(criteria)]
            self.best_fit = self.super_data[0]
            data = self.super_data[1]
            contents = [self.id,self.super_data[0], data["coverage"][0], data["coverage"][1], data["rmsd"], round(data["identity"]*100)]
            contents.extend(data["t_matrix"].values())
            #print3(contents)
            vars.monomers_df.loc[self.id] = contents
            self.monomers_entries.append(contents)
            self.super_path = data["out_path"]
        else:
            self.super_path = ""
            vars.failed_df.loc[len(vars.failed_df)] = [self.id, "monomer", "No coverage", "No reference meets coverage (80%)" + str(coverages)]
            self.failed_entries.append([self.id, "monomer", "No coverage", "No reference meets coverage (80%)" + str(coverages)])



class Mate(BioObject):
    pickle_extension = ".mate"
    pickle_folder = "mates"

    def __init__(self, op_n, operation, params, fixed_chain, moving_chain):
        self.structure = None
        self.op_n = op_n
        self.operation = operation
        self.key = None
        self.params = params
        self.positions = {}
        self.f_chain = fixed_chain
        self.m_chain = moving_chain
        self.update_id()
        self.dimers = []

    def update_id(self):
        self.id = "{}_mate_{}_{}_{}".format(self.name, self.f_chain.id, self.op_n, self.m_chain.id)

    def process(self, parent):
        self.name = parent.name
        self.key = parent.space_group[1]
        self.update_id()
        print1("Processing mate:", self.id)
        
        self.reconstruct_mates()
        print2(self.dimers)
        
    def reconstruct_mates(self, min_contacts = 0):
        from symmetries import generate_displaced_copy, entity_to_orth, get_operation
        dimers = []
        for position in self.positions.values():
            if position["n_contacts"] >= min_contacts:
                fixed_monomer = Monomer(self.name, self.f_chain.id, entity_to_orth(self.f_chain, self.params))
                moved_mate = generate_displaced_copy(self.m_chain, distance = position["position"], rotation = get_operation(self.key, self.op_n))
                moved_mate = entity_to_orth(moved_mate, self.params)
                moved_mate = Monomer(self.name, self.m_chain.id, moved_mate, is_mate = True)
                dimers.append(Dimer(fixed_monomer, moved_mate))
        self.dimers = dimers
        return self.dimers


        
    



class Dimer(BioObject):
    pickle_extension = '.dimer'
    pickle_folder = "dimers"

    def __init__(self, monomer1, monomer2, position = None):
        self.monomer1 = monomer1
        self.monomer2 = monomer2
        self.name = monomer1.name
        self.position = position
        if position is None:
            p = ""
        else:
            p = "_" + str(position)
        self.id = "{}_{}{}{}".format(self.name, monomer1.chain, monomer2.chain, p)
        self.incomplete = False

        self.failed_entries = []

        if monomer1.best_fit == monomer2.best_fit:
            self.best_fit = monomer1.best_fit
        else:
            self.best_fit = "Missmatch"
        self.best_match = None

        if monomer1.super_path is None or monomer2.super_path is None or monomer1.super_path == "" or monomer2.super_path == "":
            if "failed_df" in vars:
                vars.failed_df.loc[len(vars.failed_df)] = [self.id, "dimer", "Missing superposition", "At least one superposition is missing, {}:{}, {}:{}".format(monomer1.chain,monomer1.super_path,monomer2.chain,monomer2.super_path)]
            self.failed_entries.append([self.id, "dimer", "Missing superposition", "At least one superposition is missing, {}:{}, {}:{}".format(monomer1.chain,monomer1.super_path,monomer2.chain,monomer2.super_path)])
            self.incomplete = True
        else:
            self.original_structure, self.replaced_structure, self.merged_structure = self.merge_structures()

    def merge_structures(self):

        monomer1 = self.monomer1
        monomer2 = self.monomer2

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

        monomer1.replaced = PDBParser(QUIET=True).get_structure(monomer1.id+"_r", monomer1.super_path)[1].get_list()[0]
        monomer2.replaced = PDBParser(QUIET=True).get_structure(monomer2.id+"_r", monomer1.super_path)[1].get_list()[0]

        monomer1.replaced.detach_parent()
        monomer2.replaced.detach_parent()


        return original_structure, replaced_structure, merged_structure

    def export(self):
        if not self.incomplete:
            self.original_path = super().export(subfolder="dimers_original", in_structure=self.original_structure)
            self.replaced_path = super().export(subfolder="dimers_replaced", in_structure=self.replaced_structure)
            self.merged_path = super().export(subfolder="dimers_merged", in_structure=self.merged_structure)
            return self.merged_path


    def summary(self):
        try:
            data = {
                "id": self.id,
                "name": self.name,
                "chain1": self.monomer1.chain,
                "chain2": self.monomer2.chain,
                "best_fit": self.best_fit,
                "best_match": self.best_match,
                "link": "https://files.rcsb.org/download/{}.pdb".format(self.name),
            }
        except:
            data = {"id": self.id,
                    "name": self.name}
        return data





class Reference(Monomer):
    pickle_extension = '.reference'
    pickle_folder = "refs"
