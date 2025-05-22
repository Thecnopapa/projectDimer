import os

from pyMol import pymol_colour
from surface import get_dimer_sasa, get_monomer_sasa
from symmetries import entity_to_orth, print_all_coords, convertFromFracToOrth
from utilities import *
from Globals import root, local, vars
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, StructureBuilder, Structure,  Model
import numpy as np
import pandas as pd
import string
from maths import *



class BioObject:
    pickle_extension = '.pickle'
    pickle_folder = "other"
    name = "BioObject"
    path = None


    def pickle(self):
        import pickle
        local["pickles"] = "pickles"
        #print(local.pickles, self.pickle_folder)
        pickle_folder = os.path.join(local.pickles, self.pickle_folder)
        #print(pickle_folder)

        os.makedirs(pickle_folder, exist_ok=True)
        local[self.pickle_folder] = "pickles/{}".format(self.pickle_folder)
        file_name = "{}{}".format(self.id, self.pickle_extension)
        self.pickle_path = os.path.join(pickle_folder, file_name)
        #print(self.pickle_path)
        with open(self.pickle_path, 'wb') as f:
            pickle.dump(self, f)


    def restore_dfs(self):
        for key, entries in self.__dict__.items():
            if "_entries" in key:
                df_name = key.split("_entries")[0] + "_df"
                for contents in entries:
                    if df_name in vars:
                        vars[df_name].loc[len(vars[df_name])] = contents

    def restore_reference_dfs(self, reset = False):
        assert "Reference" in str(self.__class__)
        if reset:
            return
        if "clustering" not in vars:
            vars["clustering"] = {}
        for key, df in self.__dict__.items():
            if "df" in key:
                try:
                    df_folder = key.split("_")[0]
                    if df_folder not in vars["clustering"]:
                        vars["clustering"][df_folder] = {}
                    vars["clustering"][df_folder][self.name] = df
                except:
                    print(key)
                    print(df_folder)
                    quit()



    def __repr__(self):
        if "best_fit" in self.__dict__:
            return "{} (Best Fit: {})".format(self.id, self.best_fit)
        else:
            return "{} ({} at {})".format(self.id, self.__class__.__name__, id(self))

    def parse_structure(self, parse_original = False, calculate_sasa = False, n_points=100, radius=3):
        if self.path is None or parse_original:
            self.path = self.o_path
        try:
            if self.path.endswith('.cif'):
                self.structure = MMCIFParser(QUIET=True).get_structure(self.name[:4], self.path)
            else:
                self.structure = PDBParser(QUIET=True).get_structure(self.name[:4], self.path)
        except:
            from imports import download_pdbs
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
                if chain.id not in string.ascii_uppercase:
                    try:
                        chain.id = string.ascii_uppercase.index(chain.id)
                    except:
                        chain.id = [letter for letter in string.ascii_uppercase if letter not in [c.id for c in self.structure.get_chains()]][0]
                if calculate_sasa:
                    from Bio.PDB.SASA import ShrakeRupley
                    sr = ShrakeRupley(n_points=n_points, probe_radius=radius)
                    sr.compute(chain, level="R")
                    self.sasa_info = dict(n_points=n_points, probe_radius=radius)
                for residue in chain.get_list():
                    if residue.id[0] != " ":
                        chain.__delitem__(residue.id)
                        continue
                    for atom in residue.get_list():
                        if atom.name != "CA" or atom.bfactor > 90:
                            residue.__delitem__(atom.id)
                    if len(residue.get_list()) == 0:
                        chain.__delitem__(residue.id)
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

    def copy(self):
        print("Copying {}".format(self))
        copied =self.__class__.__new__(self.__class__)
        import copy
        for key, value in self.__dict__.items():
            try: copied.__setattr__(key, value.copy())
            except:
                copied.__setattr__(key, copy.copy(value))
            #try: copy.__setattr__(key, value.copy())
            #except: setattr(copy, key, value); print("Copying {} failed".format(key))
        return copied

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
        self.setup()

    def setup(self):
        self.read_card()


    
    def read_card(self):
        from symmetries import get_crystal, calculate_parameters, get_space_group
        try:
           self.card = get_crystal(self.o_path)
           self.params = calculate_parameters(self.card)
           self.space_group = get_space_group(self.card)
           self.key = self.space_group[1]
           print2("Crystal card read:", self.space_group, self.key)
        except:
            self.card = None
            self.params = None
            self.space_group = None
            self.key = None

        return self.card

    def get_monomers(self, as_reference=False):
        print("Using a deprecated feature, please don't")
        quit()
        if len(self.monomers) == 0 or True:
            self.monomers = []
            for chain in self.structure.get_chains():
                if len(chain) > 200:
                    if as_reference:
                        return Reference(self.name, chain.id, chain)
                    self.monomers.append(Monomer(self.name, chain.id, chain, self))
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


    def get_all_dimers(self, force = False, minimum_chain_length = 100, contact_distance = 8, min_contacts = 0):
        if len(self.dimers) > 0 and not force:
            return self.dimers
        self.dimers = []
        self.read_card()
        from symmetries import find_relevant_mates
        self.mates = find_relevant_mates(self, self.structure, self.params, self.key,
                                         minimum_chain_length=minimum_chain_length,
                                         contact_distance=contact_distance,
                                         min_contacts=min_contacts)
        if self.mates is None: 
            return []

        self.mate_paths = []
        self.contacts = []
        for mate in self.mates:
            if mate is None:
                continue
            mate.process(self, min_contacts=min_contacts)
            self.dimers.extend(mate.dimers)
            self.contacts.extend(mate.contacts)
            #self.mate_paths.extend(mate.paths)

        dimer_paths = []
        for dimer in self.dimers:
            dimer_paths.append(dimer.export())
        self.dimer_paths = dimer_paths
        return self.dimers




class Monomer(BioObject):
    pickle_extension = '.monomer'
    pickle_folder = "monomers"

    def __init__(self, name, chain, frac_structure, parent, is_mate = False, op_n = None, position = None, parent_monomer = None, sasa=False):
        self.name = name
        self.is_mate = is_mate
        self.fractional_structure = frac_structure
        self.params = parent.params
        self.key = parent.key
        self.op_n = op_n
        self.position = position
        self.structure = entity_to_orth(self.fractional_structure.copy(), self.params)
        self.parent_monomer = parent_monomer
        self.face = None
        self.faces = None
        #self.parent = parent

        if self.is_mate:
            self.chain = chain.lower()
            #print(self.position["position"], self.op_n, self.name)
            self.extra_id = "_{}_{}".format(self.op_n, clean_string(self.position["position"], allow=["-"]))
            self.contacts = parent.contacts

        else:
            self.chain = chain.upper()
            self.extra_id = ""

        self.structure.id = self.chain
        self.id = "{}_{}{}".format(name, chain, self.extra_id)
        self.export()



        self.previews = None
        self.raw_monomers_entries = []
        self.failed_entries = []
        self.monomers_entries = []
        self.sequence = None
        self.scores = None
        self.sasas = None
        self.rmsds = {}
        self.super_path = None
        self.super_data = None
        self.best_fit = None
        self.replaced = None


        if self.is_mate:
            self.rmsds = self.parent_monomer.rmsds
            self.super_data = self.parent_monomer.super_data
            self.best_fit = self.parent_monomer.best_fit
            self.super_path = self.move_parent_superposition(self.parent_monomer.super_path)
            if self.super_path is not None:
                self.replaced = PDBParser(QUIET=True).get_structure(self.id, self.super_path).get_list()[1].get_list()[0]
                if sasa:
                    self.sasas = self.parent_monomer.sasas.copy()

        else:
            self.superpose()
            if sasa:
                from surface import get_monomer_sasa
                get_monomer_sasa(self)

        from faces import get_face_coms
        if self.replaced is not None:
            if self.best_fit == "GR":
                self.face_coms = get_face_coms(self)
            self.get_monomer_pca()

        self.pickle()


    def get_monomer_pca(self):
        from maths import find_com
        from faces import get_pca
        self.com = find_com(self.replaced.get_atoms())
        self.pca = dict(pca=get_pca(self.replaced, com=self.com, closer_to="N"),
                        com=self.com,
                        chain=self.chain)
        return self.pca

    def move_parent_superposition(self,super_path):
        from symmetries import entity_to_orth, entity_to_frac, coord_operation_entity, print_all_coords
        if super_path is None:
            return super_path

        original_chain = self.structure.copy()
        #print_all_coords(original_chain)
        #print("original_chain", original_chain, len(list(original_chain.get_residues())))
        parent_structure = PDBParser(QUIET=True).get_structure(self.id, super_path)
        #print("parent_structure", parent_structure, parent_structure.get_list())
        #print_all_coords(parent_structure)
        new_structure = Structure.Structure(self.id)
        model0 = Model.Model(0)
        model1 = Model.Model(1)
        model0.add(original_chain)
        assert len(parent_structure.get_list()[1].get_list()) == 1
        moving_superposition = parent_structure.get_list()[1].get_list()[0]
        moving_superposition.id = self.chain

        #print("moving_superposition", moving_superposition, len(list(moving_superposition.get_residues())))
        #print_all_coords(moving_superposition)
        moving_superposition = entity_to_frac(moving_superposition, self.params)
        #print_all_coords(moving_superposition)
        #print(self.key, self.op_n, self.position["position"])
        moving_superposition = coord_operation_entity(moving_superposition, self.key, self.op_n, self.position["position"])
        #print_all_coords(moving_superposition)
        moving_superposition = entity_to_orth(moving_superposition, self.params)
        #print_all_coords(moving_superposition)
        #print("moved_superposition", moving_superposition)

        model1.add(moving_superposition)
        #print("model0", model0, model0.get_list())
        #print("model1", model1, model1.get_list())
        new_structure.add(model0)
        new_structure.add(model1)
        #print("new_structure", new_structure, new_structure.get_list())

        exporting = PDBIO()
        exporting.set_structure(new_structure)
        path = os.path.join(local.super_raw, os.path.basename(super_path).split(".")[0] + self.extra_id + "_merged.pdb")
        exporting.save(path)
        return path






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
        print1("Superposing:", self.id)
        if self.super_path is None or force_superpose:
            #print(self.super_path)
            from superpose import superpose_single
            contents = [self.id, self.name, self.chain]
            self.superpositions = {}
            if references is None:
                references = vars.references

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
                        if "failed_df" in vars:
                            vars.failed_df.loc[len(vars.failed_df)] = [self.id, "gesamt", "gesamt error", "are DISSIMILAR and cannot be reasonably aligned"]
                        contents.extend([99,0])
            if "raw_monomers_df" in vars:
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
            if  data["coverage"] >= (90, 90):
                finalists.append((ref_name,data))
                criteria.append(data["identity"])
            else:
                coverages.append(data["coverage"])
                pass
        if len(finalists) > 0:
            self.super_data = finalists[np.argmax(criteria)]
            self.best_fit = self.super_data[0]
            data = self.super_data[1]
            contents = [self.id, self.super_data[0], round(data["coverage"][0], 1), round(data["coverage"][1],1), data["rmsd"], round(data["identity"]*100)]
            contents = [str(x) for x in contents]
            contents.extend(data["t_matrix"].values())
            #print3(contents)
            if "monomers_df" in vars:
                #print(self.id, type(self.id))
                #print(contents, type(contents))
                #print(vars.monomers_df)
                vars.monomers_df.loc[self.id] = contents
            self.monomers_entries.append(contents)
            self.super_path = data["out_path"]
            self.replaced = PDBParser(QUIET=True).get_structure(self.id, self.super_path).get_list()[1].get_list()[0]
            self.replaced.id = self.chain
            from superpose import create_maps
            self.map_to_ref, self.map_to_ori = create_maps(self.structure, self.replaced, data["map"])
        else:
            self.super_path = None
            if "failed_df" in vars:
                vars.failed_df.loc[len(vars.failed_df)] = [self.id, "monomer", "No coverage", "No reference meets coverage (80%)" + str(coverages)]
            self.failed_entries.append([self.id, "monomer", "No coverage", "No reference meets coverage (80%)" + str(coverages)])




class Mate(BioObject):
    pickle_extension = ".mate"
    pickle_folder = "mates"

    def __init__(self, op_n, operation, params, fixed_monomer, moving_monomer):
        self.structure = None
        self.op_n = op_n
        self.operation = operation
        self.key = None
        self.params = params
        self.positions = {}
        self.fixed_monomer = fixed_monomer
        self.moving_monomer = moving_monomer
        self.update_id()
        self.dimers = []
        self.structures = []
        self.paths = None
        self.c_lines = []
        self.contacts = []

    def update_id(self):
        self.id = "{}_mate_{}_{}_{}".format(self.name, self.fixed_monomer.chain, self.op_n, self.moving_monomer.chain)

    def process(self, parent, min_contacts=0):
        self.name = parent.name
        self.key = parent.key
        self.update_id()

        print2("Processing mate:", self.id)
        self.check_redundancy()
        self.unpack_contacts()
        self.reconstruct_mates(min_contacts=min_contacts)
        #self.export()
        print3(self.dimers)


    def check_redundancy(self):
        if self.fixed_monomer.chain.upper() == self.moving_monomer.chain.upper() and len(list(self.positions.keys())) > 1:
            best_pos = None
            for key, pos in self.positions.items():
                if best_pos is None:
                    #print(pos["value"])
                    best_pos = key, pos
                    continue
                if pos["value"] > best_pos[1]["value"] :
                    #print(pos["value"])
                    best_pos = key, pos
            self.positions = {best_pos[0]: best_pos[1]}


    def unpack_contacts(self):
        from symmetries import convertFromFracToOrth
        for pos in self.positions.values():
            cl = []
            for contact in pos["contacts"]:
                self.contacts.append(contact)
                for c in contact.all_contacts:
                    cl.append(c["line"])
            self.c_lines.extend(cl)


    def reconstruct_mates(self, min_contacts = 0):
        from symmetries import generate_displaced_copy, entity_to_orth, get_operation, print_all_coords
        dimers = []
        import copy
        for position in self.positions.values():
            from symmetries import print_all_coords
            if position["n_contacts"] >= min_contacts:
                #print(position["position"], self.op_n)

                moved_mate = generate_displaced_copy(self.moving_monomer.fractional_structure, distance = position["position"], key =self.key, op_n =self.op_n)
                
                #moved_mate = entity_to_orth(moved_mate, self.params)
                #print_all_coords(moved_mate)
                #self.structures.append(moved_mate)

                moved_monomer = Monomer(self.moving_monomer.name,
                                        self.moving_monomer.chain,
                                        moved_mate,
                                        self,
                                        parent_monomer= self.moving_monomer,
                                        op_n = self.op_n,
                                        position = position,
                                        is_mate = True)

                dimer = Dimer(self.fixed_monomer, moved_monomer)
                dimer.export()
                dimer.pickle()
                dimers.append(dimer)
                
        self.dimers = dimers
        return self.dimers

    def export(self):
        return None
        paths = []
        for structure, position in zip(self.structures, self.positions.values()):
            
            paths.append(super().export(subfolder="mates", in_structure = convertFromFracToOrth(structure.copy(), self.params), extra_id = "_" + clean_string(position["position"], allow = ["-"])))
        self.paths = paths
        return self.paths
    



class Dimer(BioObject):
    pickle_extension = '.dimer'
    pickle_folder = "dimers"

    def __init__(self, monomer1, monomer2, sasa = False):
        self.monomer1 = monomer1
        self.monomer2 = monomer2

        self.position = monomer2.position
        self.op_n = monomer2.op_n

        self.name = monomer1.name
        self.extra_id = monomer2.extra_id
        self.id = "{}_{}{}{}".format(self.name, monomer1.chain, monomer2.chain, self.extra_id, sasa = False)
        print4("Generating dimer:", self.id)

        self.incomplete = True
        self.failed_entries = []
        self.contacts_sasa = []
        self.contacts_symm = []
        self.contacts = []
        self.full_array = None
        self.face1 = None
        self.face2 = None
        self.contact_face1 = None
        self.contact_face2 = None
        self.contacts_faces1 = None
        self.contacts_faces2 = None
        self.m1faces = None
        self.m2faces = None
        self.used_backup = False

        self.com1 = None
        self.com2 = None
        self.pca1 = None
        self.pca2 = None

        self.outer_ids = None
        self.contact_surface = None


        self.process(sasa = sasa)

    def get_contact_surface(self):
        print5("Calculating contact surface")
        from faces import ContactSurface
        if self.best_fit is not None and self.best_fit != "Missmatch":
            self.outer_ids = [ref for ref in vars.references if ref.name == self.best_fit][0].outer_ids
            self.contact_surface = ContactSurface(self.monomer1.replaced, self.monomer2.replaced, outer_ids=self.outer_ids)

    def process(self, sasa = False, pickle=True):
        print4("Processing dimer:", self)
        self.validate()
        if self.incomplete:
            return
        if sasa:
            from surface import get_dimer_sasa
            get_dimer_sasa(self)
        from maths import find_com
        self.com1 = find_com(self.replaced_structure.get_list()[0].get_list()[0].get_atoms())
        self.com2 = find_com(self.replaced_structure.get_list()[0].get_list()[1].get_atoms())
        self.pca1 = self.monomer1.get_monomer_pca()
        self.pca2 = self.monomer2.get_monomer_pca()
        #print(self.com)
        self.get_contact_surface()
        if pickle:
            self.pickle()

    def reprocess(self, by_com = False, contacts=True, faces=True, pickle=True):
        self.process(pickle=pickle)
        if contacts:
            self.get_contacts(force = True)
        if faces:
            self.get_faces(by_com=by_com)
        if pickle:
            self.pickle()



    def get_contacts(self, max_distance= 16, force = False):
        self.contacts = []
        if len(self.contacts) != 0 and not force:
            return
        print2("Calculating contacts")
        from symmetries import Contact
        mon1_atoms = list(self.monomer1.replaced.get_atoms())
        if self.best_fit == "GR":
            from faces import GR_dict
            ref_dict = GR_dict
        else:
            ref_dict = None
        for atom in self.monomer2.replaced.get_atoms():
            #sprint(atom.parent.id[1])
            contact = Contact(atom, self.position, mon1_atoms, max_distance=max_distance, ref_dict = ref_dict)
            #print1(contact.is_contact, contact.shortest_contact)
            if contact.is_contact:
                self.contacts.append(contact)
                #print(contact)

    def get_dihedrals(self, reverse = False):
        #from matplotlib import pyplot as plt
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')

        if not reverse:
            pca_data1 = self.pca1
            pca_data2 = self.pca2
        else:
            pca_data1 = self.pca2
            pca_data2 = self.pca1

        pca1 = pca_data1["pca"]
        pca2 = pca_data2["pca"]

        com1 = pca_data1["com"]
        com2 = pca_data2["com"]

        c0 = add(com1, pca1.components_[0]), add(com2, pca2.components_[0])
        c1 = add(com1, pca1.components_[1]), add(com2, pca2.components_[1])
        c2 = add(com1, pca1.components_[2]), add(com2, pca2.components_[2])

        d = vector(com1, com2)
        dist = length(d)
        #ax.scatter(0,0,0, color = 'yellow')
        #ax.scatter(*com1, color = 'purple')
        d0 = dihedral_angle2(c0[0], com1, com2, c0[1])
        d1 = dihedral_angle2(c1[0], com1, com2, c1[1])
        d2 = dihedral_angle2(c2[0], com1, com2, c2[1])


        #ax.plot(*points_to_line(com1, com2), color="black")
        #ax.plot(*points_to_line(com1, c0[0]), color = "red")
        #ax.plot(*points_to_line(com1, c1[0]), color='green')
        #ax.plot(*points_to_line(com1, c2[0]), color='blue')

        a0 = angle_3_points(c0[0], com1, com2)
        a1 = angle_3_points(c1[0], com1, com2)
        a2 = angle_3_points(c2[0], com1, com2)

        #print(d, c0[0], c1[0], c2[0])
        #print(a0)
        #print(a1)
        #print(a2)

        #ax.set_aspect('equal')
        #plt.show()
        #quit()
        return [d0, d1, d2, a0, a1, a2, dist]


    def _count_contacts(self, backup1=False, backup2=False, use_backups = False):
        print3("Counting Contacts, backup:", backup1, backup2)
        m1faces = None
        m2faces = None
        # print(self.contacts)
        for contact in self.contacts:
            #print(contact)
            #print("normal:", contact.face, contact.face_opposite)
            #print("backup:",contact.backup_face , contact.backup_face_opposite)
            contact.backup_distance = contact.backup_distance - 2
            contact.reprocess_contacts()
            if backup1:
                face_opposite = contact.backup_face_opposite
            else:
                face_opposite = contact.face_opposite

            if backup2:
                face = contact.backup_face
            else:
                face = contact.face

            # print(contact)
            '''if face is None and face_opposite is None:
                contact.reprocess_contacts()'''
            # print_dict(contact.__dict__)

            if face is not None:
                if m2faces is None:
                    m2faces = {face: 1}
                elif face in m2faces.keys():
                    m2faces[face] += 1
                else:
                    m2faces[face] = 1

            if face_opposite is not None:
                if m1faces is None:
                    m1faces = {face_opposite: 1}
                elif face_opposite in m1faces.keys():
                    m1faces[face_opposite] += 1
                else:
                    m1faces[face_opposite] = 1

        # m1 and m2 swapped as face 1 is actually monomer2 # NOT ANYMORE
        need_backup1 = False
        need_backup2 = False
        if m1faces is not None:
            print("#1", m1faces)
            if len(m1faces) > 0:
                self.m1faces = sort_dict(m1faces, as_list=True)
                if len(self.m1faces) > 1:
                    if self.m1faces[0][1] == self.m1faces[1][1]:
                        need_backup1 = True
                self.contact_face1 = self.m1faces[0][0]
            elif not backup1:
                self.contact_face1 = None
                self.m1faces = None
        elif not backup1:
            self.contact_face1 = None
            self.m1faces = None
        print("##1", m1faces)

        if m2faces is not None:
            if len(m2faces) > 0:
                print("#2", m2faces)
                self.m2faces = sort_dict(m2faces, as_list=True)
                if len(self.m2faces) > 1:
                    if self.m2faces[0][1] == self.m2faces[1][1]:
                        need_backup2 = True
                self.contact_face2 = self.m2faces[0][0]
            elif not backup2:
                self.contact_face2 = None
                self.m2faces = None
        elif not backup2:
            self.contact_face2 = None
            self.m2faces = None
        print("##2", m2faces)

        print("self:")
        print("#1#", self.m1faces)
        print("#2#", self.m2faces)

        print("new: (backup: ", backup1, backup2)
        print(m1faces)
        print(m2faces)
        self.used_backup = (backup1 or backup2)
        if (need_backup1 or need_backup2) and not (backup1 or backup2) and use_backups:
            print3("Backup for dimer:", self)
            self._count_contacts(backup1=need_backup1, backup2=need_backup2)
        else:
            print3("Contacts counted for dimer:", self)
            print("#1 {}#".format(backup1), self.m1faces)
            print("#2 {}#". format(backup2), self.m2faces)

    def get_faces(self, by_com = False):
        print2("Identifying faces, by COM: {}".format(by_com))

        #if by_com:
        from faces import get_dimer_faces
        self.face1, self.face2, self.interface_distance = get_dimer_faces(self)

        self._count_contacts()

        face1 = None
        face2 = None

        if not by_com:
            face1 = self.contact_face1
            face2 = self.contact_face2
        else:
            face1 = self.face1
            face2 = self.face2

        if face1 is not None:
            from faces import GR_dict, GR_colours
            self.contacts_faces1 = [GR_colours[face1]]
            for res in GR_dict[face1]:
                self.contacts_faces1.append([self.monomer1.chain, res])

        if face2 is not None:
            from faces import GR_dict, GR_colours
            self.contacts_faces2 = [GR_colours[face2]]
            for res in GR_dict[face2]:
                self.contacts_faces2.append([self.monomer2.chain, res])





    def validate(self):
        if self.monomer1.best_fit == self.monomer2.best_fit:
            self.best_fit = self.monomer1.best_fit
        else:
            self.best_fit = "Missmatch"
        self.best_match = None
        print5("Best fit:", self.best_fit)

        if self.monomer1.super_path is None or self.monomer2.super_path is None:
            if "failed_df" in vars:
                vars.failed_df.loc[len(vars.failed_df)] = [self.id, "dimer", "Missing superposition", "At least one superposition is missing, {}:{}, {}:{}".format(self.monomer1.chain, self.monomer1.super_path,self.monomer2.chain,self.monomer2.super_path)]
            self.failed_entries.append([self.id, "dimer", "Missing superposition", "At least one superposition is missing, {}:{}, {}:{}".format(self.monomer1.chain,self.monomer1.super_path,self.monomer2.chain,self.monomer2.super_path)])
            self.incomplete = True
        else:
            self.incomplete = False
            self.original_structure, self.replaced_structure, self.merged_structure = self.merge_structures()
            from faces import get_pca
            from maths import find_com
            self.com = find_com(self.replaced_structure.get_atoms())
            self.pca = get_pca(self.replaced_structure)

            #_, _, self.merged_structure = self.merge_structures()
            #[print(chain) for chain in self.merged_structure.get_chains()]

        
    def merge_structures(self):

        monomer1 = self.monomer1
        monomer2 = self.monomer2

        #print1("({}) Merging monomers:".format(self.id), monomer1, monomer1.structure.id, monomer2, monomer2.structure.id)
        #print(monomer1.super_path, monomer2.super_path)
        structureA = PDBParser(QUIET=True).get_structure(monomer1.id, monomer1.super_path)
        structureB = PDBParser(QUIET=True).get_structure(monomer2.id, monomer2.super_path)

        originalA = structureA[0]
        originalA.get_list()[0].id = originalA.get_list()[0].id.upper()
        originalB = structureB[0]
        originalB.get_list()[0].id = originalB.get_list()[0].id.lower()
        replacedA = structureA[1]
        replacedB = structureB[1]
        originalA.id = originalB.id = 0
        replacedA.id = replacedB.id = 1
        #print(originalA.get_list()[0], originalB.get_list()[0], replacedA.get_list()[0], replacedB.get_list()[0])
        
        replacedA.get_list()[0].id = originalA.get_list()[0].id
        replacedB.get_list()[0].id = originalB.get_list()[0].id
        
        
        #print(originalA.get_list()[0], originalB.get_list()[0], replacedA.get_list()[0], replacedB.get_list()[0])
        

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
        
        #print(merged_structure)
        #[print(model) for model in merged_structure.get_models()]
        #[print(chain) for chain in merged_structure.get_chains()]

        
        #monomer1.replaced = PDBParser(QUIET=True).get_structure(monomer1.id+"_r", monomer1.super_path)[1].get_list()[0]
        #monomer2.replaced = PDBParser(QUIET=True).get_structure(monomer2.id+"_r", monomer1.super_path)[1].get_list()[0]

        #monomer1.replaced.detach_parent()
        #monomer2.replaced.detach_parent()
        

        assert len([1 for _ in merged_structure.get_chains()]) == 4   

        
        return original_structure, replaced_structure, merged_structure



    def export(self):
        if not self.incomplete:
            #self.original_path = super().export(subfolder="dimers_original", in_structure=self.original_structure, extra_id = "_dimer_ori")
            self.replaced_path = super().export(subfolder="dimers_replaced", in_structure=self.replaced_structure, extra_id = "_dimer_rep")
            self.merged_path = super().export(subfolder="dimers_merged", in_structure=self.merged_structure, extra_id = "_dimer_merged")
            return self.merged_path
        return None


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

    def __init__(self, path):

        self.o_path = path
        self.parse_structure(parse_original=True, calculate_sasa=True, n_points=200, radius=6)
        self.outer_ids = self.get_outer_res_list()
        self.name = os.path.basename(path).split(".")[0]
        self.best_fit = self.name
        self.structure = [chain for chain in self.structure.get_chains()][0]
        self.replaced = self.structure
        self.chain = self.structure.id
        self.structure.id = self.chain
        self.id = "reference_{}_{}".format(self.name, self.chain)
        self.is_reference = True
        self.export()

        self.raw_monomers_entries = []
        self.failed_entries = []
        self.monomers_entries = []

        self.sequence = None
        self.face_dict = None
        if self.name == "GR":
            self.face_dict = None
        from faces import get_pca, get_terminals, find_com, get_face_coms
        self.terminals = get_terminals(self.structure)
        #self.pca = get_pca(self.structure)
        self.com = find_com(self.structure.get_atoms())
        if self.name == "GR":
            self.face_coms = get_face_coms(self)
        self.pca = dict(pca=get_pca(self.structure, com=self.com, closer_to="N"),
                        com=self.com,
                        chain=self.chain)


    def reshape_face_dict(self):
        from faces import GR_dict
        self.face_dict = GR_dict.copy()

    def get_outer_res_list(self, threshold=10, inner=False, id_only=True, complete_list = False):
        res_list = []
        """[print(res.sasa) for res in self.structure.get_residues()]
        for res in self.structure.get_residues():
            res.get_list()[0].bfactor = res.sasa
        from pyMol import pymol_temp_show, pymol_start, pymol_save_temp_session, pymol_open_session_terminal
        #pymol_start(show=False)
        pymol_temp_show(self.structure)
        pymol_colour("rainbow", spectrum="b")
        session = pymol_save_temp_session()
        pymol_open_session_terminal(session)
        quit()"""

        for residue in self.structure.get_residues():
            is_outer = False
            if residue.sasa >= threshold:
                is_outer = True
            if inner:
                is_outer = not is_outer
            if is_outer:
                if id_only:
                    res_list.append(residue.id[1])
                else:
                    res_list.append(residue)
            elif complete_list:
                res_list.append(None)
        return res_list