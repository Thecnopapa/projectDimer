import os
from email.encoders import encode_noop

from utilities import *
from Globals import root, local, vars
import numpy as np
import pandas as pd

#import Bio.PDB
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB import SASA, Structure, Model


def get_sasa_contacts(self):
    pass

def get_monomer_sasa(self, n_points =100, radius = 1.6, use_replaced = True):
    print2("Calculating SASA for {}, replaced = ".format(self), use_replaced)
    if use_replaced:
        structure = self.replaced
        #print(structure)
        #print(structure.parent)
        structure.detach_parent()
        #print(structure.parent)
        from maths import find_com
        #print(find_com(structure.get_atoms()))
    else:
        print("why not?")
        quit()
        structure = self.structure
    if self.best_fit is None:
        return None

    sr = SASA.ShrakeRupley(n_points=n_points, probe_radius=radius)
    sr.compute(structure, level="A")
    sasas = []
    for atom in structure.get_atoms():
        sasas.append(round(atom.sasa,4))
    from pyMol import pymol_temp_show, pymol_start, pymol_colour, pymol_reset
    #pymol_start(show=True)

    for atom in structure.get_atoms():
        atom.bfactor = atom.sasa

    #pymol_temp_show(structure)

    print(sum(sasas))
    if self.sasas is None:
        self.sasas = sasas
    else:
        self.sasas.extend(sasas)
        quit()

    #print(len(list(structure.get_atoms())))



def get_dimer_sasa(self, n_points =100, radius = 1.6, use_replaced = True):
    print3("Calculating SASA for {}, replaced = ".format(self), use_replaced)

    if self.monomer1.best_fit != self.monomer2.best_fit or self.monomer1.best_fit is None:
        print1("{}: Best fits are not the same: {} / {}".format(self.id, self.monomer1.best_fit, self.monomer2.best_fit,
                                                                "\n"))
        self.sasa_df = None
        return None
    from pyMol import pymol_temp_show, pymol_start, pymol_colour, pymol_reset
    #pymol_start(show=True)

    #pymol_temp_show(self.monomer1.replaced)
    #pymol_temp_show(self.monomer2.replaced)
    from maths import find_com
    if use_replaced:
        structure = self.replaced_structure.copy()
        structure.detach_parent()
        #print(structure)
        #print(structure.get_list()[0].get_list())
        #print([find_com(chain.get_atoms()) for chain in structure.get_chains()])
    else:
        print("why not?")
        quit()
        structure = self.original_structure.copy()

    assert len(list(structure.get_chains())) == 2

    sr = SASA.ShrakeRupley(n_points=n_points, probe_radius=radius)
    sasas1D = []
    sasas2D = []

    for res in structure.get_residues():
        res.sasa = None
    chain1, chain2 = list(structure.get_chains())
    #print(chain1,chain2)
    chain1 = chain1.copy()
    chain2 = chain2.copy()
    chain1.detach_parent()
    chain2.detach_parent()

    new_structure = Structure.Structure("temp_dimer_{}".format(self.id))
    new_structure.add(Model.Model(0))
    new_model = new_structure[0]
    #print(new_structure)
    #print(new_model)
    new_model.add(chain1)
    new_model.add(chain2)


    #print(structure)
    #print(self.monomer1, self.monomer1.parent, chain1.id, "//", self.monomer2, self.monomer2.parent_monomer, chain2.id)
    #print("COM at sasa", find_com(self.monomer1.replaced.get_atoms()), find_com(self.monomer2.replaced.get_atoms()))
    #print("COM at sasa (saved)",self.monomer1.com, self.monomer2.com)
    '''sr.compute(chain1, level="R")
    for res in chain1.get_residues():
        sasas1D.append(res.sasa)
    sr.compute(chain2, level="R")
    for res in chain2.get_residues():
        sasas2D.append(res.sasa)'''
    #print(len(list(chains[0].get_residues())),  len(list(chains[1].get_residues())))

    sr.compute(new_structure, level="A")
    #print()
    for atom in new_structure.get_atoms():
        atom.bfactor = atom.sasa

    #pymol_temp_show(chain1)
    #pymol_temp_show(chain2)
    #print(chain1.parent)
    #print(chain2.parent)
    #pymol_temp_show(self.monomer1.replaced)
    #pymol_temp_show(self.monomer2.replaced)

    #pymol_colour("red_yellow_green",spectrum="b" )
    #input()
    #pymol_reset()
    assert len(list(chain1.get_atoms())) == len(list(chain2.get_atoms()))
    #print(len(list(new_structure.get_atoms())))
    for res1D, res2D in zip(chain1.get_atoms(), chain2.get_atoms()):
        #print(res1D.get_full_id(), res2D.get_full_id(), res1D.sasa, res2D.sasa)
        sasas1D.append(round(res1D.sasa,4))
        sasas2D.append(round(res2D.sasa,4))
        #print(res1D.sasa, res2D.sasa)

    print4(sum(sasas1D))
    print4(sum(sasas2D))
    self.sasas1D = sasas1D
    self.sasas2D = sasas2D




def build_contact_arrays(self, sasa = True):

    ################### SASA contacts ##############################
    if sasa:
        sasas1 = self.monomer1.sasas
        sasas2 = self.monomer2.sasas
        sasas1D = self.sasas1D
        sasas2D = self.sasas2D
        print(self.monomer1.chain)
        print(self.monomer2.chain)
        residues = list(self.monomer1.replaced.get_residues())
        sasa_array = []
        print(len(residues))
        print(len(sasas1))
        assert len(residues) == len(list(sasas1))
        assert len(residues) == len(list(sasas2))
        assert len(residues) == len(list(sasas1D))
        assert len(residues) == len(list(sasas2D))

        for i, res in enumerate(residues):
            #print(sasas1[i], sasas2[i], sasas1D[i], sasas2D[i], sasas1[i]-sasas1D[i], sasas2[i]-sasas2D[i])
            if bool(sasas1[i]-sasas1D[1] > 0):
                sasa_array.append([self.monomer1.chain, res.id[1]])
                #print(sasa_array[-1], sasa1, sasa1D,"\t", sasa1-sasa1D)
            if bool(sasas2[i]-sasas2D[0] > 0):
                sasa_array.append([self.monomer2.chain,res.id[1]])
                #print(sasa_array[-1], sasa2, sasa2D, "\t", sasa2-sasa2D)


        '''for sasa1, sasa2, sasa1D, sasa2D, residue in zip(sasas1, sasas2, sasas1D, sasas2D, residues):
            #print([type(s) for s in (sasa1, sasa2, sasa1D, sasa2D)])
            print(sasa1, sasa1D, "\t", sasa2, sasa2D, "\t", round(sasa1-sasa1D, 2), round(sasa2-sasa2D, 2))
            if bool(sasa1-sasa1D > 0):
                sasa_array.append([self.monomer1.chain, residue.id[1]])
                #print(sasa_array[-1], sasa1, sasa1D,"\t", sasa1-sasa1D)
            if bool(sasa2-sasa2D > 0):
                sasa_array.append([self.monomer2.chain,residue.id[1]])
                #print(sasa_array[-1], sasa2, sasa2D, "\t", sasa2-sasa2D)'''

        self.contacts_sasa = sasa_array
        print2("Number of contacts by SASA:", len(sasa_array))

    ############### Symmetry contacts ##############################
    print(self.contacts)
    contact_array = []
    self.c_lines = []
    empty_array = {res.id[1]:[False, False] for res in self.monomer1.replaced.get_residues()}
    for contact in self.contacts:
        contact_array.append([self.monomer2.chain, contact.atom.get_full_id()[-2][1]])
        #print(contact_array[-1])
        sc = contact.shortest_contact
        cl = [self.monomer1.chain, sc["target_atom"].get_full_id()[-2][1]]
        if cl not in contact_array:
            contact_array.append(cl)
        self.c_lines.append(sc["line"])
        print(contact_array[-1])
        empty_array[contact.atom.parent.id[1]][1] = True
        for c in contact.all_contacts:
            empty_array[c["target_atom"].parent.id[1]][0] = True

    #print(contact_array)
    #[print(row) for row in empty_array.items()]
    print2("Number of contacts by symmetry:", len(contact_array))
    self.contacts_symm = contact_array

    if self.best_fit is not None and self.best_fit != "Missmatch":
        full_array = [value for value in empty_array.values()]
        contact_df = vars["clustering"]["contacts"][self.best_fit]
        contact_df[self.id] = full_array
        print(contact_df)







########## OLD #########################################################################################################

def get_contact_res(self, use_replaced = True, radius=1.6, n_points=100):


    if use_replaced:
        if not "replaced_structure" in self.__dict__.keys():
            print1("{}: has no replaced_structure".format(self.id))
            return None
        #print(self.replaced_structure.get_list()[0].get_list())
        structure1 = self.monomer1.replaced
        structure2 = self.monomer1.replaced
        dimer_structure = self.replaced_structure.get_list()[0]
        #print(structure1, structure2, dimer_structure)
        #print(structure1.get_full_id(), structure2.get_full_id(), dimer_structure.get_full_id())

    if self.monomer1.best_fit == self.monomer2.best_fit and self.monomer1.best_fit is not None:
        sr = SASA.ShrakeRupley(n_points=n_points, probe_radius=radius)

        sasas1 = []
        sasas2 = []
        sasas1D = []
        sasas2D = []

        sr.compute(dimer_structure.get_list()[0],level="R")
        for res1 in dimer_structure.get_list()[0].get_list():
            sasas1.append(res1.sasa)
        #print(sasas1)

        sr.compute(dimer_structure.get_list()[1],level="R")
        for res2 in dimer_structure.get_list()[1].get_list():
            sasas2.append(res2.sasa)
        #print(sasas2)

        sr.compute(dimer_structure,level="R")
        for res1D, res2D in zip(dimer_structure.get_list()[0].get_list(), dimer_structure.get_list()[1].get_list()):
            sasas1D.append(res1D.sasa)
            sasas2D.append(res2D.sasa)

        #print(sasas1)
        #print(sasas1D)
        #print(sasas2)
        #print(sasas2D)

        sasa_df = pd.DataFrame(columns=["ID", "name", "sasa1", "sasa2", "sasa1_dimer", "sasa2_dimer", "diff1", "diff2", "is_contact1", "is_contact2"])


        #print(len(sasas1), len(sasas2), len(sasas1D), len(sasas2D))
        for i, (sasa1, sasa2, sasa1D,sasa2D) in enumerate(zip(sasas1, sasas2, sasas1D, sasas2D)):
            #print(i, sasa1, sasa2, sasa1D, sasa2D)
            sasa_df.loc[i] = [self.replaced_structure.get_list()[0].get_list()[0].get_list()[i].id[1], structure1.get_list()[i].resname, sasa1, sasa2, sasa1D, sasa2D, sasa1-sasa1D, sasa2-sasa2D,(sasa1-sasa1D) >0, (sasa2-sasa2D) >0 ]
            #print(sasa_df.loc[i])

    else:
        print1("{}: Best fits are not the same: {} / {}".format(self.id, self.monomer1.best_fit, self.monomer2.best_fit, "\n"))
        self.sasa_df = None
        return None

    if len(sasa_df) > 0:
        self.sasa_df = sasa_df
        if True in sasa_df["is_contact1"].values or True in sasa_df["is_contact2"].values:
            self.no_contacts = False
            print1("{}: No contacts found".format(self.id))
        else:
            self.no_contacts = True
        return sasa_df

    else:
        self.sasa_df = None
        return None


def import_X_df(path, name):
    df = pd.read_csv(path, header=0)
    df.replace("X", True, inplace=True)
    df.fillna(False, inplace=True)
    #print(df)
    new_path = os.path.join(root.dataframes, name+".csv")
    df.to_csv(new_path, index=False)
    return new_path








def surface(FORCE_SASA = True, FORCE_SIMILARITY = True, BALL_SIZE = 1.6):


    tprint("Loading dimers")
    from imports import load_dimers
    dimers = load_dimers()
    eprint("Dimers loaded")

    from dataframes import load_failed_dfs
    load_failed_dfs()

    tprint("Calculating SASAs")
    progress = ProgressBar(len(dimers))
    for dimer in dimers:
        if "sasa_df" not in dimer.__dict__.keys() or FORCE_SASA:
            get_contact_res(dimer, radius=BALL_SIZE)
        progress.add(info = dimer.id)
    eprint("SASAs calculated")


    sprint("loading references")
    from imports import load_from_files
    references = load_from_files(root.references,
                                 pickle_extension=".reference",
                                 is_reference=True,
                                 ignore_selection=True,
                                 force_reload=False)


    tprint("Merging SASAs")
    for reference in references:
        print("Reference:", reference.name)
        ref_data={"ResNum": [],
                  "ResName": []}
        for res in reference.structure.get_list():
            ref_data["ResNum"].append(res.id[1])
            ref_data["ResName"].append(res.resname)
        ref_df = pd.DataFrame(data=ref_data)

        progress = ProgressBar(len(dimers))
        for dimer in dimers:
            if dimer.best_fit == reference.name:
                if not dimer.no_contacts:
                    print(len(ref_df),len(dimer.sasa_df["diff1"]), len(dimer.sasa_df["diff2"]))
                    ref_df[dimer.id] = list(zip(dimer.sasa_df["is_contact1"], dimer.sasa_df["is_contact2"]))
                else:
                    print1("{} has no contacts". format(dimer.id))
            progress.add(info=dimer.id)

        print(ref_df)
        ref_df.to_csv(os.path.join(root.dataframes, "{}_sasas.csv".format(reference.name)))
    eprint("SASAs merged")


    tprint("Comparing contacts")
    classified_df = pd.DataFrame(columns = ["ID","Best_Fit", "Best_Match", "Similarity", "Inverse"])
    for sasa_df in os.listdir(root.dataframes):
        if "_sasas" in sasa_df:
            ref_name = sasa_df.split("_sasas")[0]
            print1("Reference:", ref_name, "({})". format(sasa_df))

            ref_df_name = "{}_reference_clusters.csv".format(ref_name)
            ref_df_path = os.path.join(root.dataframes, ref_df_name)
            if not ref_df_name in os.listdir(root.dataframes):
                try:
                    import_X_df(os.path.join(root.dataframes, "{}_reference_clusters_raw.csv".format(ref_name)),
                                "{}_reference_clusters".format(ref_name))
                    print2("Reference df:", ref_df_path)
                except:
                    print2("Reference df (raw or processed) not found")
                    print3(ref_df_name)
                    print3("processed:", ref_df_path)
                    print3("raw:", os.path.join(root.dataframes, "{}_reference_clusters_raw.csv".format(ref_name)))
                    continue
            ref_df = pd.read_csv(ref_df_path)

            sasa_df = pd.read_csv(os.path.join(root.dataframes, sasa_df), index_col=0)
            n_dimers = len(sasa_df.columns) - 2
            print2("Number of dimers:", n_dimers)

            groups = []
            n_groups = (len(ref_df.columns) - 1) / 2
            assert n_groups % 2 == 0
            ref_nums = ref_df["ResNum"]

            progress = ProgressBar(n_dimers)
            from maths import difference_between_boolean_pairs
            for c in range(n_dimers):
                dimer_sasas = sasa_df.iloc[:, [0, c + 2]]
                dimer_id = sasa_df.columns[c+2]
                similarities = []
                for group in range(int(n_groups)):
                    ref_sasaA = ref_df.iloc[:, [0, group * 2 + 1]]
                    ref_sasaB = ref_df.iloc[:, [0, group * 2 + 2]]
                    total_len = len(ref_sasaA)

                    diffX = [0, 0]
                    diffx = [0, 0]

                    for ref_num in ref_nums:
                        print(group, ref_num, end = "\r")
                        #print(dimer_sasas.loc[dimer_sasas["ResNum"] == ref_num].values)
                        sA, sB = clean_list([dimer_sasas.loc[dimer_sasas["ResNum"]==ref_num].values[0,1]], delimiter=",", format="bool")
                        rsA = ref_sasaA.loc[ref_sasaA["ResNum"] == ref_num].iloc[:,1].values[0]
                        rsB = ref_sasaB.loc[ref_sasaB["ResNum"] == ref_num].iloc[:,1].values[0]

                        resX, resx = difference_between_boolean_pairs(sA, sB, rsA, rsB)
                        diffX[0] += resX[0]
                        diffX[1] += resX[1]
                        diffx[0] += resx[0]
                        diffx[1] += resx[1]

                    if diffX[0] != 0:
                        diffX = diffX[0] / diffX[1]
                    else:
                        diffX = 0
                    if diffx[0] != 0:
                        diffx = diffx[0] / diffx[1]
                    else:
                        diffx = 0


                    inverse = False
                    if diffx > diffX:
                        inverse = True
                    similarities.append((group+1, max([diffX,diffx]), inverse))
                    print1(similarities[-1][0], round(similarities[-1][1],2))
                best_match = max(similarities, key= lambda x: x[1])
                print("Best match for {}: {}, with {}% similarity, inverse: {}\n".format(dimer_id, best_match[0],
                                                                                         round(100 * best_match[1]),
                                                                                         best_match[2]))

                classified_df.loc[len(classified_df)] = [dimer_id, ref_name, best_match[0],
                                                         round(best_match[1] * 100), best_match[2]]
                if "dimers" in locals() or "dimers" in globals():
                    for dimer in dimers:
                        if dimer.id == dimer_id:
                            dimer.best_match = best_match
                progress.add(info=dimer_id)
            classified_df.to_csv(os.path.join(root.dataframes, "classified_df.csv"))
    eprint("Contacts compared")

    if "dimers" in locals() or "dimers" in globals():
        from imports import pickle
        print("Pickling...")
        pickle(dimers)
        print("Done")











if __name__ == "__main__":

    import setup

    from Globals import root, local, vars

    '''#### BALL SIZE TESTING ###
    ball_sizes = [1.3, 1.4, 1.5, 1.6,1.7, 1.8, 1.9, 2.0]
    for n in ball_sizes:
        vars["BALL_SIZE"] = n
        surface(FORCE_SASA=True, FORCE_SIMILARITY=True, BALL_SIZE=n)
        from clustering import clustering
        clustering(FORCE_ALL = True)
    #### BALL SIZE TESTING ###'''

    surface(FORCE_SASA=True, FORCE_SIMILARITY=True, BALL_SIZE =1.6)

