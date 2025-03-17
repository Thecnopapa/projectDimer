import os


from utilities import *
from Globals import root, local, vars
import numpy as np
import pandas as pd

import Bio.PDB
from Bio.PDB import SASA




def get_monomer_sasa(self, n_points =100, radius = 1.6):
    if self.best_fit is None:
        return None
    sr = SASA.ShrakeRupley(n_points=n_points, probe_radius=radius)
    sr.compute(self.structure, level="R")
    sasas = []
    for res in self.structure.get_residues():
        print(res.__dict__)
        sasas.append(res.sasa)
    self.sasas = sasas.copy()


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

