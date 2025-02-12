import os


from utilities import *
from globals import root, local, vars
import numpy as np
import pandas as pd

import Bio.PDB
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, StructureBuilder, Structure, SASA


def get_contact_res(self, use_replaced = True, radius=1.4, n_points=100):


    if use_replaced:
        if not "replaced_structure" in self.__dict__.keys():
            print1(self.id, "has no replaced_structure")
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
        print(sasas1)

        sr.compute(dimer_structure.get_list()[1],level="R")
        for res2 in dimer_structure.get_list()[1].get_list():
            sasas2.append(res2.sasa)
        print(sasas2)

        sr.compute(dimer_structure,level="R")
        for res1D, res2D in zip(dimer_structure.get_list()[0].get_list(), dimer_structure.get_list()[1].get_list()):
            sasas1D.append(res1D.sasa)
            sasas2D.append(res2D.sasa)
        #print(sasas1)
        #print(sasas1D)
        #print(sasas2)
        #print(sasas2D)

        sasa_df = pd.DataFrame(columns=["ID", "name", "sasa1", "sasa2", "sasa1_dimer", "sasa2_dimer", "diff1", "diff2", "is_contact1", "is_contact2"])


        print(len(sasas1), len(sasas2), len(sasas1D), len(sasas2D))
        for i, (sasa1, sasa2, sasa1D,sasa2D) in enumerate(zip(sasas1, sasas2, sasas1D, sasas2D)):
            #print(i, sasa1, sasa2, sasa1D, sasa2D)
            sasa_df.loc[i] = [self.replaced_structure.get_list()[0].get_list()[0].get_list()[i].id[1], structure1.get_list()[i].resname, sasa1, sasa2, sasa1D, sasa2D, sasa1-sasa1D, sasa2-sasa2D,(sasa1-sasa1D) >0, (sasa2-sasa2D) >0 ]
            #print(sasa_df.loc[i])

    else:
        print("Best fits are not the same:",self.monomer1.best_fit, self.monomer1.best_fit)
        self.sasa_df = None
        return None

    if len(sasa_df) > 0:
        self.sasa_df = sasa_df
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

def compare_to_reference(sasa_df_path, reference_df_path):
    sasa_df = pd.read_csv(sasa_df_path, index_col=0)
    reference_df = pd.read_csv(reference_df_path,index_col=0)
    iterate_over_dimers(sasa_df, start_col=2)

def iterate_over_dimers(df, start_col = 2):
    resnums = df["ResNum"]
    df_dimersA = df.iloc[:,start_col::2]
    df_dimersB = df.iloc[:,(start_col+1)::2]
    for i in range(len(df_dimersB.columns)):
        dimers.append(df_dimersA.iloc[:,i], df_dimersB.iloc[:,i])



















if __name__ == "__main__":

    import globals

    globals.set_root("../")
    if os.name == "nt":
        globals.set_local("C:/Users/iainv/localdata/_local/projectB")
    elif os.name == "posix":
        globals.set_local("/localdata/iain/_local/projectB")
    from globals import root, local, vars


    FORCE_SASA = False
    FORCE_SIMILARITY = False

    tprint("Loading dimers")
    from imports import load_dimers
    dimers = load_dimers()
    eprint("Dimers loaded")



    tprint("Calculating SASAs")
    progress = ProgressBar(len(dimers))
    for dimer in dimers:
        if "sasa_df" not in dimer.__dict__.keys() or FORCE_SASA:
            get_contact_res(dimer)
            print("")
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
    progress = ProgressBar(len(dimers))
    for reference in references:
        print("Reference:", reference.name)
        ref_data={"ResNum": [],
                  "ResName": []}
        for res in reference.structure.get_list():
            ref_data["ResNum"].append(res.id[1])
            ref_data["ResName"].append(res.resname)
        ref_df = pd.DataFrame(data=ref_data)

        for i,dimer in enumerate(dimers):
            '''print(dimer.__dict__.keys())
            if "best_fit" in dimer.__dict__.keys():
                print(dimer.best_fit == reference.id)
                if dimer.best_fit == reference.id:
                    print(dimer.sasa_df["diff1"])
                    print(dimer.sasa_df["diff2"])
            try:
                ref_df.at[i, dimer.monomer1.id] = dimer.sasa_df["diff1"].tolist()
            except:
                pass
            try:
                ref_df.at[i, dimer.monomer2.id] = dimer.sasa_df["diff2"].tolist()
            except:
                pass
            '''
            if dimer.best_fit == reference.name:
                print(len(ref_df),len(dimer.sasa_df["diff1"]), len(dimer.sasa_df["diff2"]))
                ref_df[dimer.monomer1.id] = dimer.sasa_df["is_contact1"]
                ref_df[dimer.monomer2.id] = dimer.sasa_df["is_contact2"]
            progress.add(info=dimer.id)

        print(ref_df)
        ref_df.to_csv(os.path.join(root.dataframes, "{}_sasas.csv".format(reference.name)))
    eprint("SASAs merged")



    tprint("Comparing contacts")
    progress = ProgressBar(len(dimers))
    for dimer in dimers:
        progress.add(info = dimer.id)
        if dimer.best_fit != "GR":
            continue
        if "best_match" in dimer.__dict__.keys() and not FORCE_SIMILARITY:
            continue
        sprint(dimer.id, dimer.best_fit)
        ref_df_name = "{}_reference_clusters.csv".format(dimer.best_fit)
        ref_df_path = os.path.join(root.dataframes, ref_df_name)
        if not ref_df_name in os.listdir(root.dataframes):
            try:
                import_X_df(os.path.join(root.dataframes,"{}_reference_clusters_raw.csv".format(dimer.best_fit)),
                            "{}_reference_clusters".format(dimer.best_fit))
            except:
                sprint("Reference df (raw or processed) not found")
                print1(ref_df_name)
                print2(ref_df_path)
                print2(os.path.join(root.dataframes,"{}_reference_clusters_raw.csv".format(dimer.best_fit)))
                continue

        ref_df = pd.read_csv(ref_df_path)

        groups =[]
        #print(ref_df.columns)
        #print(len(ref_df.columns))
        n_groups = (len(ref_df.columns)-1) / 2
        #print(n_groups % 2)
        assert n_groups % 2 == 0

        sasaA = dimer.sasa_df[["ID", "is_contact1"]]
        sasaB = dimer.sasa_df[["ID", "is_contact2"]]

        similarities = []

        ref_nums = ref_df["ResNum"]
        for group in range(int(n_groups)):
            ref_sasaA = ref_df.iloc[:,[0,group*2+1]]
            ref_sasaB = ref_df.iloc[:,[0,group*2+2]]

            total_len= len(ref_sasaA)

            AA = [0,0]
            AB = [0,0]
            BA = [0,0]
            BB = [0,0]
            for ref_num in ref_nums:
                print(group, ref_num, end = "\r")
                sA = sasaA.loc[sasaA["ID"]==ref_num].iloc[:,1].values[0]
                sB = sasaB.loc[sasaB["ID"]==ref_num].iloc[:,1].values[0]
                rsA = ref_sasaA.loc[ref_sasaA["ResNum"] == ref_num].iloc[:,1].values[0]
                rsB = ref_sasaB.loc[ref_sasaB["ResNum"] == ref_num].iloc[:,1].values[0]
                if sA and rsA:
                    AA[0] +=1
                elif sA or rsA:
                    AA[1] +=1
                if sB and rsB:
                    BB[0] +=1
                elif sB or rsB:
                    BB[1] +=1
                if sA and rsB:
                    AB[0] +=1
                elif sA or rsB:
                    AB[1] +=1
                if sB and rsA:
                    BA[0] +=1
                elif sB or rsA:
                    BA[1] +=1

            per1 = (AA[0]/(AA[0]+AA[1]) + BB[0]/(BB[0]+BB[1]))/2
            per2 = (AB[0]/(AA[0]+AA[1]) + BA[0]/(BB[0]+BB[1]))/2
            inverse = False
            if per2 > per1:
                inverse = True
            similarities.append((group+1, max([per1,per2]), inverse))
            print1(similarities[-1][0], round(similarities[-1][1]*100) ,"%")
        best_match = max(similarities, key= lambda x: x[1])
        print("Best match for {}: {}, with {}% similarity, inverse: {}\n".format(dimer.id, best_match[0], round(100*best_match[1]), best_match[2]))
        dimer.best_match = best_match
    eprint("Contacts compared")


    from imports import pickle
    print("Pickling...")
    pickle(dimers)
    print("Done")

    tprint("Classification to df")
    classified_df = pd.DataFrame(columns = ["ID","Best_Fit", "Best_Match", "Similarity", "Inverse"])
    for dimer in dimers:
        if "best_match" in dimer.__dict__.keys():
            classified_df.loc[len(classified_df)] = [dimer.id, dimer.best_fit, dimer.best_match[0], round(dimer.best_match[1]*100),dimer.best_match[2]]
    classified_df.to_csv(os.path.join(root.dataframes, "classified_df.csv"))
    eprint("Classification to df done")



