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
            print(self.id, "has no replaced_structure")
            return None
        print(self.replaced_structure.get_list()[0].get_list())
        structure1 = self.monomer1.replaced
        structure2 = self.monomer1.replaced
        dimer_structure = self.replaced_structure.get_list()[0]
        print(structure1, structure2, dimer_structure)
        print(structure1.get_full_id(), structure2.get_full_id(), dimer_structure.get_full_id())

    if self.monomer1.super_data[0] == self.monomer2.super_data[0]:
        self.best_fit = self.monomer1.super_data[0]
        sr = SASA.ShrakeRupley(n_points=n_points, probe_radius=radius)


        sr.compute(structure1,level="R")
        sr.compute(structure2,level="R")
        sasas1 = []
        sasas2 = []
        sasas1D = []
        sasas2D = []
        for res1 in structure1.get_list():
            sasas1.append(res1.sasa)
        #print(sasas1)
        for res2 in structure2.get_list():
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

        sasa_df = pd.DataFrame(columns=["ID", "name", "sasa1", "sasa2", "sasa1_dimer", "sasa2_dimer", "diff1", "diff2"])


        print(len(sasas1), len(sasas2), len(sasas1D), len(sasas2D))
        for i, (sasa1, sasa2, sasa1D,sasa2D) in enumerate(zip(sasas1, sasas2, sasas1D, sasas2D)):
            #print(i, sasa1, sasa2, sasa1D, sasa2D)
            sasa_df.loc[i] = [structure1.get_list()[i].id[1], structure1.get_list()[i].resname, sasa1, sasa2, sasa1D, sasa2D, sasa1-sasa1D, sasa2-sasa2D]
            #print(sasa_df.loc[i])

    else:
        print("Best fits are not the same:",self.monomer1.super_data[0], self.monomer1.super_data[0])
        self.sasa_df = None
        self.best_fit = None
        return None

    if len(sasa_df) > 0:
        self.sasa_df = sasa_df
        return sasa_df
    else:
        self.sasa_df = None
        return None



if __name__ == "__main__":

    import globals

    globals.set_root("../")
    if os.name == "nt":
        globals.set_local("C:/Users/iainv/localdata/_local/projectB")
    elif os.name == "posix":
        globals.set_local("/localdata/iain/_local/projectB")
    from globals import root, local, vars


    FORCE_SASA = False

    tprint("Loading dimers")
    from imports import load_dimers
    dimers = load_dimers()
    eprint("Dimers loaded")



    tprint("Calculating SASAs")
    progress = ProgressBar(len(dimers))
    for dimer in dimers:
        if "sasa_df" not in dimer.__dict__.keys() or FORCE_SASA:
            get_contact_res(dimer)
        progress.add(info = dimer.id)
    eprint("SASAs calculated")

    from imports import pickle
    print("Pickling...")
    pickle(dimers)
    print("Done")


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
                    print(dimer.sasa_df["diff2"])'''
            try:
                ref_df.at[i, dimer.monomer1.id] = dimer.sasa_df["diff1"].tolist()
            except:
                pass
            try:
                ref_df.at[i, dimer.monomer2.id] = dimer.sasa_df["diff2"].tolist()
            except:
                pass
            progress.add(info=dimer.id)

        print(ref_df)
        ref_df.to_csv(os.path.join(root.dataframes, "{}_sasas.csv".format(reference.name)))
    eprint("SASAs merged")


    from imports import pickle
    print("Pickling...")
    pickle(dimers)
    print("Done")