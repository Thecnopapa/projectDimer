import os

from utilities import *
from Globals import root, local, vars
import numpy as np
import pandas as pd




def save_dfs(force = False, general= True, clustering = False):
    sprint("Saving dataframes...")
    if "do_only" in vars:
        if len(vars.do_only) == 0 or vars.do_only is None or force:
            for key, value in vars.copy().items():
                if "df" in key and general:
                    print1("Saving {}.csv".format(key))
                    value.to_csv(os.path.join(root.dataframes,f"{key}.csv"), header = True, index=False)
                if "clustering" == key and clustering:
                    for folder, dfs in value.items():
                        os.makedirs(os.path.join(root.clustering,folder), exist_ok=True)
                        vars[folder] = "clustering/{}".format(folder)
                        for name, df in dfs.items():
                            path = os.path.join(root.clustering,folder,name+".csv")
                            print1("Saving {}.csv at {}".format(name, path))
                            df.to_csv(path, header = True, index=False)

            try:
                from visualisation import generate_charts
                generate_charts()
            except:
                sprint("Charts not generated")


def create_dfs(references):
    root["dataframes"] = "dataframes"
    columns_raw = ["ID", "name", "chain"]
    for reference in references:
        columns_raw.extend(["rmsd_" + reference.name, "align_len_" + reference.name])
    vars["raw_monomers_df"] = pd.DataFrame(columns=columns_raw)

    colums_filtered = ["ID", "best_fit", "coverage (%query)", "coverage (%reference)", "rmsd", "identity (%)", "Rx",
                       "Ry", "Rz", "T"]
    vars["monomers_df"] = pd.DataFrame(columns=colums_filtered)

    vars["failed_df"] = pd.DataFrame(columns=["ID", "stage", "error", "details"])

    '''columns_alignment = ["ID"]
    for reference in references:
        columns_alignment.extend(["score_" + reference.name])
    vars["alignments_df"] = pd.DataFrame(columns=columns_alignment)'''

def create_clustering_dfs(references):
    vars["clustering"] = {}
    root["clustering"] = "dataframes/clustering"
    vars["clustering"]["contacts"] = {}
    root["contacts"] = "dataframes/clustering/contacts"
    vars["clustering"]["faces"] = {}
    root["faces"] = "dataframes/clustering/faces"
    vars["clustering"]["classified"] = {}
    root["classified"] = "dataframes/clustering/classified"
    for reference in references:
        print("Reference:", reference.name)
        ref_data = {"ResNum": [],
                    "ResName": []}
        for res in reference.structure.get_list():
            ref_data["ResNum"].append(res.id[1])
            ref_data["ResName"].append(res.resname)
        vars["clustering"]["contacts"][reference.name] = pd.DataFrame(data=ref_data)
        vars["clustering"]["faces"][reference.name] = pd.DataFrame(columns=["ID", "face1", "face2", "contact_face1", "contact_face2"])





def load_failed_dfs():
    sprint("Loading failed dataframes...")
    for file in os.listdir(root.dataframes):
        if "failed" in file:
            vars[file.split(".")[0]] = pd.read_csv(os.path.join(root.dataframes, file))










