import os

from utilities import *
from Globals import root, local, vars
import numpy as np
import pandas as pd




def save_dfs(force = [], general= True, clustering = True):
    excluded = ["classified"]
    sprint("Saving dataframes...")
    if "do_only" in vars:
        if len(vars.do_only) == 0 or vars.do_only is None:
            for key, value in vars.items():
                if "df" in key and general:
                    for e in excluded:
                        if e in key and e not in [f for f in force]:
                            continue
                    print1("Saving {}.csv".format(key))
                    value.to_csv(os.path.join(root.dataframes,f"{key}.csv"), header = True, index=False)
                if "clustering" == key and clustering:
                    for folder, dfs in value.items():
                        os.makedirs(os.path.join(root.clustering,folder), exist_ok=True)
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
    root["clustering"] = "dataframes/clustering"
    vars["clustering"] = {}
    vars["clustering"]["contacts"] = {}
    for reference in references:
        print("Reference:", reference.name)
        ref_data = {"ResNum": [],
                    "ResName": []}
        for res in reference.structure.get_list():
            ref_data["ResNum"].append(res.id[1])
            ref_data["ResName"].append(res.resname)
        vars["clustering"]["contacts"][reference.name] = pd.DataFrame(data=ref_data)
    vars["classified_df"] = pd.DataFrame(columns=["ID", "Best_Fit", "Best_Match", "Similarity", "Inverse"])



def load_failed_dfs():
    sprint("Loading failed dataframes...")
    for file in os.listdir(root.dataframes):
        if "failed" in file:
            vars[file.split(".")[0]] = pd.read_csv(os.path.join(root.dataframes, file))










