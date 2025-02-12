import os

from utilities import *
from globals import root, local, vars
import numpy as np
import pandas as pd

from visualisation import generate_charts


def save_dfs():
    sprint("Saving dataframes...")
    if "do_only" in vars:
        if len(vars.do_only) == 0 or vars.do_only is None:
            for key, value in vars.items():
                if "df" in key:
                    print1("Saving {}.csv".format(key))
                    value.to_csv(os.path.join(root.dataframes,f"{key}.csv"), header = True, index=False)
            generate_charts()


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

    columns_alignment = ["ID"]
    for reference in references:
        columns_alignment.extend(["score_" + reference.name])
    vars["alignments_df"] = pd.DataFrame(columns=columns_alignment)