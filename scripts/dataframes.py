import os

from utilities import *
from globals import root, local, vars
import numpy as np
import pandas as pd

from visualisation import generate_charts


def save_dfs():
    sprint("Saving dataframes...")
    for key, value in vars.items():
        if "df" in key:
            print1("Saving {}".format(key))
            value.to_csv(os.path.join(root.dataframes,f"{key}.csv"), header = True, index=False)
    generate_charts()