import os
from globals import root, local
from utilities import *
import pandas as pd
import matplotlib.pyplot as plt


def generate_charts():
    root["charts"] = "charts"
    if "super_filtered.csv" in os.listdir(root.dataframes):
        data = pd.read_csv(os.path.join(root.dataframes,"super_filtered.csv"))

        labels = data["best_fit"].unique()
        sizes = []
        for label in labels:
            sizes.append(len(data[data["best_fit"] == label]))

        fig, ax = plt.subplots()
        ax.pie(sizes, labels=labels)
        fig.savefig(os.path.join(root.charts,"super_filtered.png"))

    if "failed_df.csv" in os.listdir(root.dataframes):
        data = pd.read_csv(os.path.join(root.dataframes,"failed_df.csv"))

        labels = data["reason"].unique()
        sizes = []
        for label in labels:
            sizes.append(len(data[data["reason"] == label]))

        fig, ax = plt.subplots()
        ax.pie(sizes, labels=labels)
        fig.savefig(os.path.join(root.charts,"failed_df.png"))
