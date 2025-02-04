import os
from globals import root, local, vars
from utilities import *
import pandas as pd
import matplotlib.pyplot as plt


def generate_charts():
    root["charts"] = "charts"
    sprint("Generating charts")
    if "monomers_df.csv" in os.listdir(root.dataframes):
        data = pd.read_csv(os.path.join(root.dataframes,"monomers_df.csv"))
        if len(data) > 0:

            labels = list(data["best_fit"].unique())
            sizes = []
            large_labels = []
            for label in labels:
                sizes.append(len(data[data["best_fit"] == label]))
            print(sizes)
            labels += ["Missing"]
            failed_df = pd.read_csv(os.path.join(root.dataframes,"failed_df.csv"))
            sizes.append(len(failed_df[failed_df["stage"] == "monomer"]))

            print(sizes)
            total = sum(sizes)
            for index, (label, size) in enumerate(zip(labels, sizes)):
                labels[index] = "{} ({}%)".format(label, round(size/total*100))
                if round(size/total*100) >= 3:
                    large_labels.append(label)
                else:
                    large_labels.append("")
            print(labels)
            print(sizes)
            fig, ax = plt.subplots()
            ax.pie(sizes, labels=large_labels, startangle=90)
            ax.set_title("Best fit reference to each monomer (N = {})".format(total))
            fig.legend(title = "Best Fit:", labels = labels, loc = "lower right")
            fig.savefig(os.path.join(root.charts,"monomers_df.png"))
            print1("monomers_df.png generated")






    if "failed_df.csv" in os.listdir(root.dataframes):
        data = pd.read_csv(os.path.join(root.dataframes,"failed_df.csv"))
        if len(data) > 0:
            labels = data["error"].unique()
            sizes = []
            for label in labels:
                sizes.append(len(data[data["error"] == label]))

            fig, ax = plt.subplots()
            ax.pie(sizes, labels=labels)
            fig.savefig(os.path.join(root.charts,"failed_df.png"))

if __name__ == "__main__":
    import globals

    globals.set_root("../")
    if os.name == "nt":
        globals.set_local("C:/Users/iainv/localdata/_local/projectB")
    elif os.name == "posix":
        globals.set_local("/localdata/iain/_local/projectB")
    from globals import root, local, vars
    generate_charts()
