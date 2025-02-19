import os
import sys

from Globals import root, local, vars
from utilities import *
import pandas as pd
import matplotlib.pyplot as plt




def generate_piechart(df_name:str, column, extra_data:dict[str, int] = None):
    if df_name in os.listdir(root.dataframes):
        data = pd.read_csv(os.path.join(root.dataframes, df_name))
        if len(data) > 0:
            labels = list(data[column].unique())
            sizes = []
            large_labels = []
            for label in labels:
                sizes.append(len(data[data[column] == label]))
            if extra_data is not None:
                for key, value in extra_data.items():
                    labels += [key]
                    sizes += [value]

            total = sum(sizes)
            for index, (label, size) in enumerate(zip(labels, sizes)):
                labels[index] = "{} ({}%)".format(label, round(size / total * 100))
                if round(size / total * 100) >= 3:
                    large_labels.append(label)
                else:
                    large_labels.append("")
            fig, ax = plt.subplots()
            ax.pie(sizes, labels=large_labels, startangle=90)
            ax.set_title("{} in {} (N = {})".format(column, df_name, total))
            fig.legend(title="Best Fit:", labels=labels, loc="lower right")
            fig_name = "{}.png".format(df_name.split(".")[0])
            fig_path = os.path.join(root.charts, fig_name)
            fig.savefig(fig_path)
            print1("{} generated".format(fig_name))





def generate_charts():
    root["charts"] = "charts"
    sprint("Generating charts")


    # Successful monomers piechart:
    failed_df = pd.read_csv(os.path.join(root.dataframes, "failed_df.csv"))
    extra_data ={"Missing": (len(failed_df[failed_df["stage"] == "monomer"]))}
    generate_piechart("monomers_df.csv", "best_fit", extra_data)

    # Errors piechart:
    generate_piechart("failed_df.csv", "error")

def generate_html():
    pass


if __name__ == "__main__":
    tprint("Visualising data")
    import setup
    from Globals import root, local, vars
    from imports import *

    import sys

    sprint("Argvs:")
    print(sys.argv)


    if ("dimer" in sys.argv[1] or "dimers" in sys.argv[1]) and len(sys.argv[2:]) != 0:
        vars["do_only"] = sys.argv[2:]
        dimers = load_dimers()
        tprint("Showing dimers")
        for dimer in dimers:
            sprint(dimer.id)
            for key, item in dimer.__dict__.items():
                if type(item) == list :
                    print1(key,":")
                    for i in item:
                        print2(i)
                if type(item) == (str or int or float or bool):
                    print1(key,":",item)
                if type(item) == dict:
                    print1(key, ":")
                    for k, v in item.items():
                        print2(k,":",v)
                if "molecules" in str(type(item)):
                    print1(key, ":", item)




    elif ("cluster" in sys.argv[1] or "clusters" in sys.argv[1]):
        tprint("Showing clusters")
        classified_df = pd.read_csv(os.path.join(root.dataframes, "classified_df.csv"),index_col=0)
        if len(sys.argv[2:]) == 0:
            print(classified_df.sort_values("Best_Match").to_string())


        elif len(sys.argv[2:]) == 1:
            print(type(sys.argv[2]))
            filtered_df =classified_df[classified_df["Best_Match"] == int(sys.argv[2])]
            print(filtered_df.to_string())

            selection = input("Select similarity threshold to display (\"all\" or int), not inclusive:\n >>")
            from pyMol import *

            pymol_start(show=True)
            if len(selection) > 0:
                if selection == "all":
                    threshold = -1
                else:
                    threshold = int(selection)

                filtered_df = filtered_df[filtered_df["Similarity"] > threshold]
                for name in filtered_df["ID"].values:
                    print(name)
                    file_path = os.path.join(local.dimers_merged, name+"_merged.pdb")
                    pymol_load_path(file_path)
                pymol_align_all()





    eprint("Done visualising")


























