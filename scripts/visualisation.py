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


def show_objects(obj_list, args):
    for obj in obj_list:
        sprint(obj.id)
        print(obj.__dict__)
        for key, item in obj.__dict__.items():
            #print(str(type(item)))
            if type(item) == list:
                print1(key,":")
                for i in item:
                    print2(i)
            if type(item) in (str, int, float,  bool):
                print1(key,":",item)
            if type(item) == dict:
                print1(key, ":")
                for k, v in item.items():
                    print2(k,":",v)
            if "molecules" in str(type(item)):
                print1(key, ":", item)
            if "Bio" in str(type(item)):
                print1(key, ":", item)
            if item is None:
                print1(key, ":", "None")
        if "pymol" in args:
            from pyMol import pymol_start, pymol_load_path, pymol_format, pymol_set_state
            pymol_start(show=True)
            for key, item in obj.__dict__.items():
                if type(item) == str:
                    if item.endswith(".pdb") and not "fractional" in item:
                        if "many_pdbs" in item:
                            pymol_load_path(item, os.path.basename(item)+"_original")
                        elif "pdb_" in item:
                            pymol_load_path(item, os.path.basename(item) + "_processed")
                        else:
                            pymol_load_path(item)
                        pymol_format("surface", "neighbour", "all", colour="rainbow", spectrum="b")
                        pymol_format("mesh", "original", "all", colour="white")
                        pymol_format("mesh", "processed", "all", colour="white")
                        pymol_set_state(1)


def print_available_commands():
    sprint("Available commands:")
    print1("To show saved data, add \"pymol\" to load structures in pymol")
    print2("molecule + ID e.g 1M2Z")
    print2("monomer + ID e.g 1M2Z_A")
    print2("dimer + ID e.g 1M2Z_AD")
    print1("To show clustering data:")
    print2("clusters-eva")
    print2("clusters-cc + [reference_name] e.g. GR (includes all data from clusters-eva)")
    print2("clusters-score (only for GR)")



if __name__ == "__main__":
    tprint("Visualising data")


    if len(sys.argv) < 2:
        sprint("Command not provided")
        print_available_commands()
        eprint("Done visualising")
        quit()



    import setup
    from Globals import root, local, vars
    from imports import *

    import sys

    sprint("Argvs:")
    print(sys.argv)



    if "dimer" in sys.argv[1] and len(sys.argv[2:]) != 0:
        vars["do_only"] = sys.argv[2:]
        dimers = load_dimers()
        tprint("Showing dimers")
        show_objects(dimers, sys.argv[2:])

    elif "monomer" in sys.argv[1] and len(sys.argv[2:]) != 0:
        vars["do_only"] = sys.argv[2:]
        monomers = load_monomers()
        tprint("Showing monomers")
        show_objects(monomers, sys.argv[2:])

    elif ("molecule" in sys.argv[1] or "pdb" in sys.argv[1]) and len(sys.argv[2:]) != 0:
        vars["do_only"] = sys.argv[2:]
        molecules = load_from_files(local.many_pdbs)
        tprint("Showing molecules")
        show_objects(molecules, sys.argv[2:])


    elif ("clusters-eva" in sys.argv[1]):
        tprint("Showing clusters vs eva data")
        classified_df = pd.read_csv(os.path.join(root.dataframes, "classified_df.csv"),index_col=0)
        if len(sys.argv[2:]) == 0:
            print(classified_df.sort_values("Best_Match").to_string())


        elif len(sys.argv[2:]) == 1:
            #print(type(sys.argv[2]))
            filtered_df =classified_df[classified_df["Best_Match"] == int(sys.argv[2])]
            print(filtered_df.to_string())

            selection = input("Select similarity threshold to display (\"all\" or int), not inclusive:\n >>")

            if len(selection) > 0:


                if selection == "all":
                    threshold = -1
                else:
                    threshold = int(selection)

                filtered_df = filtered_df[filtered_df["Similarity"] > threshold]
                if len(filtered_df) > 0:
                    from pyMol import *
                    pymol_start(show=True)
                    chains_to_align = []
                    for row in filtered_df.itertuples():
                        print1(row.ID)
                        file_path = os.path.join(local.dimers_merged, row.ID+"_merged.pdb")
                        pymol_load_path(file_path)
                        if row.Inverse:
                            chains_to_align.append(row.ID[-1])
                        else:
                            chains_to_align.append(row.ID[-2])
                    pymol_set_state(2)
                    pymol_align_chains(chains_to_align)
                    pymol_align_all()
                    #pymol_colour_all("chainbows")

                else:
                    print1("No matches found")

    elif ("clusters-cc" in sys.argv[1]):
        tprint("Showing clusters CC + Kmeans")
        print(sys.argv)
        if len(sys.argv[2:]) >0:
            reference_name = sys.argv[2]
        else:
            reference_name = "GR"


        sprint("Showing results for: {}".format(reference_name))

        clustered_df = pd.read_csv(os.path.join(root.dataframes, "{}_cc_clustered.csv".format(reference_name)), index_col=0).sort_values("cluster")

        print(clustered_df.to_string())
        selection = input("Select cluster to display [cc/eva] + int:\n >>")

        if type(selection) is list:
            if "eva" in selection[0]:
                filtered_df = clustered_df[clustered_df["group"] == int(selection[1])]
            elif "cc" in selection[0]:
                filtered_df = clustered_df[clustered_df["cluster"] == int(selection[1])]
            else:
                filtered_df = clustered_df[clustered_df["cluster"] == int(selection[1])]
        elif len(selection) > 0:
            filtered_df = clustered_df[clustered_df["cluster"] == int(selection)]
        else:
            quit()
        print(selection)
        print(filtered_df.to_string())


        if len(filtered_df) > 0:
            from pyMol import *

            pymol_start(show=True)
            for row in filtered_df.itertuples():
                print1(row.id)
                file_path = os.path.join(local.dimers_merged, row.id + "_merged.pdb")
                pymol_load_path(file_path)
            pymol_set_state(2)
            pymol_align_all()
            # pymol_colour_all("chainbows")

        else:
            print1("No matches found")

    elif "cluster-score" in sys.argv[1] or "clusters-score" in sys.argv[1]:

        from clustering import calculate_scores_GR
        calculate_scores_GR(pd.read_csv(os.path.join(root.dataframes, "GR_cc_clustered.csv"), index_col=0).sort_values("cluster"), save=False)


    else:
        print_available_commands()


    eprint("Done visualising")


























