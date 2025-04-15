import os
import sys

from Globals import root, local, vars
from pyMol import pymol_start, pymol_load_name, pymol_load_path, pymol_set_state, pymol_align_chains, pymol_align_all, \
    pymol_paint_contacts, pymol_colour, pymol_draw_line, pymol_move, pymol_orient, pymol_group
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
        #print(obj.__dict__)
        for key, item in obj.__dict__.items():
            if key in ["lines", "c_lines", "sasas1D", "sasas2D", "full_array","contacts_faces1", "contacts_faces2",
                       "contacts", "contacts_symm", "contacts_sasa" ]:
                try:
                    print1(key, ": OMITTED (len: {})".format(len(item)))
                except:
                    print1(key, ": OMITTED (no len)")
                continue
            #print(str(type(item)))
            elif key == "pca":
                print1("pca:")
                pca = item
                for n, component in enumerate(pca.components_):
                    print2("Component {}: {} / Value: {} --> Vector: {}".format(n,
                                                                                pca.components_[n],
                                                                                pca.explained_variance_[n],
                                                                                pca.components_[n] *
                                                                                pca.explained_variance_[n]))

            elif type(item) in (list, tuple, set):
                print1(key,":")
                for i in item:
                    print2(i)
            elif type(item) in (str, int, float,  bool):
                print1(key,":",item)
            elif type(item) == dict:
                print1(key, ":")
                for k, v in item.items():
                    print2(k,":",v)
            elif "molecules" in str(type(item)):
                print1(key, ":", item, id(item))
            elif "Bio" in str(type(item)):
                print1(key, ":", item, id(item))
            elif "pandas" in str(type(item)):
                print1(key, ":\n", item.iloc[:, 0:2])#list(item.columns.to_series()))
            elif item is None:
                print1(key, ":", "None")
            else:
                print1(key, ":", item)
        if "pymol" in args:
            from pyMol import pymol_start, pymol_load_path, pymol_symmetries, pymol_group, pymol_paint_contacts, pymol_format, pymol_draw_line
            pymol_start(show=True)
            for key, item in obj.__dict__.items():
                if type(item) == str:
                    if item.endswith(".pdb") and not "fractional" in item:
                        if "many_pdbs" in item:
                            og = pymol_load_path(item, os.path.basename(item)+"_original")
                        elif "pdb_" in item:
                            pymol_load_path(item, os.path.basename(item) + "_processed")
                        else:
                            pymol_load_path(item)
                        pymol_format("surface", os.path.basename(item), colour= "gray")
                        if "faces" in args or True:
                            print("Painting faces")
                            if "contacts_faces1" in obj.__dict__.keys():
                                pymol_paint_contacts(os.path.basename(item), obj.contacts_faces1[1:],
                                                     colour=obj.contacts_faces1[0])
                            if "contacts_faces2" in obj.__dict__.keys():
                                pymol_paint_contacts(os.path.basename(item), obj.contacts_faces2[1:],
                                                     colour=obj.contacts_faces2[0])
                            print("Faces painted")
                        if "contacts_sasa" in obj.__dict__.keys():
                            pymol_paint_contacts(os.path.basename(item), obj.contacts_sasa, colour ="red")
                            pass
                        if "contacts_symm" in obj.__dict__.keys():
                            pymol_paint_contacts(os.path.basename(item), obj.contacts_symm)
                            pass
                        if "pca" in obj.__dict__.keys():
                            from faces import pca_to_lines
                            if type(obj.pca) is dict:
                                pca = obj.pca["pca"]
                            else:
                                pca = obj.pca
                            print(obj.com)
                            point_list = pca_to_lines(pca, com=obj.com, just_points=True)
                            for p in point_list:
                                pymol_draw_line(coord1=p[0], coord2=p[1], name="pca")

                        elif "pca1" in obj.__dict__.keys():
                            from faces import pca_to_lines
                            for pca in [obj.pca1, obj.pca2]:
                                point_list = pca_to_lines(pca["pca"], com=pca["com"], just_points=True)
                                print("point list:", point_list)
                                for p in point_list:

                                    pymol_draw_line(coord1=p[0], coord2=p[1], name="pca", quiet=False)


                ### Development
                if key == "mate_paths" or key == "dimer_paths":
                    for mate in item:
                        #print(mate)
                        pymol_load_path(mate, os.path.basename(mate))

                if key == "lines" and False:
                    from pyMol import pymol_draw_line
                    l = 2
                    for line_set in item:
                        for line in line_set:
                            #print(l,line)
                            pymol_draw_line(line[0], line[1], state=l)
                        l += 1
                if key == "c_lines":
                    from pyMol import pymol_draw_line
                    for line in item:
                        pymol_draw_line(line[0], line[1], name  = "c")
                if key == "contacs_sasa":
                    pass

                ###

            from pyMol import pymol_format, pymol_set_state, pymol_orient, pymol_show_cell, pymol_hide
            pymol_format("spheres", "neighbour", "all", colour="rainbow", spectrum="b")
            pymol_format("mesh", "original", "all", colour="white")
            pymol_format("mesh", "processed", "all", colour="blue")
            pymol_set_state(1)
            pymol_orient()
            pymol_show_cell()
            try:
                pymol_hide("c", "label")
            except:
                pass
            pymol_group(identifier="mate", name="mates")
            #pymol_group(identifier= "dimer")
            pymol_group(identifier="rep", name="replaced")
            pymol_group(identifier="merged", name="merged")
            try:
                pymol_symmetries(og)
                pymol_group(identifier = "sym")
            except:
                pass



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
    print2("clusters-faces (only for GR)")



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

    elif "ref" in sys.argv[1] and len(sys.argv[2:]) != 0:
        refs = load_references(identifier = sys.argv[2])
        tprint("Showing references")
        show_objects(refs, sys.argv[2:])



    elif "clusters-face" in sys.argv[1] or "clusters-pca" in sys.argv[1]:
        if "clusters-pca" in sys.argv[1]:
            pca = True
        else:
            pca = False
        tprint("Showing clustered data by interaction faces")
        if pca:
            faces = os.listdir(root.pca_figs_GR)
        else:
            faces = os.listdir(root.cc_figs_GR)
        sprint("Clustered face-face combinations for GR:")
        for n, f in enumerate(faces):
            print1("{}: {}".format(n,f.split(".")[0]))
        while True:
            face = input("\n # Please select cluster to display (using the associated number):\n >> ")
            try:
                face = faces[int(face)].split(".")[0]
                break
            except:
                pass
        sprint("Selected interaction: {}".format(face))
        if pca:
            clustered_df = pd.read_csv(os.path.join(root.clustered_pcas_GR, face + ".csv"), index_col=0)
            clustered_df.sort_values(by=["cluster"], inplace=True)
        else:
            clustered_df = pd.read_csv(os.path.join(root.clustered_GR, face+".csv"), index_col=0)
            clustered_df.sort_values(by=["cluster", "similarity"], inplace=True)
        print(clustered_df.to_string(index=False))

        while True:
            c = input("\n # Please select cluster to display (int):\n >> ")
            try:
                if c == "all":
                    break
                c = int(c)
                break
            except: pass
        if c != "all":
            subset = clustered_df[clustered_df["cluster"] == int(c)]
        else:
            subset = clustered_df
        #subset.sort_values(by = "similarity", inplace = True)
        print(subset.to_string(index=False))
        if not pca:
            threshold = input("\n # Please select minimum similarity threshold (int):\n >> ")
            if threshold == "":
                threshold = 0
            subset = subset[subset["similarity"] >= float(threshold)]
        print(subset.to_string(index=False))




        if "pymol" in sys.argv:
            from pymol import *
            pymol_start(show=True)

            if c == "all":
                l = list(subset["cluster"].unique())
            else:
                l = [c]
            for n in range(len(l)):
                chains_to_align = []
                first_to_align = None
                print("CLUSTER", n, l[n])
                if c == "all":
                    subset = clustered_df[clustered_df["cluster"] == l[n]]
                print(subset)
                for row in subset.itertuples():
                    dimers = load_single_pdb(identifier=row.id, pickle_folder=local.dimers)
                    for dimer in dimers:
                        # TODO: delete after full run
                        #dimer.get_contacts(force=True)
                        #dimer.get_faces()
                        #############################
                        pymol_load_path(dimer.merged_path, dimer.id)
                        if first_to_align is None:
                            first_to_align = dimer.face1
                        if dimer.face1 == first_to_align:
                                chains_to_align.append((row.id, row.id.split("_")[-3][-2]))
                        else:
                            chains_to_align.append((row.id, row.id.split("_")[-3][-1]))


                        pymol_colour("gray", dimer.id)
                        pymol_paint_contacts(os.path.basename(dimer.id), dimer.contacts_faces1[1:],
                                                 colour=dimer.contacts_faces1[0])
                        pymol_paint_contacts(os.path.basename(dimer.id), dimer.contacts_faces2[1:],
                                                 colour=dimer.contacts_faces2[0])



                        '''from faces import pca_to_lines
                        for pca in [dimer.pca1, dimer.pca2]:
                            point_list = pca_to_lines(pca["pca"], com=pca["com"], just_points=True)
                            #print("point list:", point_list)
                            for p in point_list:
                                pymol_draw_line(coord1=p[0], coord2=p[1], name="pca", quiet=False)'''


                print(chains_to_align)
                pymol_align_chains(chains_to_align)
                sele = "({})".format(" or ".join(chain[0] for chain in chains_to_align))
                print(sele)
                pymol_move(sele=sele, distance=[150*n, 0, 0])
                pymol_group([chain[0] for chain in chains_to_align], name=str(n))
            pymol_set_state(2)
            pymol_orient()


        if "plot" in sys.argv or "mpl" in sys.argv:
            from clustering import plot_cc
            reference = load_references(identifier = "GR")[0]
            #print(reference)
            in_path = os.path.join(root.clustered_GR, face+".csv")
            #print(in_path)
            #print(pd.read_csv(in_path))
            plot_cc(reference=reference, subset=c, subfolder=face, in_path=in_path, labels=True)

        if "pca" in sys.argv:
            from faces import plot_pcas
            pcas = []
            for row in subset.itertuples():
                dimers = load_single_pdb(identifier=row.id, pickle_folder=local.dimers)
                for dimer in dimers:
                    pcas.append(dimer.pca)
            plot_pcas(pcas, title= "GR:({} : cluster {} / N = {})".format(face, c, len(pcas)))







    elif ("clusters-eva" in sys.argv[1]):
        tprint("Showing clusters vs eva data")
        classified_df = pd.read_csv(os.path.join(root.classified, "GR.csv")).sort_values("Similarity")#,index_col=0)
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

                filtered_df = filtered_df[filtered_df["Similarity"] > threshold].sort_values("Similarity")
                if len(filtered_df) > 0:
                    from pyMol import *
                    pymol_start(show=True)
                    chains_to_align = []
                    for row in filtered_df.itertuples():
                        print1(row.ID)
                        file_path = os.path.join(local.dimers_merged, row.ID+"_dimer_merged.pdb")
                        pymol_load_path(file_path)
                        if row.Inverse:
                            chains_to_align.append((row.ID,row.ID.split("_")[-3][-1]))
                        else:
                            chains_to_align.append((row.ID,row.ID.split("_")[-3][-2]))
                    pymol_set_state(2)
                    pymol_align_chains(chains_to_align)
                    #pymol_align_all()
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

        from clustering import calculate_scores_GR, plot_cc

        calculate_scores_GR(pd.read_csv(os.path.join(root.dataframes, "GR_cc_clustered.csv"), index_col=0).sort_values("cluster"), save=False)


    else:
        print_available_commands()


    eprint("Done visualising")


























