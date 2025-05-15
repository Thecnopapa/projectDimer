import os
import sys

from Globals import root, local, vars
from maths import angle_between_vectors
from pyMol import pymol_load_name, pymol_start, pymol_load_path, pymol_align_chains

from utilities import *
import pandas as pd
import matplotlib.pyplot as plt




def generate_piechart(df_name:str = None, column = None, extra_data:dict[str, int] = None, name= None):
    labels = []
    sizes = []
    large_labels = []
    if df_name is not None:
        assert column is not None
        if df_name in os.listdir(root.dataframes):
            data = pd.read_csv(os.path.join(root.dataframes, df_name))
            if len(data) > 0:
                labels = list(data[column].unique())
                for label in labels:
                    sizes.append(len(data[data[column] == label]))
    else:
        assert extra_data is not None
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
    if df_name is not None:
        ax.set_title("{} in {} (N = {})".format(column, df_name, total))
        fig_name = "{}.png".format(df_name.split(".")[0])
    else:
        ax.set_title("{} (N = {})".format(name, total))
        fig_name = "{}.png".format(name)
    fig.legend(title="Best Fit:", labels=labels, loc="lower right")

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

    # Classified piechart
    classified_chart()


def classified_chart(ref_name="GR"):
    root["charts"] = "charts"
    classified_df = pd.read_csv(os.path.join(root.classified, ref_name+".csv"))
    data = {}
    for n in range(11):
        subset = classified_df[classified_df["Similarity"] <= n*10 ]
        subset = subset[subset["Similarity"] > (n-1)*10]
        data[str((n-1)*10+1)+"-"+str(n*10)] =  len(subset)
    generate_piechart(extra_data=data, name = "Classified Similarities of GR")
    print(data)


def generate_html():
    pass


def show_objects(obj_list, args, mates = False, merged = False, paint_all_faces =True):

    merged = "merged" in sys.argv


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
                if type(item) is dict:
                    pca = item["pca"]
                else:
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
            from pyMol import pymol_start, pymol_load_path, pymol_symmetries, pymol_group, pymol_paint_contacts, pymol_format, pymol_draw_line, pymol_set_state
            pymol_start(show=True)
            for key, item in obj.__dict__.items():
                if type(item) == str:
                    #print(item.endswith(".pdb"), not "fractional" in item, not "merged" in item, not merged)
                    if item.endswith(".pdb"):
                        if "fractional" in item or ("merged" in item and not merged):
                            continue
                        if "many_pdbs" in item:
                            og = pymol_load_path(item, os.path.basename(item)+"_original")
                        elif "pdb_" in item:
                            pymol_load_path(item, os.path.basename(item) + "_processed")
                        else:
                            pymol_load_path(item, os.path.basename(item))
                        pymol_format("surface", os.path.basename(item), colour= "gray")
                        if "is_reference" in obj.__dict__.keys() or paint_all_faces:
                            if "best_fit" in obj.__dict__.keys():
                                if obj.best_fit == "GR":
                                    from pyMol import pymol_paint_all_faces
                                    pymol_paint_all_faces(obj)

                        elif "faces" in args and not paint_all_faces:
                            #print("Painting faces")
                            if "contacts_faces1" in obj.__dict__.keys():
                                pymol_paint_contacts(os.path.basename(item), obj.contacts_faces1[1:],
                                                     colour=obj.contacts_faces1[0])
                            if "contacts_faces2" in obj.__dict__.keys():
                                pymol_paint_contacts(os.path.basename(item), obj.contacts_faces2[1:],
                                                     colour=obj.contacts_faces2[0])
                            #print("Faces painted")
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
                            #print(obj.com)
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
                        if "face_coms" in obj.__dict__.keys():
                            from faces import GR_colours
                            from pyMol import pymol_sphere
                            for face, com in obj.face_coms.items():
                                pymol_sphere(com, colour=GR_colours[face], name="face_" + obj.id)
                if key == "monomer1" or key == "monomer2":
                    from faces import GR_colours
                    from pyMol import pymol_sphere
                    if "face_coms" in obj.__getattribute__(key).__dict__.keys():
                        for face, com in obj.__getattribute__(key).face_coms.items():
                            pymol_sphere(com, colour=GR_colours[face], name="face_"+obj.id)


                ### Development
                if key == "mate_paths" or key == "dimer_paths" and mates:
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

            from pyMol import pymol_format, pymol_orient, pymol_show_cell, pymol_hide, pymol_disable
            if paint_all_faces:
                pymol_format("spheres", "neighbour", "all", colour="rainbow", spectrum="b")
                pymol_format("mesh", "original", "all")
                pymol_format("mesh", "processed", "all")
            else:
                pymol_format("spheres", "neighbour", "all", colour="rainbow", spectrum="b")
                pymol_format("mesh", "original", "all", colour="white")
                pymol_format("mesh", "processed", "all", colour="blue")
            pymol_set_state(0)
            pymol_orient()
            pymol_show_cell()
            pymol_disable("pca")
            pymol_group(identifier="face", name="faces")
            try:
                pymol_hide("c", "label")
                pymol_disable("c")
            except:
                pass
            #pymol_group(identifier="rep", name="replaced")
            # pymol_group(identifier= "dimer")
            if mates:
                pymol_group(identifier="mate", name="mates")
            if merged:
                pymol_group(identifier="merged", name="merged")
            try:
                pymol_symmetries(og)
                pymol_group(identifier = "sym")
            except:
                pass
            pymol_set_state(2)



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
    print2("clusters-pca [pymol] [pca] [filter-threshold]")
    print2("clusters-pca-global [pymol] [pca] (filter-threshold !last)")



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



    elif any(arg in sys.argv[1] for arg in ["clusters-face", "clusters-pca", "clusters-pca-global"]):
        global_pca = False
        cluster_colname = "cluster"
        if "pca" in sys.argv[1]:
            pca = True
            if "global" in sys.argv[1]:
                global_pca = True
                cluster_colname = "global_cluster"
            else:
                global_pca = False

        else:
            pca = False
        tprint("Showing clustered data by interaction faces")

        filtered = None
        if len(sys.argv) >= 3:
            filtered = sys.argv[-1]
            print("Subset:", filtered)

        if pca:

            if global_pca:
                faces = [face for face in os.listdir(root.clustered_pcas) if ("clustered" not in face and "centres" not in face)]
            else:
                faces = os.listdir(root.clustered_pcas_GR)
            print("PCA clustering")
            print(faces)
        else:
            faces = os.listdir(root.cc_figs_GR)

        try:
            int(filtered)
            print("Filtering")
            for f in faces.copy():
                if not str(filtered) in f:
                    faces.remove(f)
        except:
            pass


        sprint("Clustered face-face combinations for GR:")
        for n, f in enumerate(faces):
            print1("{}: {}".format(n,f.split(".")[0]))
        if len(faces) > 1:
            while True:
                face = input("\n # Please select cluster to display (using the associated number):\n >> ")
                try:
                    face = faces[int(face)].split(".")[0]
                    break
                except:
                    pass
        elif len(faces) == 1:
            face = faces[0].split(".")[0]
        else:
            print(faces)
            print("No clustered dataframes found")
            print(faces[0])
        sprint("Selected interaction: {}".format(face))
        if pca:
            if global_pca:
                clustered_df = pd.read_csv(os.path.join(root.clustered_pcas, face + ".csv"), index_col=0)
            else:
                clustered_df = pd.read_csv(os.path.join(root.clustered_pcas_GR, face + ".csv"), index_col=0)
            #print(clustered_df)
            clustered_df.sort_values(by=[cluster_colname], inplace=True)

        else:
            clustered_df = pd.read_csv(os.path.join(root.clustered_GR, face+".csv"), index_col=0)
            clustered_df.sort_values(by=[cluster_colname, "similarity"], inplace=True)
        print(clustered_df.to_string(index=False))

        while True:
            c = input("\n # Please select cluster to display (int):\n >> ")
            try:
                if c == "all":
                    break
                c = [int(n) for n in c.split("/")]
                print(c)
                break
            except: pass
        if c != "all":
            print(clustered_df.query(" | ".join(["{} == {}".format(cluster_colname, n) for n in c])))
            subset = clustered_df.query(" | ".join(["{} == {}".format(cluster_colname, n) for n in c]))
        else:
            subset = clustered_df
        #subset.sort_values(by = "similarity", inplace = True)
        #print(subset.to_string(index=False))
        if not pca:
            threshold = input("\n # Please select minimum similarity threshold (int or all):\n >> ")
            if threshold == "":
                threshold = 0
            subset = subset[subset["similarity"] >= float(threshold)]
        print(subset.to_string(index=False))



        if "plot" in sys.argv or "mpl" in sys.argv:
            from clustering import plot_cc
            reference = load_references(identifier = "GR")[0]
            #print(reference)
            if pca:
                if global_pca:
                    in_path = os.path.join(root.clustered_pcas, face+".csv")
                else:
                    in_path = os.path.join(root.clustered_pcas_GR, face + ".csv")

            #print(in_path)
            #print(pd.read_csv(in_path))
            plot_cc(reference=reference, subset=c, subfolder=face, in_path=in_path, labels=True, pca=pca)




        if "pca" in sys.argv:
            print1("Plotting cluster PCAs")
            from faces import plot_pcas
            pcas = []
            progress = ProgressBar(len(subset))

            subset.sort_values("id", inplace=True)
            subset.reset_index(inplace=True)
            print(subset)
            for row in subset.itertuples():
                dimers = load_single_pdb(identifier=row.id, pickle_folder=local.dimers, quiet=True)
                for dimer in dimers:
                    dimer.pca.parent_id = dimer.id
                    pcas.append(dimer.pca)
                    progress.add(info=dimer.id)
            dimensions = [0,1,2]
            comps = [0,1,2]
            mode = "variance"
            cluster = None
            bandwidth=None
            for arg in sys.argv:
                if "d=" in arg:
                    dimensions = [int(i) for i in arg.split("=")[1]]
                if "c=" in arg:
                    comps = [int(i) for i in arg.split("=")[1]]
                if "mode=" in arg:
                    mode = arg.split("=")[1]
                if "cluster=" in arg:
                    cluster = int(arg.split("=")[1])
                if "bandwidth=" in arg:
                    bandwidth = float(arg.split("=")[1])
            print("dimension:", dimensions)
            subcluster_df = plot_pcas(pcas, title= "GR:({} : cluster {} / N = {})".format(face, c, len(pcas)),
                      dimensions=dimensions, mode=mode, bandwidth=bandwidth, comps=comps, cluster=cluster)

            print(subcluster_df)
            if subcluster_df is not None:
                subset["subcluster"] = subcluster_df["cluster"]
                cluster_colname = "subcluster"
                print(subset)
                while True:
                    c = input("\n # Please select cluster to display (int):\n >> ")
                    try:
                        if c == "all":
                            break
                        c = [int(n) for n in c.split("/")]
                        break
                    except:
                        pass
                if c != "all":
                    print(clustered_df.query(" | ".join(["{} == {}".format(cluster_colname, n) for n in c])))
                    subset = clustered_df.query(" | ".join(["{} == {}".format(cluster_colname, n) for n in c]))
                else:
                    subset = subset


        if "gesamt" in sys.argv:
            if c == "all":
                clusters_to_display = sorted(list(subset[cluster_colname].unique()))
            else:
                clusters_to_display = [c]
            sprint("Displaying clusters:",cluster_colname, clusters_to_display)
            for n in clusters_to_display:
                chains_to_align = []
                first_to_align = None
                print("CLUSTER", n)
                print(subset)
                pymol_subset = subset[subset[cluster_colname] == n].sort_values(by="id", ascending=True)
                print(pymol_subset)
                for row in pymol_subset.itertuples():
                    dimers = load_single_pdb(identifier=row.id, pickle_folder=local.dimers, quiet=True)
                    for dimer in dimers:
                        chains_to_align.append(dimer.replaced_path)

                from superpose import superpose_multiple
                data = superpose_multiple(chains_to_align)
                print(data)
                quit()



        elif "pymol" in sys.argv:
            from pyMol import  *
            pymol_start(show=False)
            pymol_set_state(2)
            cluster_paths = []
            if c == "all":
                clusters_to_display = sorted(list(subset[cluster_colname].unique()))
            else:
                clusters_to_display = [*c]
            sprint("Displaying clusters:",cluster_colname, clusters_to_display)
            for n in clusters_to_display:
                chains_to_align = []
                first_to_align = None
                print("CLUSTER", n)
                print(subset)
                pymol_subset = subset[subset[cluster_colname] == n].sort_values(by="id", ascending=True)

                print(pymol_subset)
                for row in pymol_subset.itertuples():
                    dimers = load_single_pdb(identifier=row.id, pickle_folder=local.dimers, quiet=True)
                    for dimer in dimers:

                        pymol_load_path(dimer.merged_path, dimer.id)
                        '''if first_to_align is None:
                            first_to_align = dimer.face1
                        if dimer.face1 == first_to_align:
                                chains_to_align.append((row.id, row.id.split("_")[-3][-2]))
                        else:
                            chains_to_align.append((row.id, row.id.split("_")[-3][-1]))'''

                        chains_to_align.append((row.id, *row.id.split("_")[-3][-2:]))

                        pymol_colour("gray", dimer.id)

                        # pymol_paint_contacts(os.path.basename(dimer.id), dimer.contacts_faces1[1:], colour=dimer.contacts_faces1[0])
                        # pymol_paint_contacts(os.path.basename(dimer.id), dimer.contacts_faces2[1:], colour=dimer.contacts_faces2[0])

                        pymol_paint_all_faces(dimer)

                        '''from faces import pca_to_lines
                        for pca in [dimer.pca1, dimer.pca2]:
                            point_list = pca_to_lines(pca["pca"], com=pca["com"], just_points=True)
                            #print("point list:", point_list)
                            for p in point_list:
                                pymol_draw_line(coord1=p[0], coord2=p[1], name="pca", quiet=False)'''

                print(chains_to_align)
                # pymol_align_chains(chains_to_align)
                pymol_align_all([obj[0] for obj in chains_to_align])

                pymol_align_chains_best(chains_to_align, double_best=True, cluster=n)

                sele = "({})".format(" or ".join(chain[0] for chain in chains_to_align))
                pymol_move(sele=sele, distance=[150 * n, 0, 0])
                pymol_group([chain[0] for chain in chains_to_align], name=str(n))

                print(pymol_subset.to_string())
                cluster_path = pymol_save_cluster(chains_to_align, name="CLUSTER_{}.pdb".format(n) )
                cluster_paths.append([cluster_path, chains_to_align])


            pymol_set_state(2)
            pymol_orient()
            print(subset.to_string())




            if "post-pca" in sys.argv:
                from Bio.PDB import PDBParser
                from faces import *
                print(cluster_paths)



                for data in cluster_paths:
                    print(data)
                    cluster_data = {"names": [],
                                    "models": [],
                                    "pcas": [],
                                    "path": "",
                                    "corners": [],
                                    "coms": []}
                    cluster_objects = [chain_info[0] for chain_info in data[1]]
                    cluster_data["names"] = cluster_objects
                    path = data[0]
                    cluster_data["path"] = path
                    cluster_structure = PDBParser(QUIET=True).get_structure(os.path.basename(path.split(".")[0]), path)
                    #print(cluster_structure.__dict__)
                    [print(model) for model in cluster_structure.get_models()]
                    models = [model for model in cluster_structure.get_models() if model.id % 2 == 0]
                    [print(model, name) for model, name in zip(models, cluster_objects)]
                    spheres = []
                    for model, name in zip(models, cluster_objects):
                        print(model, list(model.get_chains()))
                        com = find_com(model.get_atoms())
                        pca = get_pca(model, com=com)
                        cluster_data["models"].append(model)
                        cluster_data["pcas"].append(pca)
                        cluster_data["coms"].append(com)

                        components = pca.components_
                        variances = pca.explained_variance_
                        ratios = pca.explained_variance_ratio_
                        if pca.inverse:
                            components[1], components[2] = components[2].copy(), components[1].copy()
                            variances[1], variances[2] = variances[2].copy(), variances[1].copy()
                        sphere_coords = com
                        corner = [0,0,0]
                        for n, (component, variance, ratio) in enumerate(zip(components, variances, ratios)):
                            print("sphere-coords:", sphere_coords)
                            print(com, component, variance)
                            pymol_draw_line(com, tuple([c+(co*variance) for c, co in zip(com, component)]), name="{}_component_{}".format(name, n), quiet=False)
                            sphere_coords = add(sphere_coords, tuple([co*variance for co in component]))
                            corner = add(corner, tuple([co*ratio for co in component]))
                        print("sphere-coords:",sphere_coords)
                        spheres.append(sphere_coords)
                        cluster_data["corners"].append(corner)
                        pymol_sphere(sphere_coords, name=name+"_corner")
                    pymol_group("component", name="components")
                    print(cluster_data)
                    df = pd.DataFrame(cluster_data)
                    print(df)

                    df["id"] = cluster_data["names"]
                    df[["_0","_1","_2"]] = cluster_data["corners"]
                    if "sm" in sys.argv:
                        sm = pd.DataFrame(columns=["id1", "id2", "angle"])
                        done= []
                        for row1 in df.itertuples():
                            for row2 in df.itertuples():
                                if row2.names in done or row1.names == row2.names:
                                    continue
                                angle = angle_between_vectors(vector(row1.coms, row1.corners), vector(row2.coms, row2.corners))
                                print(angle)
                                sm.loc[len(sm)] = [row1.names, row2.names, angle]
                            done.append(row1.names)
                        print(sm)





                    else:
                        #sele = "({})".format(" or ".join(chain[0] for chain in cluster_data["names"]))
                        #pymol_move(sele=sele, distance=[150 * n, 0, 0])
                        #points=[add_multiple(pca.components_[1,2] for pca in cluster_data["pcas"]]
                        #print(points)

                        print(df)
                        df["cluster"] = quick_cluster(df[["_0","_1", "_2"]], bandwidth=0.08)
                        df.sort_values(by=["cluster"], ascending=True, inplace=True)
                        print(df.to_string())


                    for row in df.itertuples():
                        pymol_group(row.id, path.split(".")[0][-1]+"_sc_"+str(row.cluster))

                    session_path = pymol_save_temp_session(name=os.path.basename(path).split(".")[0]+".pse")
                    print(df)
                    open_session_terminal(session_path)
                    plot_points(df)




            else:
                session_path = pymol_save_temp_session()
                open_session_terminal(session_path)
                print("Session temporarily saved at:")
                print(session_path)
                pymol_close()







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







    elif "clusters2" in sys.argv[1]:
        sprint("Showing clusters v2")
        for n, file in enumerate(os.listdir(root.dihedral_clusters)):
            print1(n, ":", file)

        i = int_input("Select df to display:\n")
        file = os.listdir(root.dihedral_clusters)[i]
        dihedrals_path = os.path.join(root.dihedral_clusters, file)
        df = pd.read_csv(dihedrals_path, index_col=0)
        print(df)
        c = int_input("Select cluster to display {}:\n".format([int(a) for a in set(df["angle_cluster"].values)]))
        df = df[df["angle_cluster"] == c]
        print(df.to_string())

        if "pymol" in sys.argv:
            from pymol import *
            pymol_start(show=True)
            chains_to_align = []
            chains_to_align_reverse = []
            i = 0
            for row in df.itertuples():
                if i ==10:
                    break
                dimer = load_single_pdb(identifier=row.id, pickle_folder=local.dimers)[0]
                name = pymol_load_path(dimer.replaced_path, row.id + str(row.is1to2))
                if row.is1to2:
                    chains_to_align.append([name, row.mon1])
                    if i == 0:
                        chains_to_align_reverse.append([name, row.mon1])
                else:
                    chains_to_align_reverse.append([name, row.mon2])
                    if i == 0:
                        chains_to_align.append([name, row.mon1])
                i+=1
            pymol_align_chains(chains_to_align)
            pymol_align_chains(chains_to_align_reverse)

        elif "plot" in sys.argv:
            from clustering import plot_dihedrals
            plot_dihedrals(dihedrals_path, subset_col="angle_cluster", subset=c, save=False, label_col="id", only_first=10)





            pass






    else:
        print_available_commands()


    eprint("Done visualising")




























