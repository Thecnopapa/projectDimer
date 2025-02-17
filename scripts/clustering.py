import os
from utilities import *
from Globals import root, local, vars
import numpy as np
import pandas as pd





def generate_sm(reference, force=False):
    sprint("Generating SM for {}".format(reference.name))
    sasa_name = '{}_sasas.csv'.format(reference.name)
    if '{}_sm_ssd.csv'.format(reference.name) in os.listdir(root.dataframes) and not force:
        print1("Skipping SM generation for {}".format(reference.name))
        return
    if not sasa_name in os.listdir(root.dataframes):
        print1("Dataframe not found: {}".format(os.path.join(root.dataframes, sasa_name)))
        return
    sasas_df = pd.read_csv(os.path.join(root.dataframes, sasa_name), index_col=0)
    print(sasas_df)

    sm_ssd = pd.DataFrame(columns=["dimer1", "dimer2", "index1", "index2", "similarity"])
    n_dimers = len(sasas_df.columns) - 2
    n_res = len(sasas_df)
    progress = ProgressBar(n_dimers * n_dimers)
    index1 = 0
    for id1, contacts1 in zip(sasas_df.columns, sasas_df._iter_column_arrays()):
        if id1 in ["ResNum", "ResName"]:
            continue
        index1 += 1
        index2 = 0
        for id2, contacts2 in zip(sasas_df.columns, sasas_df._iter_column_arrays()):
            if id2 in ["ResNum", "ResName"]:
                continue
            index2 += 1
            if id1 == id2:
                continue
            similarity1 = 0
            similarity2 = 0
            for res in range(n_res):
                c1a, c1b = clean_list([contacts1[res]], delimiter=",", format="bool")
                c2a, c2b = clean_list([contacts2[res]], delimiter=",", format="bool")

                if c1a == c2a:
                    similarity1 += 0.5
                if c1b == c2b:
                    similarity1 += 0.5
                if c1a == c2b:
                    similarity2 += 0.5
                if c1b == c2a:
                    similarity2 += 0.5

            similarity = max(similarity1, similarity2)
            similarity = similarity / n_res
            #print(id1, id2, similarity)
            sm_ssd.loc[len(sm_ssd)] = id1, id2,index1,index2, similarity
            progress.add()
    print(sm_ssd)
    sm_ssd.to_csv(os.path.join(root.dataframes, '{}_sm_ssd.csv'.format(reference.name)),header=False, index=False)




def cc_analysis(reference, dimensions=3, force =False   ):
    sprint("CC analysis for {}".format(reference.name))

    if "{}_cc_output.csv".format(reference.name) in os.listdir(root.dataframes) and not force:
        print1("Skipping CC analysis for {}".format(reference.name))
        return

    sm_ssd_path = os.path.join(root.dataframes, "{}_sm_ssd.csv".format(reference.name))
    sm_ssd = pd.read_csv(sm_ssd_path, index_col=0)
    print(sm_ssd)
    if len(sm_ssd) >= 6:
        ## CC analysis
        import subprocess

        cc_path = os.path.join(root.scripts, "cc_analysis.py")
        cc_line = ["python", cc_path, str(dimensions), sm_ssd_path]  # Currently uses correlation matrix
        print1("cc_line:")
        print2(" ".join(cc_line))
        cc_std = subprocess.run(cc_line,
                                capture_output=True,
                                text=True)
        print(cc_std.stdout)
        print(cc_std.stderr)
        # print(cc_std.stdout.split("\n"))
        out = []
        for l in cc_std.stdout.split("\n"):
            l = l.strip().split(" ")
            for i in l:
                if i == "":
                    l.remove(i)
            # print(pd.to_numeric(l))
            if len(l) == 3:
                l.append("0")
                out.append(l)
            elif len(l) == 4:
                out.append(l)
            # print(l)
        # print("-")
        # print(out)
        # print("-")

        cc_out = pd.DataFrame(out, columns=["index", "1", "2", "3"])
        # print(cc_out)


        for i in cc_out.columns:
            if i == "index":
                cc_out["index"] = pd.to_numeric(cc_out["index"], downcast='integer', errors='coerce')
            else:
                cc_out[i] = pd.to_numeric(cc_out[i], downcast='float', errors='coerce')

        cc_out_path = os.path.join(root.dataframes, "{}_cc_output.csv".format(reference.name))
        cc_out.to_csv(cc_out_path, index=False, header=True)

        print1("CC output:")
        print2(cc_out_path)
        print(cc_out, "\n")



def plot_cc(reference, force=False, dimensions = 3):
    print1("Plotting 2D: {}".format(reference.name))
    import matplotlib.pyplot as plt
    from adjustText import adjust_text

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    root["cc"] = "images/cc"
    figure_path = os.path.join(root.cc, "{}_cc.png")

    if dimensions == 2:
        ax.scatter(cc_out["1"], cc_out["2"], c=cc_out["colour"])
    elif dimensions == 3:
        ax.scatter(cc_out["3"], cc_out["2"],c=cc_out["colour"])
    ax.scatter(0, 0, color='red')
    for file in os.listdir(folder_path):
        #print(file)
        if "cluster_centres" in file:
            cluster_centres = pd.read_csv(os.path.join(folder_path, file))
            #print(cluster_centres)
            break
        else:
            continue
        print("centres (.csv) not found")
        quit()

    ax.scatter(cluster_centres["2"],cluster_centres["1"], color="black")

    centres = []
    for centre in cluster_centres.itertuples():
        #print(centre[0])
        #print(centre[1:])
        centres.append(centre[1:])
        #print(centre[0],centre[2],centre[3])


    lines = []
    #print("points")
    n = 0

    for point in cc_out[["1","2","3"]].itertuples():
        #print(point)
        point = point[1:]
        #print("point:", point)
        closest = get_closest_point(point,centres)
        #print("closest:",closest)


        lines.append(points_to_line(closest,point))
        #print("point to line:", points_to_line(closest,point))
        #lines[str(n)] = line

    #print("lines")
    #print(lines)
    for line in lines:
        #print(line)
        ax.plot(line[2],line[1], c="black")
    #scripts.Mpl.plot_lines(lines, ax)
    #ax.plot(data=lines(lines[2],lines[1]))


    method = kwargs["method"]
    mode = kwargs["plot_mode"]

    texts = []
    for centre in cluster_centres.itertuples():
        texts.append(ax.annotate(centre[0], (centre[3],centre[2]), size=10))
    #print(cc_out.iloc[0])
    for n in range(len(cc_out)):
        index = cc_out["index"][n]
        # print(cc_out["index"])
        if dimensions == 3:
            X = cc_out["3"][n]
            Y = cc_out["2"][n]
        else:
            X = cc_out["1"][n]
            Y = cc_out["2"][n]
        if kwargs["labels"]:
            texts.append(
                ax.annotate(cc_out.loc[:, "cluster"][n], (X, Y)))  # ax.annotate(labels[n]
        '''else:
            texts.append(
                ax.annotate(df.loc[:, "ID"][n], (X, Y)))  # ax.annotate(labels[n]'''

    n_clusters_str = add_front_0(kwargs["n_clusters"])

    ax.set_title("number of clusters: {}".format(n_clusters_str))
    adjust_text(texts, autoalign='y',
                only_move={'points': 'y', 'text': 'y'}, force_points=0.15,
                arrowprops=dict(arrowstyle="->", color='blue', lw=0.5))
    fig.tight_layout()
    print2("Saving at {}".format(folder_path))
    fig_name = "cc_plot_2D_{}_{}_{}.png".format(len(df), method, n_clusters_str )
    print2(fig_name)
    if "show" in mode:
        fig.show(block=True)
    elif "skip" not in mode:
        fig.show(block=False)
    if "save"in mode:
        fig.savefig(os.path.join(folder_path, fig_name))
        fig.savefig(os.path.join(figure_dir_folder,fig_name))
    #plt.close()











if __name__ == "__main__":

    FORCE_SM = False
    FORCE_CC = True
    FORCE_PLOT = False

    FORCE_ALL = False



    import setup
    from Globals import root, local, vars

    if FORCE_ALL:
        FORCE_SM = True
        FORCE_CC = True

    from imports import import_references
    references = import_references()

    tprint("Similarity analysis")
    for reference in references:
        generate_sm(reference, force=FORCE_SM)
    eprint("Similarity analysis")

    tprint("CC analysis")
    for reference in references:
        cc_analysis(reference, force=FORCE_CC)
        #plot_cc(reference, force=FORCE_PLOT)
    eprint("CC analysis")


    eprint("Done")


                
            