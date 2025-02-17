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
    if n_dimers <2:
        print1("Not enough dimers in {} dataframe".format(reference.name))
        return
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
            if index2 <= index1:
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
            progress.add(info="time")
    print(sm_ssd)
    sm_ssd.to_csv(os.path.join(root.dataframes, '{}_sm_ssd.csv'.format(reference.name)),header=False, index=False)




def cc_analysis(reference, dimensions=3, force =False):
    sprint("CC analysis for {}".format(reference.name))

    if "{}_cc_output.csv".format(reference.name) in os.listdir(root.dataframes) and not force:
        print1("Skipping CC analysis for {}".format(reference.name))
        return

    sm_ssd_path = os.path.join(root.dataframes, "{}_sm_ssd.csv".format(reference.name))
    try:
        sm_ssd = pd.read_csv(sm_ssd_path, index_col=None, header=None)
    except:
        print1("Error parsing {}, might be empy".format(sm_ssd_path))
        return
    #print(sm_ssd)
    if len(sm_ssd.columns) != 5:
        print1("SM Must be 5 columns wide (id1, id2, index1, index2, similarity)")
        print2("Current:", len(sm_ssd.columns))
        return

    if len(sm_ssd) < 6:
        print1("Not enough values for cc analysis, min: 6, current: {}".format(len(sm_ssd)))
        return
    ## CC analysis
    import subprocess

    cc_path = os.path.join(root.scripts, "cc_analysis.py")
    cc_line = ["python", cc_path, str(dimensions), sm_ssd_path]  # Currently uses correlation matrix
    print1("cc_line:")
    print2(" ".join(cc_line))

    cc_std = subprocess.run(cc_line,
                            capture_output=True,
                            text=True)
    #print(cc_std.stdout)
    #print(cc_std.stderr)
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

    #print(out)
    cc_out = pd.DataFrame(out, columns=["index", "1", "2", "3"])


    for i in cc_out.columns:
        if i == "index":
            cc_out["index"] = pd.to_numeric(cc_out["index"], downcast='integer', errors='coerce')
        else:
            cc_out[i] = pd.to_numeric(cc_out[i], downcast='float', errors='coerce')
    print(cc_out)
    cc_out["id"] = pd.read_csv(os.path.join(root.dataframes, "{}_sasas.csv".format(reference.name)), index_col=0).columns[2:]

    classified_df = pd.read_csv(os.path.join(root.dataframes, "classified_df.csv"), index_col=0)
    classified_ids = classified_df["ID"].values
    groups = []
    for row in cc_out.itertuples():
        if row.id in classified_ids:
            groups.append(classified_df[classified_df["ID"]==row.id].Best_Match.values[0])
        else:
            groups.append("na")
    cc_out["group"] = groups

    cc_out_path = os.path.join(root.dataframes, "{}_cc_output.csv".format(reference.name))
    cc_out.to_csv(cc_out_path, index=False, header=True)

    print1("CC output:")
    print2(cc_out_path)
    print(cc_out, "\n")


def clusterize_cc(reference, force=False, n_clusters = 20):
    cc_clustered_name = "{}_cc_clustered.csv".format(reference.name)
    if cc_clustered_name in os.listdir(root.dataframes) and not force:
        print1("Skipping CC clustering for {}".format(reference.name))
        return

    sprint("Clustering {}".format(reference.name))
    from sklearn.cluster import KMeans

    cc_out_path = os.path.join(root.dataframes, "{}_cc_output.csv".format(reference.name))
    try:
        cc_out = pd.read_csv(cc_out_path, index_col=0)
    except:
        print1("Error parsing {}, might be empy".format(cc_out_path))
        return


    model = KMeans(n_clusters=n_clusters)
    model.fit(cc_out.loc[:, ["1", "2", "3"]])
    pred = model.fit(cc_out.loc[:, ["1", "2", "3"]])

    print(model.cluster_centers_)
    for n in range(len(cc_out)):
        cluster = model.labels_[n]
        cc_out.loc[n, "cluster"] = cluster
        new_colour = "".join(["C", str(cluster)])
        # print(new_colour)
        if new_colour == "C":
            new_colour = "black"
        cc_out.loc[n, "colour"] = new_colour

    cluster_centres_path = os.path.join(root.dataframes, "{}_cluster_centres.csv".format(reference.name))
    cluster_centres_df = pd.DataFrame(model.cluster_centers_)
    cluster_centres_df.to_csv(cluster_centres_path, index=False)
    cc_out.to_csv(os.path.join(root.dataframes,cc_clustered_name))





def plot_cc(reference, force=False, dimensions = 3, labels = True):
    print1("Plotting 2D: {}".format(reference.name))

    root["cc"] = "images/cc"
    figure_path = os.path.join(root.cc, "{}_cc.png")
    if figure_path in os.listdir(root.cc) and not force:
        return

    import matplotlib.pyplot as plt
    from adjustText import adjust_text

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    cc_out = pd.read_csv(os.path.join(root.dataframes, "{}_cc_output.csv".format(reference.name)))


    if dimensions == 2:
        ax.scatter(cc_out["1"], cc_out["2"], c=cc_out["colour"])
    elif dimensions == 3:
        ax.scatter(cc_out["3"], cc_out["2"],c=cc_out["colour"])
    ax.scatter(0, 0, color='red')

    texts = []

    if "{}_cluster_centres.csv".format(reference.name) in os.listdir(root.dataframes):
        from maths import get_closest_point, points_to_line
        cluster_centres = pd.read_csv(os.path.join(root.dataframes, "{}_cluster_centres.csv".format(reference.name)))
        #print(cluster_centres)
        ax.scatter(cluster_centres["2"], cluster_centres["1"], color="black")

        centres = []
        for centre in cluster_centres.itertuples():
            # print(centre[0])
            # print(centre[1:])
            centres.append(centre[1:])
            # print(centre[0],centre[2],centre[3])

        lines = []
        # print("points")
        n = 0

        for point in cc_out[["1", "2", "3"]].itertuples():
            # print(point)
            point = point[1:]
            # print("point:", point)
            closest = get_closest_point(point, centres)
            # print("closest:",closest)

            lines.append(points_to_line(closest, point))
            # print("point to line:", points_to_line(closest,point))
            # lines[str(n)] = line

        # print("lines")
        # print(lines)
        for line in lines:
            # print(line)
            ax.plot(line[2], line[1], c="black")
        # scripts.Mpl.plot_lines(lines, ax)
        # ax.plot(data=lines(lines[2],lines[1]))
        for centre in cluster_centres.itertuples():
            texts.append(ax.annotate(centre[0], (centre[3],centre[2]), size=10))

    else:
        print2("centres ({}_cluster_centres.csv) not found".format(reference.name))





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
        if labels:
            texts.append(
                ax.annotate(cc_out.loc[:, "cluster"][n], (X, Y)))  # ax.annotate(labels[n]
        '''else:
            texts.append(
                ax.annotate(df.loc[:, "ID"][n], (X, Y)))  # ax.annotate(labels[n]'''


    ax.set_title("CC analysis + Kmeans for {}".format(reference.name))
    adjust_text(texts, autoalign='y',
                only_move={'points': 'y', 'text': 'y'}, force_points=0.15,
                arrowprops=dict(arrowstyle="->", color='blue', lw=0.5))
    fig.tight_layout()
    print2("Saving at {}".format(figure_path))
    fig.savefig(figure_path, dpi=300)











if __name__ == "__main__":

    FORCE_SM = False
    FORCE_CC = False
    FORCE_CLUSTER = True
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
        clusterize_cc(reference, force=FORCE_CLUSTER)
        #plot_cc(reference, force=FORCE_PLOT)
    eprint("CC analysis")


    eprint("Done")


                
            