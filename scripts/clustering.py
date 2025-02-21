import os
from itertools import count

from numpy.ma.extras import average

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

    progress = ProgressBar(n_dimers)
    index1 = 0
    from maths import difference_between_boolean_pairs
    for id1, contacts1 in zip(sasas_df.columns, sasas_df._iter_column_arrays()):
        if id1 in ["ResNum", "ResName"]:
            continue
        progress.add(info="{}". format(id1), show_time = True)
        index1 += 1
        index2 = 0
        for id2, contacts2 in zip(sasas_df.columns, sasas_df._iter_column_arrays()):
            if id2 in ["ResNum", "ResName"]:
                continue
            index2 += 1
            if index2 <= index1:
                continue
            diffX = [0,0]
            diffx = [0,0]
            for res in range(n_res):
                c1a, c1b = clean_list([contacts1[res]], delimiter=",", format="bool")
                c2a, c2b = clean_list([contacts2[res]], delimiter=",", format="bool")

                resX, resx =difference_between_boolean_pairs(c1a,c1b,c2a,c2b)
                diffX[0] += resX[0]
                diffX[1] += resX[1]
                diffx[0] += resx[0]
                diffx[1] += resx[1]

            if diffX[0] != 0:
                diffX = diffX[0]/diffX[1]
            else:
                diffX = 0

            if diffx[0] != 0:
                diffx = diffx[0]/diffx[1]
            else:
                diffx = 0

            similarity = max(diffX, diffx)
            #similarity = similarity / n_res
            #print(id1, id2, round(similarity,2))
            sm_ssd.loc[len(sm_ssd)] = id1, id2,index1,index2, similarity
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
        if len(l)-1  == 2:
            l.append("0")
            out.append(l)
        elif len(l) >= 3:
            out.append(l)
        # print(l)
    # print("-")
    # print(out)
    # print("-")

    #print(out)
    cols = ["index"]
    for n in range(dimensions):
        cols.append(str(n+1))
    cc_out = pd.DataFrame(out, columns=cols)


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


def clusterize_cc(reference, force=False, n_clusters = 20, dimensions=3):
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

    if len(cc_out) < 30 and len(cc_out) > 6:
        n_clusters = 5
    model = KMeans(n_clusters=n_clusters, random_state=6, algorithm="elkan")
    cols = []
    for n in range(dimensions):
        cols.append(str(n + 1))
    model.fit(cc_out.loc[:, cols])
    #pred = model.fit(cc_out.loc[:, ["1", "2", "3"]])

    #print(model.cluster_centers_)
    for n in range(len(cc_out)):
        cluster = model.labels_[n]
        cc_out.loc[n+1, "cluster"] = cluster
        new_colour = "".join(["C", str(cluster)])

        if new_colour == "C":
            new_colour = "black"
        cc_out.loc[n+1, "colour"] = new_colour

    cluster_centres_path = os.path.join(root.dataframes, "{}_cluster_centres.csv".format(reference.name))
    cluster_centres_df = pd.DataFrame(model.cluster_centers_)
    #print(cluster_centres_df)
    cluster_centres_df.to_csv(cluster_centres_path,header=None)
    cc_out.to_csv(os.path.join(root.dataframes,cc_clustered_name))





def plot_cc(reference, force=False, dimensions = 3, labels = False, labels_centres=True, adjust=False):
    print1("Plotting 2D: {}".format(reference.name))

    root["cc"] = "images/cc"
    figure_name ="{}_cc.png".format(reference.name)
    figure_path = os.path.join(root.cc,figure_name )
    if figure_name in os.listdir(root.cc) and not force:
        print2("Figure {} already exists".format(figure_name))
        return
    if dimensions > 3:
        print2("Plotting more than 3 dimensions not currently supported")
        return
    import matplotlib.pyplot as plt


    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    try:
        cc_out = pd.read_csv(os.path.join(root.dataframes, "{}_cc_clustered.csv".format(reference.name)))
    except:
        print1("Error parsing {}, might be empy".format(figure_path))
        return

    if dimensions == 2:
        ax.scatter(cc_out["1"], cc_out["2"], c=cc_out["colour"])
    elif dimensions == 3:
        ax.scatter(cc_out["3"], cc_out["2"],c=cc_out["colour"])
    ax.scatter(0, 0, color='red')

    texts = []

    if "{}_cluster_centres.csv".format(reference.name) in os.listdir(root.dataframes):
        print2("Plotting cluster centres")
        from maths import get_closest_point, points_to_line
        cluster_centres = pd.read_csv(os.path.join(root.dataframes, "{}_cluster_centres.csv".format(reference.name)),index_col=0,header=None)
        print(cluster_centres.columns)
        ax.scatter(cluster_centres.loc[:,3].values, cluster_centres.loc[:,2].values, color="black", marker=".")

        centres = []
        for centre in cluster_centres.itertuples():
            print(centre[0])
            print(centre[1:])
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
        if labels_centres:
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
    if adjust:
        from adjustText import adjust_text
        adjust_text(texts, autoalign='y',
                    only_move={'points': 'y', 'text': 'y'}, force_points=0.15,
                    arrowprops=dict(arrowstyle="->", color='blue', lw=0.5))
    fig.tight_layout()
    print2("Saving at {}".format(figure_path))
    fig.savefig(figure_path, dpi=300)





def clustering(FORCE_ALL = False, FORCE_SM = True, FORCE_CC = True, FORCE_CLUSTER = True, FORCE_PLOT = True, DIMENSIONS = 3, ONLY_GR=False):


    if FORCE_ALL:
        FORCE_SM = True
        FORCE_CC = True
        FORCE_CLUSTER = True
        FORCE_PLOT = True

    from imports import import_references
    references = import_references()
    for reference in references:
        if ONLY_GR and reference.name == "GR":
            references = [reference]

    from dataframes import load_failed_dfs
    load_failed_dfs()

    tprint("Similarity analysis")
    for reference in references:
        generate_sm(reference, force=FORCE_SM)
    eprint("Similarity analysis")

    tprint("CC analysis")
    for reference in references:
        cc_analysis(reference, force=FORCE_CC, dimensions=DIMENSIONS)
        clusterize_cc(reference, force=FORCE_CLUSTER, dimensions=DIMENSIONS)
        plot_cc(reference,labels=False, labels_centres=True, force=FORCE_PLOT, dimensions=DIMENSIONS)
    eprint("CC analysis")

    tprint("Comparing to Eva")
    gr_df = pd.read_csv(os.path.join(root.dataframes, "GR_cc_clustered.csv"))
    if "BALL_SIZE" in vars:
        scores = calculate_scores_GR(gr_df, vars.BALL_SIZE+"_"+DIMENSIONS)
    else:
        scores = calculate_scores_GR(gr_df,str(DIMENSIONS))
    print1("Scores: cc: {}, eva: {}".format(scores[0], scores[1]))
    eprint("Compared successfully")


    eprint("Done")

def get_cluster_score(df, primary, secondary):
    scores = []
    for cluster in df[primary].unique():
        print1(cluster)
        f_df = df[df[primary] == cluster]
        cluster_total = len(f_df)
        counts = []
        group_list = f_df[secondary].tolist()
        groups = set(group_list)

        for unique_group in groups:
            c = group_list.count(unique_group)
            counts.append(c)
            print2(unique_group, c)
        maximum = max(counts)
        score = maximum / cluster_total
        scores.append(score)
        print3("Score:", round(score, 2))
    av_score = sum(scores) / len(scores)
    return av_score, scores



def calculate_scores_GR(df, name="undefined"):
    if "scores_df.csv" in os.listdir(root.dataframes):
        scores_df = pd.read_csv(os.path.join(root.dataframes, "scores_df.csv"), index_col=0)
    else:
        scores_df = pd.DataFrame(columns = ["label", "cc_score", "eva_score", "cc_values", "eva_values" ])

    print(scores_df)

    cc_scores = get_cluster_score(df, primary="cluster", secondary="group")
    eva_scores = get_cluster_score(df, primary="group", secondary="cluster")
    scores_df.loc[len(scores_df)]= [name, cc_scores[0], eva_scores[0], cc_scores[1], eva_scores[1]]
    print(scores_df)
    scores_df.to_csv(os.path.join(root.dataframes, "scores_df.csv"))
    return cc_scores, eva_scores






if __name__ == "__main__":




    import setup
    from Globals import root, local, vars




    '''#### DIMENSION TESTING ###
    dimensions = [3,4,5,6,7,8,9,10]
    for n in dimensions:
        clustering(FORCE_SM = False,
                   FORCE_CC = True,
                   FORCE_CLUSTER = True,
                   FORCE_PLOT=False,
                   DIMENSIONS=n,
                   ONLY_GR = True
                   )
    #### DIMENSION TESTING ###'''

    clustering(FORCE_ALL=False,
               FORCE_SM=False,
               FORCE_CC=True,
               FORCE_CLUSTER=True,
               FORCE_PLOT=True,
               DIMENSIONS=5,
               )
    from github import automatic_push_to_branch
    #automatic_push_to_branch(target="auto")
            