import os, sys

from pyMol import pymol_delete
from utilities import *
from Globals import root, local, vars
import numpy as np
import pandas as pd



class Cluster:
    def __init__(self, ref_name, df, n, column):

        self.name = ref_name
        self.df = df
        self.n = n
        self.column = column
        self.eva = False
        if column == "Best_Match":
            self.eva = True

    def __repr__(self):
        return f"Cluster({self.n}, {self.name})"


def get_clusters(df, column, ref_name, **kwargs):
    #df = pd.read_csv(df_path)
    clusters = []
    for c in set(df[column]):
        subset = df[df[column] == c]
        clusters.append(Cluster(ref_name, subset, c, column, **kwargs))
    return clusters



def generate_sm(reference, force=False, subfolder=None, in_path = None, use_csv = True):
    print1("Generating SM for {}".format(reference.name))

    root["sms"] = "dataframes/clustering/sms"
    if subfolder is None:
        sm_path = os.path.join(root.sms, '{}.csv'.format(reference.name))
        if "{}.csv".format(reference.name) in os.listdir(root.sms) and not force:
            print2("Skipping SM generation for {}".format(reference.name))
            return sm_path
        if use_csv:
            contacts_df = pd.read_csv(os.path.join(root.contacts, "{}.csv".format(reference.name)))
        else:
            contacts_df = reference.contacts_df
    else:
        assert in_path is not None
        name = os.path.basename(in_path).split(".")[0]
        subfolder = subfolder.format("sms")
        os.makedirs(os.path.join(root.sms, subfolder), exist_ok=True)
        sm_path = os.path.join(root.sms, subfolder, name + ".csv")
        if name+".csv" in os.listdir(os.path.join(root.sms, subfolder)) and not force:
            print2("Skipping SM generation for {}".format(name))
            return sm_path
        contacts_df = pd.read_csv(in_path)

    print(contacts_df)


    sm_ssd = pd.DataFrame(columns=["dimer1", "dimer2", "index1", "index2", "similarity", "diffX", "diffx"])
    n_dimers = len(contacts_df.columns)
    #print(contacts_df.columns)
    if n_dimers <2:
        print1("Not enough dimers in {} dataframe".format(reference.name))
        return None
    n_res = len(contacts_df)

    progress = ProgressBar(n_dimers)
    index1 = 0
    from maths import difference_between_boolean_pairs
    for id1, contacts1 in zip(contacts_df.columns, contacts_df._iter_column_arrays()):
        if id1 in ["ResNum", "ResName"]:
            continue
        progress.add(info="{}". format(id1), show_time = True)
        index1 += 1
        index2 = 0
        for id2, contacts2 in zip(contacts_df.columns, contacts_df._iter_column_arrays()):
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

            dX = diffX
            dx = diffx

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
            sm_ssd.loc[len(sm_ssd)] = id1, id2,index1,index2, similarity, dX, dx
    #print(sm_ssd)

    if subfolder is None:
        reference.sms_df = sm_ssd

    sm_ssd.to_csv(sm_path, index=False)

    print(sm_ssd)
    return sm_path




def cc_analysis(reference, dimensions=3, force =False, subfolder = None, in_path = None, use_csv = True):
    print1("CC analysis for {}".format(reference.name))

    root["ccs"] = "dataframes/clustering/ccs"
    if subfolder is None:
        ccs_path = os.path.join(root.ccs, '{}.csv'.format(reference.name))
        if "{}.csv".format(reference.name) in os.listdir(root.ccs) and not force:
            print2("Skipping CC analysis for {}".format(reference.name))
            return ccs_path
        if use_csv:
            sm_ssd = pd.read_csv(os.path.join(root.sms, "{}.csv".format(reference.name)))
        else:
            sm_ssd = reference.sms_df

    else:
        assert in_path is not None
        name = os.path.basename(in_path).split(".")[0]
        subfolder = subfolder.format("ccs")
        os.makedirs(os.path.join(root.ccs, subfolder), exist_ok=True)
        ccs_path = os.path.join(root.ccs, subfolder, name + ".csv")
        if name+".csv" in os.listdir(os.path.join(root.ccs, subfolder)) and not force:
            print2("Skipping CC analysis for {}".format(name))
            return ccs_path
        sm_ssd = pd.read_csv(in_path)

    print(sm_ssd)

    sm_path = os.path.join(local.temp, "sm_ssd_for_cc.csv")
    sm_ssd.iloc[:, 2:5].to_csv(sm_path, index=False, header=False)

    if len(sm_ssd.columns) != 7:
        print1("SM Must be 7 columns wide (id1, id2, index1, index2, similarity, diffX, diffx)")
        print2("Current:", len(sm_ssd.columns))
        return None

    if len(sm_ssd) < 6:
        print1("Not enough values for cc analysis, min: 6, current: {}".format(len(sm_ssd)))
        return None
    ## CC analysis
    import subprocess

    cc_path = os.path.join(root.scripts, "cc_analysis.py")
    cc_line = ["python", cc_path, str(dimensions), sm_path]  # Currently uses correlation matrix
    print1("cc_line:")
    print2(" ".join(cc_line))
    progress = ProgressBar(1)
    cc_std = subprocess.run(cc_line,
                            capture_output=True,
                            text=True)
    progress.add()
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
    #print("---")
    #print(out)
    #print("---")
    if subfolder is None:
        reference.cc_raw = out


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
    print("cc_out:")
    print(cc_out)
    #print("ids:")
    if subfolder is None:
        cc_out["id"] = pd.read_csv(os.path.join(root.contacts, "{}.csv".format(reference.name))).columns[2:]
    else:
        print(name)
        cc_out["id"] = pd.read_csv(os.path.join(root.contacts, "contacts_{}/{}.csv".format(reference.name, name))).columns[2:]
    print(cc_out)

    if "{}.csv".format(reference.name) in os.listdir(root.classified):
        classified_df = pd.read_csv(os.path.join(root.classified, "{}.csv".format(reference.name)))
        #print(classified_df)
        classified_ids = classified_df["ID"].values
        groups = []
        similarities = []
        inverses = []
        for row in cc_out.itertuples():
            if row.id in classified_ids:
                groups.append(classified_df[classified_df["ID"]==row.id].Best_Match.values[0])
                similarities.append(classified_df[classified_df["ID"]==row.id].Similarity.values[0])
                inverses.append(classified_df[classified_df["ID"] == row.id].Inverse.values[0])
            else:
                groups.append("na")
                similarities.append("na")
                inverses.append("na")
        cc_out["group"] = groups
        cc_out["similarity"] = similarities
        cc_out["inverse"] = inverses

    #cc_out_path = os.path.join(root.ccs, "{}.csv".format(reference.name))
    cc_out.to_csv(ccs_path, index=False, header=True)
    if subfolder is not None:
        reference.ccs_df = cc_out

    print1("CC output:")
    print2(ccs_path)
    print(cc_out, "\n")
    return ccs_path



def clusterize_cc(reference, force=False, n_clusters = 20, dimensions=3, subfolder=None, in_path=None, use_csv =True):
    print1("Clustering of {}".format(reference.name))


    root["clustered"] = "dataframes/clustering/clustered"
    if subfolder is None:
        clustered_path = os.path.join(root.ccs, '{}.csv'.format(reference.name))
        centres_path = os.path.join(root.ccs, '{}_centres.csv'.format(reference.name))
        if "{}.csv".format(reference.name) in os.listdir(root.clustered) and not force:
            print2("Skipping clustering for {}".format(reference.name))
            return clustered_path
        if use_csv:
            cc_out = pd.read_csv(os.path.join(root.ccs, "{}.csv".format(reference.name)))
        else:
            cc_out = reference.ccs_df
    else:
        assert in_path is not None
        name = os.path.basename(in_path).split(".")[0]
        subfolder = subfolder.format("clustered")
        os.makedirs(os.path.join(root.clustered, subfolder), exist_ok=True)
        clustered_path = os.path.join(root.clustered, subfolder, name + ".csv")
        centres_path = os.path.join(root.clustered, subfolder, name + "_centres.csv")
        if name+".csv" in os.listdir(os.path.join(root.clustered, subfolder)) and not force:
            print2("Skipping clustering for {}".format(name))
            return clustered_path
        cc_out = pd.read_csv(in_path)

    print(cc_out)
    from sklearn.cluster import KMeans
    #if len(cc_out) < 30 and len(cc_out) > 6:
        #n_clusters = 5
    model = KMeans(n_clusters=n_clusters, random_state=6, algorithm="elkan")
    cols = []
    for n in range(dimensions):
        cols.append(str(n + 1))
    model.fit(cc_out.loc[:, cols])
    #pred = model.fit(cc_out.loc[:, ["1", "2", "3"]])

    #print(model.cluster_centers_)
    for n in range(len(cc_out)):
        cluster = model.labels_[n]
        cc_out.loc[n, "cluster"] = cluster
        new_colour = "".join(["C", str(cluster)])

        if new_colour == "C":
            new_colour = "black"
        cc_out.loc[n, "colour"] = new_colour

    cluster_centres = pd.DataFrame(model.cluster_centers_)
    print(cluster_centres)
    print(cc_out)
    cluster_centres.to_csv(centres_path, index=False)
    cc_out.to_csv(clustered_path)

    if subfolder is None:
        reference.cluster_centres = cluster_centres
        reference.clustered_df = cc_out
    return clustered_path





def plot_cc(reference, force=True, dimensions = 3, labels = False, labels_centres=True, adjust=False, subfolder=None,
            in_path=None, use_csv = True, plot_centres=True, subset = None, pca=False):
    print1("Plotting: {}, Dimensions: {}".format(reference.name, dimensions))

    if dimensions > 3:
        print2("Plotting more than 3 dimensions not currently supported")
        return None

    root["cc_figs"] = "images/cc_figs"
    root["pca_figs"] = "images/pca_figs"
    if pca:
        fig_folder = root.pca_figs
        fig_folder_name = "pca_figs"
        clustered_folder = root.clustered_pcas
    else:
        fig_folder = root.cc_figs
        fig_folder_name = "cc_figs"
        clustered_folder = root.clustered
    if subfolder is None:
        fig_path = os.path.join(fig_folder, '{}.png'.format(reference.name))
        if "{}.png".format(reference.name) in os.listdir(fig_folder) and not force:
            print2("Skipping plotting for {}".format(reference.name))
            return fig_path
        if use_csv:
            cc_out = pd.read_csv(os.path.join(clustered_folder, "{}.csv".format(reference.name)))
        else:
            cc_out = reference.clustered_df

    else:
        assert in_path is not None
        name = os.path.basename(in_path).split(".")[0]
        clustered_subfolder = subfolder.format("clustered")
        subfolder = subfolder.format(fig_folder_name)
        os.makedirs(os.path.join(fig_folder, subfolder), exist_ok=True)
        fig_path = os.path.join(fig_folder, subfolder, name + ".png")
        if name+".png" in os.listdir(os.path.join(fig_folder, subfolder)) and not force:
            print2("Skipping plotting for {}".format(name))
            return fig_path
        cc_out = pd.read_csv(in_path, index_col=0)
        if subset is not None:
            cc_out = cc_out.query(" | ".join(["{} == {}".format("cluster", n) for n in subset]))

            #cc_out = cc_out[cc_out["cluster"] == subset]
            cc_out.reset_index(drop=True, inplace=True)

    #print(cc_out)

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)


    if dimensions == 2:
        ax.scatter(cc_out["1"], cc_out["2"], c=cc_out["colour"])
    elif dimensions == 3:
        ax.scatter(cc_out["3"], cc_out["2"],c=cc_out["colour"])
    ax.scatter(0, 0, color='red')

    texts = []

    if plot_centres and subset is None:
        print2("Plotting cluster centres")
        if subfolder is None:
            cluster_centres = reference.cluster_centres
        else:
            cluster_centres = pd.read_csv(os.path.join(clustered_folder,clustered_subfolder, "{}_centres.csv".format(name)))

        from maths import get_closest_point, points_to_line

        #print(cluster_centres)
        #print(cluster_centres.columns)
        ax.scatter(cluster_centres["2"].values, cluster_centres["1"].values, color="black", marker=".")

        centres = []
        for centre in cluster_centres.itertuples():
            #print(centre[0])
            #print(centre[1:])
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


    #print(cc_out)
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
                ax.annotate(cc_out.loc[:, "id"][n], (X, Y)))  # ax.annotate(labels[n]
        '''else:
            texts.append(
                ax.annotate(df.loc[:, "ID"][n], (X, Y)))  # ax.annotate(labels[n]'''

    if subfolder is None:
        ax.set_title("CC analysis + Kmeans for {} (n= {})".format(reference.name, len(cc_out)))
    else:
        if subset is None:
            ax.set_title("CC analysis + Kmeans for {} ({}) (n= {})".format(reference.name, name, len(cc_out)))
        else:
            ax.set_title("CC analysis + Kmeans for {} ({}) CLUSTER {} (n= {})".format(reference.name, name, subset, len(cc_out)))
    if adjust:
        from adjustText import adjust_text
        adjust_text(texts, autoalign='y',
                    only_move={'points': 'y', 'text': 'y'}, force_points=0.15,
                    arrowprops=dict(arrowstyle="->", color='blue', lw=0.5))
    fig.tight_layout()
    ax.set_aspect('equal')
    if subset is None:
        print2("Saving at {}".format(fig_path))
        fig.savefig(fig_path, dpi=300)
        return fig_path
    else:
        plt.show(block=True)


def cluster(reference, FORCE_ALL=False, DIMENSIONS = 3, score_id = "", thread = False):
    if FORCE_ALL:
        FORCE_SM = True
        FORCE_CC = True
        FORCE_CLUSTER = True
        FORCE_PLOT = True
    else:
        FORCE_SM = False
        FORCE_CC = False
        FORCE_CLUSTER = False
        FORCE_PLOT = False


    reference.sm_path = generate_sm(reference, force=FORCE_SM)
    reference.cc_path = cc_analysis(reference, force=FORCE_CC, dimensions=DIMENSIONS)
    reference.clustered_path = clusterize_cc(reference, force=FORCE_CLUSTER, dimensions=DIMENSIONS)
    reference.plot_path = plot_cc(reference, labels=False, labels_centres=True, force=FORCE_PLOT, dimensions=DIMENSIONS)
    return
    if reference.name == "GR":
        tprint("Comparing to Eva")
        df_cc = pd.read_csv(os.path.join(root.clustered, "GR_cc.csv"))
        scores = calculate_scores_GR(df_cc, score_id+str(DIMENSIONS))
        print1("Scores: cc: {}, eva: {}".format(scores[0], scores[1]))
        eprint("Compared successfully")




def get_cluster_score(df, primary, secondary):
    scores = []
    for cluster in df[primary].unique():
        #print1(cluster)
        f_df = df[df[primary] == cluster]
        cluster_total = len(f_df)
        counts = []
        group_list = f_df[secondary].tolist()
        groups = set(group_list)

        for unique_group in groups:
            c = group_list.count(unique_group)
            counts.append(c)
            #print2(unique_group, c)
        maximum = max(counts)
        score = maximum / cluster_total
        scores.append(score)
        #print3("Score:", round(score, 2))
    av_score = sum(scores) / len(scores)
    return av_score, scores



def calculate_scores_GR(df, name="undefined", save = True):
    if "scores_df.csv" in os.listdir(root.dataframes):
        scores_df = pd.read_csv(os.path.join(root.dataframes, "scores_df.csv"), index_col=0)
    else:
        scores_df = pd.DataFrame(columns = ["label", "cc_score", "eva_score", "cc_values", "eva_values" ])
    if not save:
        name = "current"
    cc_scores = get_cluster_score(df, primary="cluster", secondary="group")
    eva_scores = get_cluster_score(df, primary="group", secondary="cluster")
    scores_df.loc[len(scores_df)]= [name, cc_scores[0], eva_scores[0], cc_scores[1], eva_scores[1]]
    print(scores_df[["label", "cc_score", "eva_score"]])
    if save:
        scores_df.to_csv(os.path.join(root.dataframes, "scores_df.csv"))
    return cc_scores, eva_scores



def compare_contacts(reference, force = False):
    print1("Comparing contacts for GR / force= {}".format(force))

    classified_path = os.path.join(root.classified, "GR.csv")
    print(classified_path)
    if os.path.exists(classified_path) and not force:
        vars.clustering["classified"][reference.name] = pd.read_csv(classified_path)
        if len(vars.clustering["classified"][reference.name]) != 0:
            return vars.clustering["classified"][reference.name]

    vars["clustering"]["classified"][reference.name] = pd.DataFrame(
        columns=["ID", "Best_Fit", "Best_Match","Best_Match_String", "Similarity", "Inverse"])
    print2(reference)

    df_path = os.path.join(root.contacts, reference.name+".csv")

    contacts_df = pd.read_csv(df_path)
    print(contacts_df)
    if reference.name != "GR":
        return
    assert reference.name == "GR"
    if "GR_EVA.csv" in os.listdir(root.data):
        eva_df = pd.read_csv(os.path.join(root.data, "GR_EVA.csv"))
        print(eva_df)
    else:
        print("GR_EVA.csv not found (in data folder)")
        return


    n_dimers = len(contacts_df.columns) - 2
    print2("Number of dimers:", n_dimers)

    groups = []
    n_groups = (len(eva_df.columns) - 1) / 2
    assert n_groups % 2 == 0
    ref_nums = eva_df["ResNum"]

    progress = ProgressBar(n_dimers)
    from maths import difference_between_boolean_pairs
    for c in range(n_dimers):
        dimer_sasas = contacts_df.iloc[:, [0, c + 2]]
        dimer_id = contacts_df.columns[c + 2]
        similarities = []
        for group in range(int(n_groups)):
            ref_sasaA = eva_df.iloc[:, [0, group * 2 + 1]]
            ref_sasaB = eva_df.iloc[:, [0, group * 2 + 2]]
            total_len = len(ref_sasaA)

            diffX = [0, 0]
            diffx = [0, 0]

            for ref_num in ref_nums:
                print(group, ref_num, end="\r")
                # print(dimer_sasas.loc[dimer_sasas["ResNum"] == ref_num].values)
                sA, sB = clean_list([dimer_sasas.loc[dimer_sasas["ResNum"] == ref_num].values[0, 1]],
                                    delimiter=",", format="bool")
                rsA = ref_sasaA.loc[ref_sasaA["ResNum"] == ref_num].iloc[:, 1].values[0]
                rsB = ref_sasaB.loc[ref_sasaB["ResNum"] == ref_num].iloc[:, 1].values[0]

                resX, resx = difference_between_boolean_pairs(sA, sB, rsA, rsB)
                diffX[0] += resX[0]
                diffX[1] += resX[1]
                diffx[0] += resx[0]
                diffx[1] += resx[1]

            if diffX[0] != 0:
                diffX = diffX[0] / diffX[1]
            else:
                diffX = 0
            if diffx[0] != 0:
                diffx = diffx[0] / diffx[1]
            else:
                diffx = 0

            inverse = False
            if diffx > diffX:
                inverse = True
            similarities.append((group + 1, max([diffX, diffx]), inverse))
            print1(similarities[-1][0], round(similarities[-1][1], 2))
        best_match = max(similarities, key=lambda x: x[1])
        from faces import GR_groups
        best_match_string = GR_groups[best_match[0]]
        print("Best match for {}: {} ({}), with {}% similarity, inverse: {}\n".format(dimer_id, best_match[0],
                                                                                 best_match_string,
                                                                                 round(100 * best_match[1]),
                                                                                 best_match[2]))

        vars.clustering["classified"][reference.name].loc[dimer_id] = [dimer_id, reference.name, best_match[0],
                                                                       best_match_string,
                                                 round(best_match[1] * 100), best_match[2]]
        progress.add()
    #classified_path = os.path.join(root.dataframes, "classified_df.csv")
    #vars.classified_df.to_csv(classified_path)
    #return classified_path
    vars.clustering["classified"][reference.name].to_csv(classified_path)
    return vars.clustering["classified"][reference.name]

def add_info_to_classified(reference):
    #classified = pd.read_csv(os.path.join(root.classified, "{}.csv".format(reference.name)), index_col=0)
    classified = reference.classified_df
    print(classified)
    faces = pd.read_csv(os.path.join(root.faces, "{}.csv".format(reference.name)), index_col=0)
    try:
        assert len(classified) == len(faces)
    except AssertionError:
        print(classified),
        print(faces)
        print("Classified and faces dataframes do not have the same length:")
        print(len(classified), len(faces))
        quit()
    if "ID" in classified.columns:
        classified.sort_values("ID", ascending=True, inplace=True)
    else:
        classified.sort_values(classified.columns[0], ascending=True, inplace=True)
    print(classified)
    faces.sort_values("ID", ascending=True, inplace=True)
    print(classified)
    print(faces)

    print(list(classified.columns))
    original_cols = list(classified.columns)
    print(original_cols)
    from copy import deepcopy
    for c in deepcopy(original_cols):
        #print(c)
        #print(original_cols)
        if "face1" in c or "face2" in c or "Unnamed" in c:
            #print("removing:", c)
            original_cols.remove(c)
            #print(original_cols)

    print(original_cols)
    print([classified[original_cols], faces[["face1", "face2", "contact_face1", "contact_face2"]]])

    vars.clustering["classified"][reference.name] = classified[original_cols].join(faces[["face1", "face2", "contact_face1", "contact_face2"]])
    #reference.classified_df = vars.clustering["classified"][reference.name]
    reference.faces_df = vars.clustering["faces"][reference.name]
    vars.clustering["classified"][reference.name].to_csv(os.path.join(root.classified, reference.name + ".csv"))

def add_clusters_to_classified(reference, pca=True, splitted = True):

    classified_path = os.path.join(root.classified, reference.name + ".csv")
    classified = pd.read_csv(classified_path, index_col="ID")
    if "ID" in classified.columns:
        classified.sort_values("ID", ascending=True, inplace=True)
    else:
        classified.sort_values(classified.columns[0], ascending=True, inplace=True)
    print(classified)
    #empty_face = pd.DataFrame([None]*len(classified), columns=["face_group"])
    #empty_cluster = pd.DataFrame([None]*len(classified), columns=["cluster"])
    #classified = pd.concat([classified, empty_face, empty_cluster ], axis=1, ignore_index=True)
    original_cols = list(classified.columns)
    print(original_cols)
    from copy import deepcopy
    for c in deepcopy(original_cols):
        # print(c)
        # print(original_cols)
        if "Unnamed" in c:
            # print("removing:", c)
            original_cols.remove(c)
    classified = classified[original_cols]
    print("########################################")
    print(classified)

    splitted_in_path = "clustered_pcas_{}".format(reference.name) in root.list()

    for splitted in list(set([splitted_in_path, False])):
        print("##########", splitted)
        if pca:
            if splitted:
                clustered_folder = root["clustered_pcas_{}".format(reference.name)]
            else:
                clustered_folder = root["clustered_pcas"]
        else:
            clustered_folder = root["clustered_{}".format(reference.name)]
        print(classified)
        for path in os.listdir(clustered_folder):
            if "centres" not in path and "clustered" not in path:
                df = pd.read_csv(os.path.join(clustered_folder, path))
                #print(df)
                for row in df.itertuples():
                    if row.id in classified.index:
                        if splitted:
                            classified.loc[row.id, "face_group"]= path.split(".")[0]
                            classified.loc[row.id, "cluster"]= row.cluster
                        else:
                            classified.loc[row.id, "global_cluster"]= row.global_cluster
        print(classified)
    vars.clustering["classified"][reference.name] = classified
    reference.classified_df = classified
    classified.to_csv(classified_path)




def split_by_faces(reference, force= False, by_com = True):
    if "face_contacts" in reference.__dict__.keys() and not force:
        return
    faces_dict = {}
    for row in reference.faces_df.itertuples():
        if by_com:
            face1 = row.face1
            face2 = row.face2
        else:
            face1 = row.contact_face1
            face2 = row.contact_face2

        if face1 is None or face2 is None:
            continue
        faces = "_".join(sorted([face1, face2]))
        if faces not in faces_dict.keys():
            faces_dict[faces] = [row.ID]
        else:
            faces_dict[faces].append(row.ID)
        #print(faces)
        #print(row)

    reference.face_contacts = {}
    contacts_df = reference.contacts_df
    print(contacts_df)
    root["contacts_{}".format(reference.name)] = "dataframes/clustering/contacts/contacts_{}".format(reference.name)
    for face, names in faces_dict.items():
        face_df = pd.concat([contacts_df.iloc[:, 0:2], contacts_df[names]] ,axis = 1)
        face_df.to_csv(os.path.join(root["contacts_{}".format(reference.name)], face + ".csv"), index=False)
        reference.face_contacts[face] = face_df
        print(face_df)
    #print(pd.DataFrame(reference.faces_dfs))


def clusterize_pcas(subfolder, in_path,method="KMeans", quantile=0.1, n_sample_multiplier=0.5, force = False, n_clusters = None , dimensions = [0,1,2], splitted = True, bandwidth=0.1):
    print1("Clusterizing PCAs")
    if n_clusters is None:
        if splitted:
            from faces import GR_groups
            n_clusters = 0
            faces = os.path.basename(in_path).split(".")[0].split("_")[-2:]
            for value in GR_groups.values():
                print(faces, value)
                if faces[0] == value[0] and faces[1] == value[1]:
                    n_clusters += 1
        else:
            n_clusters = 20


    root["clustered_pcas"] = "dataframes/clustering/clustered_pcas"
    subfolder = subfolder.format("clustered_pcas")
    name = os.path.basename(in_path).split(".")[0]
    if splitted:
        os.makedirs(os.path.join(root.clustered_pcas, subfolder), exist_ok=True)
        clustered_path = os.path.join(root.clustered_pcas, subfolder, name + ".csv")
        centres_path = os.path.join(root.clustered_pcas, subfolder, name + "_centres.csv")
        if name + ".csv" in os.listdir(os.path.join(root.clustered_pcas, subfolder)) and not force:
            print2("Skipping clustering for {}".format(name))
            return clustered_path
    else:
        clustered_path = os.path.join(root.clustered_pcas, name + ".csv")
        centres_path = os.path.join(root.clustered_pcas, name + "_centres.csv")


    pca_df = pd.read_csv(in_path)
    print(pca_df)

    columns = ["variance_{}".format(d) for d in dimensions]
    if method == "KMeans":
        from sklearn.cluster import KMeans
        model = KMeans(n_clusters=n_clusters, random_state=6, algorithm="elkan")
        model.fit(pca_df.loc[:, columns])


    elif method == "MeanShift":
        X = pca_df.loc[:, columns].values
        print(X)
        from sklearn.cluster import MeanShift, estimate_bandwidth
        if bandwidth is None:
            if n_sample_multiplier is not None:
                bandwidth = estimate_bandwidth(X, quantile=quantile, n_samples=int(round(len(X)*n_sample_multiplier)))
            else:
                bandwidth = estimate_bandwidth(X, quantile=quantile)
        #TODO: define it at main
        bandwidth = bandwidth
        model = MeanShift(bandwidth=bandwidth, bin_seeding=False, cluster_all=False, n_jobs=-1)
        model.fit(X)
        labels = model.labels_
        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)

        sprint("Number of estimated clusters : {}, bandwith: {}".format(n_clusters_, bandwidth))


    elif method == "SpectralClustering":
        from sklearn.cluster import SpectralClustering
        pass

    else:
        print("Please use a valid clustering algorithm")
        print("Introduced value:", method)
        quit()





    for n in range(len(pca_df)):
        cluster = model.labels_[n]
        if splitted:
            pca_df.loc[n, "cluster"] = cluster
        else:
            pca_df.loc[n, "global_cluster"] = cluster
        new_colour = "".join(["C", str(cluster)])

        if new_colour == "C":
            new_colour = "black"
        pca_df.loc[n, "colour"] = new_colour

    cluster_centres = pd.DataFrame(model.cluster_centers_)
    print(cluster_centres)
    print(pca_df)

    print(centres_path)
    print(clustered_path)
    cluster_centres.to_csv(centres_path, index=False)
    pca_df.to_csv(clustered_path)
    return clustered_path


def plot_clustered_pcas(reference, force=True, dimensions = 3, pca_dimensions = (0,1,2), labels = False,
                        labels_centres=False, adjust=False, subfolder=None, in_path=None, use_csv = True,
                        plot_centres=True, subset = None, pca=False, splitted=True):
    print1("Plotting: {}, Dimensions: {}".format(reference.name, dimensions))

    columns = ["variance_{}".format(d) for d in pca_dimensions]

    if dimensions > 3:
        print2("Plotting more than 3 dimensions not currently supported")
        return None

    root["pca_figs"] = "images/pca_figs"
    fig_folder = root.pca_figs
    fig_folder_name = "pca_figs"
    clustered_folder = root.clustered_pcas
    name = os.path.basename(in_path).split(".")[0]
    if splitted:

        clustered_subfolder = subfolder.format("clustered_pcas")
        subfolder = subfolder.format(fig_folder_name)
        os.makedirs(os.path.join(fig_folder, subfolder), exist_ok=True)

        fig_path = os.path.join(fig_folder, subfolder, name + ".png")
        if name + ".png" in os.listdir(os.path.join(fig_folder, subfolder)) and not force:
            print2("Skipping plotting for {}".format(name))
            return fig_path

    else:
        fig_path = os.path.join(fig_folder, name + ".png")
        if name + ".png" in os.listdir(fig_folder) and not force:
            print2("Skipping plotting for {}".format(name))
            return fig_path

    pca_df = pd.read_csv(in_path)
    print(pca_df)

    if not splitted:
        cluster_column = "global_cluster"
    else:
        cluster_column = "cluster"
    import matplotlib.pyplot as plt
    title = "PCA CLUSTERING for {}: {}  (N = {})".format(reference.name,name, len(pca_df))
    fig = plt.figure()
    if len(pca_dimensions) == 3:
        ax = fig.add_subplot(111, projection='3d')
    elif len(pca_dimensions) <= 2:
        ax = fig.add_subplot(111)
    ax.set_title(title)
    # ax.scatter(0,0,0, marker= "o", c="red")
    for row in pca_df.itertuples():
        #print(row.variance_0, row.variance_1, row.variance_2)
        if row.__getattribute__(cluster_column) == -1:
            ax.scatter(*[row.__getattribute__(v) for v in columns], c="black")
        else:
            ax.scatter(*[row.__getattribute__(v) for v in columns], c=row.colour)

    texts = []
    if plot_centres and subset is None:
        print2("Plotting cluster centres")
        if subfolder is None:
            cluster_centres = reference.cluster_centres
        else:
            if splitted:
                cluster_centres = pd.read_csv(os.path.join(clustered_folder,clustered_subfolder, "{}_centres.csv".format(name)))
            else:
                cluster_centres = pd.read_csv(
                    os.path.join(clustered_folder, "{}_centres.csv".format(name)))

        from maths import get_closest_point, points_to_line
        #print(cluster_centres)
        #print(cluster_centres.values)
        #print(cluster_centres.columns)
        for centre in cluster_centres.itertuples():
            #print(centre)
            if centre[0] != -1:
                ax.scatter(*[centre.__getattribute__("_"+str(x+1)) for x in range(len(pca_dimensions))], color="black", marker=".")
                if labels_centres:
                    print3("annotating", centre[0], end="\r")#, centre[1:])
                    ax.text( *[c+.001 for c in centre[1:]], centre[0], color="C{}".format(centre[0]))
                    pass

        centres = []
        for centre in cluster_centres.itertuples():
            #print(centre[0])
            #print(centre[1:])
            #print(centre)
            #print(*[x+1 for x in pca_dimensions])
            #centres.append([centre.__getattribute__("_"+str(x+1)) for x in pca_dimensions])
            centres.append(centre[1:])
            # print(centre[0],centre[2],centre[3])

        lines = []
        #print(centres)
        n = 0

        for row in pca_df.itertuples():
            point = [row.__getattribute__(v) for v in columns]
            #print("point:", point)
            #print("centres:", centres)
            closest = get_closest_point(point, centres)
            # print("closest:",closest)
            #print(closest)
            #print(point)
            lines.append(points_to_line(closest, point))
            # print("point to line:", points_to_line(closest,point))
            # lines[str(n)] = line

        # print("lines")
        # print(lines)
        for line in lines:
            # print(line)
            ax.plot(*line, c="black")
        # scripts.Mpl.plot_lines(lines, ax)
        # ax.plot(data=lines(lines[2],lines[1]))



    fig.tight_layout()
    ax.set_aspect('equal')
    try:
        ax.view_init(elev=45, azim=45)
    except:
        pass
    print2("Saving at {}".format(fig_path))
    print(fig_path)
    fig.savefig(fig_path, dpi=300)
    plt.show(block=vars.block)
    return fig_path


def quick_cluster(coords, n_clusters=3, method ="MeanShift",bandwidth = None):

    if method=="KMeans":
        from sklearn.cluster import KMeans
        model = KMeans(n_clusters=n_clusters, random_state=6, algorithm="elkan")
    elif method=="MeanShift":
        from sklearn.cluster import MeanShift, estimate_bandwidth
        print1("Clustering with MeanShift, bandwith:", end=" ")

        if bandwidth is None:
            bandwidth = estimate_bandwidth(coords)
        print2(bandwidth)
        model = MeanShift(bandwidth=bandwidth, bin_seeding=False, cluster_all=False, n_jobs=-1)


    model.fit(coords)
    print2("Clusters:", [int(l) for l in set(model.labels_)])
    return model.labels_



def remove_redundancy(in_path, threshold=0.001):
    pca_df = pd.read_csv(in_path)
    pca_df.sort_values(by=["variance_0", "variance_1", "variance_2"], inplace=True)
    print(pca_df)

    pca_df.reset_index( inplace=True)
    print(pca_df)

    dimer_list = list(pca_df["id"].values)
    starting_len = len(dimer_list)
    print(dimer_list)
    print("Starting_len:", starting_len)

    for row1 in pca_df.itertuples():
        #print1(row1)
        sprint(row1.id, row1.Index)
        for row2 in pca_df.iloc[row1.Index+1:,:].itertuples():
            print2(row2.Index, end="\r")
            if abs(row1.variance_0-row2.variance_0) < threshold:
                if abs(row1.variance_1-row2.variance_1) < threshold:
                    if abs(row1.variance_2 - row2.variance_2) < threshold:
                        print3("Removing:", row1.id, row1.Index)
                        dimer_list.remove(row1.id)
                        break
    print(dimer_list)

    print("Removed {} redundant dimers".format(starting_len - len(dimer_list)))
    pca_df = pca_df.loc[pca_df["id"].isin(dimer_list)]
    print(pca_df)
    pca_df.to_csv(in_path)
    print("Saved at:", in_path)
    return in_path



def cluster_by_face(reference, FORCE_ALL=False, DIMENSIONS=3, n_clusters = 4, minimum_score=0, pca = True,
                    pca_dimensions = [0,1,2], splitted=True, rem_red = True, method = "KMeans", quantile=0.1, n_sample_multiplier = 0.5,
                    bandwidth=0.1):



    if FORCE_ALL:
        FORCE_SM = True
        FORCE_CC = True
        FORCE_CLUSTER = True
        FORCE_PLOT = True
    else:
        FORCE_SM = False
        FORCE_CC = False
        FORCE_CLUSTER = False
        FORCE_PLOT = True

    if splitted:
        subfolder_name = "{}_" + reference.name
    else:
        subfolder_name = "{}"

    #print(root[subfolder_name.format("contacts")])

    for file in os.listdir(root[subfolder_name.format("contacts")]):

        if "filtered" in file:
            continue
        if not splitted:
            if "GR" not in file:
                continue
            if "contacts" in file:
                continue

        sprint(file)
        contacts_path = os.path.join(root[subfolder_name.format("contacts")], file)
        #print(contacts_path)
        if minimum_score > 0:
            contacts_df = pd.read_csv(contacts_path)
            original_len = len(contacts_df)
            classified_df = pd.read_csv(os.path.join(root.classified, "{}.csv".format(reference.name)))
            print(classified_df)
            cols = list(contacts_df.columns)
            for col in cols:
                if col in ["ResNum", "ResName"]:
                    continue
                #print(col)
                #print(classified_df[classified_df["ID"] == col])
                class_row = classified_df[classified_df["ID"] == col]
                #print(class_row)
                if len(class_row) < 1 :
                    cols.remove(col)
                    continue
                class_row = class_row.iloc[0]
                #print(class_row.Similarity)
                if class_row.Similarity < minimum_score:
                    cols.remove(col)
            filtered_len = len(cols)
            #print(contacts_df[cols])
            print1("Filtered {} dimers with a threshold of {}".format(original_len-filtered_len, minimum_score))
            contacts_path = os.path.join(root[subfolder_name.format("contacts")], "filtered_{}_".format(minimum_score) + file)
            contacts_df[cols].to_csv(contacts_path)



        if not pca:
            sms_path= generate_sm(reference, force=FORCE_SM, subfolder = subfolder_name, in_path = contacts_path)
            ccs_path = cc_analysis(reference, force=FORCE_CC, dimensions=DIMENSIONS, subfolder = subfolder_name, in_path = sms_path)
            if ccs_path is None:
                print("CC analysis failed")
                return
            clustered_path = clusterize_cc(reference, force=FORCE_CLUSTER, dimensions=DIMENSIONS, n_clusters=n_clusters,
                                           subfolder = subfolder_name, in_path = ccs_path)
            plot_path = plot_cc(reference, labels=False, labels_centres=True, # force = FORCE_PLOT
                                          dimensions=DIMENSIONS, subfolder = subfolder_name, in_path = clustered_path)
        else:
            from faces import get_pca_df, plot_pcas
            pca_path = get_pca_df(in_path=contacts_path, subfolder=subfolder_name, force=FORCE_CC, splitted=splitted)
            if rem_red:
                pca_path = remove_redundancy(pca_path)
            #plot_pcas(pcas, title="GR: {}  (N = {})".format(file.split(".")[0], len(pcas)))
            clustered_path = clusterize_pcas(method=method, subfolder=subfolder_name, in_path = pca_path, force=FORCE_CLUSTER,
                                             dimensions=pca_dimensions,splitted=splitted, quantile =quantile, n_sample_multiplier=n_sample_multiplier, bandwidth=bandwidth)
            plot_path = plot_clustered_pcas(reference, labels=False, labels_centres=True,  force = FORCE_PLOT,
                                dimensions=DIMENSIONS, subfolder=subfolder_name, in_path=clustered_path, pca = True,
                                            pca_dimensions=pca_dimensions, splitted=splitted)
    add_clusters_to_classified(reference, pca = pca, splitted=splitted)
    return


def generate_dihedrals_df(dimer_list = None, force = False):
    sprint("Generating dihedrals dataframe")
    root["clustering2"] = "dataframes/clustering2"
    root["dihedrals"] = "dataframes/clustering2/dihedrals"
    if force or len(os.listdir(root.dihedrals)) == 0:
        from imports import load_single_pdb
        if dimer_list is None:
            dimer_list = os.listdir(local.dimers)
        dimer_list = sorted(dimer_list)
        progress = ProgressBar(len(dimer_list))
        dataframes = {}
        for d in dimer_list:
            dimers = load_single_pdb(d, pickle_folder=local.dimers, quiet=True)
            for dimer in dimers:

                if dimer.best_fit is None or dimer.best_fit =="Missmatch":
                    progress.add(info=dimer.id)
                    continue
                chains = [dimer.monomer1.chain, dimer.monomer2.chain]
                dihedrals_1to2 = dimer.get_dihedrals(reverse=False)
                dihedrals_2to1 = dimer.get_dihedrals(reverse=True)

                if dimer.best_fit in dataframes.keys():
                    dataframes[dimer.best_fit].append([dimer.id, True] + chains + dihedrals_1to2 + dihedrals_2to1[3:-1])
                    dataframes[dimer.best_fit].append([dimer.id, False] + chains + dihedrals_2to1 + dihedrals_1to2[3:-1])
                else:
                    dataframes[dimer.best_fit] = [[dimer.id, True] + chains + dihedrals_1to2 + dihedrals_2to1[3:-1]]
                    dataframes[dimer.best_fit] = [[dimer.id, False] + chains + dihedrals_2to1 + dihedrals_1to2[3:-1]]

                progress.add(info=dimer.id)
        for key in dataframes.keys():
            print(key)
            print(dataframes[key])
            df = pd.DataFrame(columns = ["id", "is1to2", "mon1", "mon2", "d0", "d1", "d2", "a0", "a1", "a2", "d", "b0", "b1", "b2",], data = dataframes[key])
            print(df)
            df_path = os.path.join(root.dihedrals, key + ".csv")
            df.to_csv(df_path)


def plot_dihedrals(path, clusters=None, ax_labels=["0","1","2"], subset_col = None, subset=None, include_all=True, save=True,
                   label_col=None, only_first=None, heatmap=False, hm_threshold = 10, outer_ids_complete = None,
                   gif = False, snapshot=False, first_matrix_only = True, chainbows = False, get_matrix = False):
    print1("plotting dihedrals, heatmap={}, GIF={}, snapshot={}".format(heatmap, gif, snapshot))
    print2(path)
    from matplotlib import  pyplot as plt
    from imports import load_single_pdb
    complete_df = pd.read_csv(path)

    name = os.path.basename(path).split(".")[0]

    r = []

    if subset_col is not None:
        assert subset_col in complete_df.columns
        subsets = subset
        if subsets is None:
            subsets = list(set(complete_df[subset_col].values))
        if type(subsets) is not list:
            subsets = [subsets]
        if include_all:
            subsets = ["all"] + subsets
    else:
        subsets = ["all"]

    print2("Subsets:", subsets)
    for subset in subsets:
        print3(subset)
        if subset == "all":
            df = complete_df
        else:
            df = complete_df[complete_df[subset_col] == subset]
        hm = None

        if only_first is not None:
            df = df.iloc[:only_first]
        print(df)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        print4("Dihedral plotting, heatmap = {}".format(heatmap))
        progress = ProgressBar(len(df), silent=True)
        for point in df.itertuples():
            if clusters is None:
                ax.scatter(point.a0, point.a1, point.a2)
            else:
                cl = point.__getattribute__(clusters)
                if cl == -1:
                    col = "black"
                else:
                    col = "C"+str(cl)
                ax.scatter(point.a0, point.a1, point.a2, c=col)
                if heatmap or get_matrix:
                    dimer = load_single_pdb(point.id, pickle_folder=local.dimers, first_only=True, quiet=True)
                    if hm is None:
                        hm = dimer.contact_surface.get_contact_map(threshold=hm_threshold, transposed=not point.is1to2)
                    else:
                        hm = np.add(hm, dimer.contact_surface.get_contact_map(threshold=hm_threshold, transposed=not point.is1to2))
            if label_col is not None:
                ax.text(point.a0, point.a1, point.a2, point.__getattribute__(label_col))
            progress.add(info=point.id)
        ax.set_xlabel(ax_labels[0])
        ax.set_ylabel(ax_labels[1])
        ax.set_zlabel(ax_labels[2])
        ax.set_xlim(0,180)
        ax.set_ylim(0,180)
        ax.set_zlim(0,180)
        title = "{}-{}".format(name,subset)
        ax.set_title(title+" N={}".format(len(df)))

        if save:
            root["dihedral_figs"] = "images/dihedral_figs"
            savepath = os.path.join(root.dihedral_figs, title + ".png")
            plt.savefig(savepath)
        if gif:
            root["dihedral_gifs"] = "images/dihedral_gifs"
            mpl_to_gif(fig, ax, name=title, folder= root.dihedral_gifs)
        if heatmap:
            root["heatmap_figs"] = "images/heatmap_figs"
            hm_title = title + "_heatmap"
            from faces import ContactSurface
            matrix, oneDmatrix1, oneDmatrix2 = ContactSurface.get_heat_map(hm, title=hm_title, normalize=len(df), folder=root.heatmap_figs,
                                        percentage=True, outer_ids_complete=outer_ids_complete)
            r.append([matrix, oneDmatrix1, oneDmatrix2])
        else:
            r.append([None,None,None])
        if snapshot:
            if heatmap and False:
                cluster_snapshot(file=path,clusters=["all",subset], matrix1=oneDmatrix1, matrix2 = oneDmatrix2)
            else:
                cluster_snapshot(file=path,clusters=["all",subset], chainbows =chainbows )
        if vars.block:
            plt.show(block = vars.block)
        plt.close()

    if first_matrix_only:
        return r[0]
    else:
        return r





def cluster_snapshot(file, clusters, color_clusters=False, chainbows = False, snapshot = True, matrix1=None, matrix2=None):
    from imports import load_single_pdb, load_references
    from superpose import superpose_many_chains


    sprint("Cluster snapshot")

    cluster_folders = ["angle_clusters2"]
    cluster_cols = ["angle_cluster2"]

    options = []
    sele = []
    fname = file.split(".")[0]
    cname = fname.split("-")[-1]
    for n, (cluster_col, cluster_folder) in enumerate(zip(cluster_cols, cluster_folders)):
        dihedrals_path = os.path.join(root[cluster_folders[n]], fname + ".csv")
        df = pd.read_csv(file, index_col=0)
        # print(df)
        options.append([int(a) for a in sorted(set(df[cluster_cols[n]].values))])
        s = clusters[n]
        if s == "all":
            sele.append(options[n])
        else:
            sele.append([s])
        df.query(" | ".join(["{} == {}".format(cluster_col, n) for n in sele]), inplace=True)
        # print(df.to_string())
        fname = fname + "-{}".format(s)

    # pymol_start(show=False)
    print(file)
    filename = os.path.basename(fname)
    print(filename)
    ref = load_references(identifier=filename.split("-")[0])[0]
    print("Sele:", sele)
    resids = [res.id[1] for res in ref.structure.get_residues()]

    for c in sele[-1]:
        if c == -1 or c == "-1":
            continue
        print("Cluster:", sele[:-1], c)
        subset = df[df[cluster_cols[-1]] == c]
        print(subset)
        chains_to_align = {ref.name: (ref.path, ref.chain, True)}
        for row in subset.itertuples():
            dimer = load_single_pdb(identifier=row.id, pickle_folder=local.dimers, quiet=True)[0]
            name = row.id + str(row.is1to2)
            if row.is1to2:
                chains_to_align[name] = (dimer.replaced_path, row.mon1, True)
            else:
                chains_to_align[name] = (dimer.replaced_path, row.mon2, False)
        print(chains_to_align)
        local["clusters"] = "clusters"
        subcname = "{}-CLUSTER-{}-{}".format(ref.name,cname,c)
        super_data = superpose_many_chains(chains_to_align, file_name=subcname+".pdb", save_folder=local.clusters)
        print(super_data)

        if snapshot:
            local["snapshots"] = "snapshots"
            from pyMol import pymol_start, pymol_load_path, pymol_colour, pymol_list_to_bfactors, pymol_save_snapshot, \
                mpl_colours, mpl_ncolours, pymol_save_temp_session, pymol_open_session_terminal, pymol_split_states, \
                pymol_orient, pymol_reinitialize, pymol_get_all_objects
            pymol_reinitialize()
            monster = pymol_load_path(super_data["out_path"], subcname)
            pymol_split_states(monster)
            pymol_orient()
            if chainbows:
                pymol_colour("chainbow", "(all)")
                pymol_save_snapshot(subcname + "chainbows", folder=local.snapshots)
            elif matrix1 is not None and matrix2 is not None:
                for obj, value in zip(pymol_get_all_objects(), chains_to_align.values()):
                    chain = value[1]
                    sele1 = obj + " and (c. {})".format(chain)
                    sele2 = obj + " and !(c. {})".format(chain)
                    pymol_list_to_bfactors(val_list=matrix1, obj_name=sele1, resids=resids)
                    pymol_list_to_bfactors(val_list=matrix2, obj_name=sele2, resids=resids)
                pymol_colour("blue_yellow_red", "(all)", spectrum="b")
                pymol_save_snapshot(subcname + "heat_map", folder=local.snapshots)
            else:
                pymol_colour(mpl_colours[c % mpl_ncolours], "(all)")
                pymol_save_snapshot(subcname+"cluster_cols", folder=local.snapshots)



            #session_path = pymol_save_temp_session()
            #pymol_open_session_terminal(session_path)







def cluster_snapshot_old(file, clusters, levels=None, color_clusters=False, chainbows = True):
    from imports import load_single_pdb, load_references
    from pyMol import pymol_start, pymol_load_path, pymol_colour,pymol_list_to_bfactors, pymol_align_chains, pymol_group, \
        pymol_open_saved_cluster, pymol_get_all_objects, pymol_save_temp_session, pymol_save_cluster, pymol_open_session_terminal, \
        colours,ncolours, pymol_reset, pymol_orient, pymol_save_snapshot, get_all_obj, pymol_disable, pymol_delete, \
        pymol_command_in_new_process, pymol_reinitialize
    sprint("Cluster snampshot")

    cluster_folders = ["angle_clusters2"]
    cluster_cols = ["angle_cluster2"]
    if levels is not None:
        l = levels -1
        cluster_folders = cluster_folders[l:]
        cluster_cols = cluster_cols[:l]


    options = []
    sele = []
    fname = file.split(".")[0]
    for n, (cluster_col, cluster_folder) in enumerate(zip(cluster_cols, cluster_folders)):
        dihedrals_path = os.path.join(root[cluster_folders[n]], fname + ".csv")
        df = pd.read_csv(file, index_col=0)
        #print(df)
        options.append([int(a) for a in sorted(set(df[cluster_cols[n]].values))])
        s = clusters[n]
        if s == "all":
            sele.append(options[n])
        else:
            sele.append([s])
        df.query(" | ".join(["{} == {}".format(cluster_col, n) for n in sele]), inplace=True)
        #print(df.to_string())
        fname = fname + "-{}".format(s)

    #pymol_start(show=False)
    print(file)
    filename = os.path.basename(fname)
    print(filename)
    ref = load_references(identifier=filename.split("-")[0])[0]
    pymol_load_path(ref.path, ref.name)
    pymol_colour("chainbow", ref.name)
    print("Sele:", sele)

    for c in sele[-1]:
        if c == -1 or c == "-1":
            continue
        if c != 1: #DeBUG
            continue
        print("##",get_all_obj())


        print("Cluster:", sele[:-1], c)
        subset = df[df[cluster_cols[-1]] == c]
        print(subset)
        chains_to_align = [[ref.name, ref.chain]]
        for row in subset.itertuples():
            dimer = load_single_pdb(identifier=row.id, pickle_folder=local.dimers, quiet=True)[0]
            name = pymol_load_path(dimer.replaced_path, row.id + str(row.is1to2))
            if row.is1to2:
                chains_to_align.append([name, row.mon1])
            else:
                chains_to_align.append([name, row.mon2])
            if chainbows:
                pymol_colour("chainbow", name)
            elif color_clusters:
                c = int(row.__getattribute__(cluster_cols[-1]))
                pymol_colour(colours[c % ncolours], name)
                pymol_colour("chainbow", name)
            else:
                resids = [res.id[1] for res in dimer.monomer1.replaced.get_residues()]
                sele1 = name + " and c. {}".format(dimer.monomer1.chain)
                sele2 = name + " and c. {}".format(dimer.monomer2.chain)
                list1 = [min(x) for x in dimer.contact_surface.d_s_matrix.T]
                list2 = [min(x) for x in dimer.contact_surface.d_s_matrix]
                pymol_list_to_bfactors(val_list=list1, obj_name=sele1, resids=resids)
                pymol_list_to_bfactors(val_list=list2, obj_name=sele2, resids=resids)
                pymol_colour("blue_yellow_red", name, spectrum="b")
        print(chains_to_align)
        print(get_all_obj())

        pymol_align_chains(chains_to_align)
        pymol_orient()
        local["snapshots"] = "snapshots"
        #session_path = pymol_save_temp_session()
        #pymol_open_session_terminal(session_path)
        pymol_save_snapshot(filename + "-{}".format(c), folder=local.snapshots)
        print("hi")
        for obj, _ in chains_to_align[1:]:
            pymol_delete(obj)
    pymol_delete()
    collect_garbage()
    pymol_reinitialize()
    collect_garbage()
    import time
    time.sleep(10)






def cluster_angles(dihedrals_path,
                   bandwidth = None,
                   angles=["a0", "a1", "a2"],
                   cluster_name = "angle_cluster",
                   folder="angle_clusters1",
                   split_by=None,
                   save_together = True):
    sprint("Clustering angles")
    print1("Name:", cluster_name)
    print1("Folder:", folder)
    print1("Angles:", angles)
    dihedrals_df = pd.read_csv(dihedrals_path)
    if split_by is not None:
        subsets=[]
        clusters = set(dihedrals_df[split_by].values)
        for cluster in clusters:
            print3("Cluster:", cluster)
            subset_df = dihedrals_df[dihedrals_df[split_by] == cluster]
            subset_df[cluster_name] = quick_cluster(subset_df[angles], bandwidth=bandwidth)
            root[folder] = "dataframes/clustering2/{}".format(folder)
            print(subset_df)
            if save_together:
                subsets.append(subset_df)
            else:
                subset_df.to_csv(os.path.join(root[folder], os.path.basename(dihedrals_path).split(".")[0]+"-"+str(cluster)+".csv"))
        if save_together:
            new_df = pd.concat(subsets, axis=0)
            new_df.to_csv(os.path.join(root[folder], os.path.basename(dihedrals_path).split(".")[0]+".csv"))
    else:
        dihedrals_df[cluster_name] = quick_cluster(dihedrals_df[angles], bandwidth=bandwidth)
        root[folder] = "dataframes/clustering2/{}".format(folder)
        print(dihedrals_df)
        dihedrals_df.to_csv(os.path.join(root[folder], os.path.basename(dihedrals_path).split(".")[0]+".csv"))
    return root[folder]


def create_clusters(df_path, ref, include_all=True, **kwargs):

    df = pd.read_csv(df_path, index_col=0)
    cluster_cols = ["angle_cluster1", "angle_cluster2"]
    cluster1list = list(set(df[cluster_cols[0]]))

    print(df)
    print(cluster1list)
    new_clusters = []
    if include_all:
        new_clusters.append(Cluster2(ref,df,cluster_cols, is_all=True, **kwargs))
    for c1 in cluster1list:
        subset = df[df[cluster_cols[0]] == c1]
        cluster2list = list(set(subset[cluster_cols[1]]))
        for c2 in cluster2list:
            c = Cluster2(ref,df,cluster_cols, c1=c1, c2=c2, **kwargs)
            if not c.outlier:
                new_clusters.append(c)
    print(new_clusters)


class Cluster2:
    pickle_extension = '.cluster'
    pickle_folder = "cluster_pickles"
    name = "ClusterObject"
    path = None

    def __init__(self,ref, df,cluster_cols, c1=None, c2=None, is_all=False, **kwargs):
        self.kwargs = kwargs
        self.ref_name = ref.name
        self.structure = ref.structure
        self.outer_ids_complete = ref.get_outer_res_list(complete_list=True)
        self.outer_ids_binary = ref.get_outer_res_list(complete_list=True, binary=True)
        self.cluster_cols = cluster_cols
        self.outlier = False
        self.merged = []
        self.redundant = False
        self.redundant_to = None
        if not is_all:
            assert c1 is not None and c2 is not None
            self.is_all = False
            self.c1, self.c2 = c1, c2
            self.cnums = [self.c1, self.c2]
            if self.c1 == -1 or self.c2 == -1:
                self.outlier = True
            print(self.c1, self.c2, self.cnums)
            print(df)
            print(cluster_cols)
            print(" & ".join(["{} == {}".format(self.cluster_cols[n], self.cnums[n]) for n in range(len(self.cnums))]))
            self.subset = df.query(" & ".join(["{} == {}".format(self.cluster_cols[n], self.cnums[n]) for n in range(len(self.cnums))]), inplace=False)
            self.id = "{}-{}-{}".format(self.ref_name, add_front_0(self.c1, digits=2), add_front_0(self.c2, digits=2))
        else:
            self.c1, self.c2 = "all", "all"
            self.is_all = True
            self.subset = df.copy()
            self.id = "{}-{}-{}".format(self.ref_name, self.c1, self.c2)


        self.subset["reversed"] = [False] * len(self.subset)
        self.ndimers = len(self.subset)
        print(self.id)
        print(self.subset)
        if len(self.subset) == 0:
            quit()
        self.subset.sort_values(by="id", inplace=True)


        self.matrix = None
        self.oneDmatrix1 = None
        self.oneDmatrix2 = None
        self.plot_path = None
        self.gif_path = None
        self.snapshot_path = None
        self.comA = None
        self.comB = None
        self.stdA = None
        self.stdB = None



        if not self.outlier:
            self.process_cluster(**self.kwargs)
            self.pickle()

    def __repr__(self):
        return "<Cluster:{} {}-{} /N={}>". format(self.ref_name, self.c1, self.c2, self.ndimers)

    def process_cluster(self, force = False, matrix=False, plot=False, gif=False, show=False, snapshot=False, **kwargs):
        print1("Processing matrix={}, plot={}, gif={}, show={}".format(matrix, plot, gif, show))
        if not self.is_all:
            self.remove_identical()
            if self.comA is None or force:
                self.get_com()

        if matrix:
            if self.matrix is None or force:
                self.get_matrix(threshold=10,)
                if plot:
                    self.show_mpl(gif=gif, show=show)
        if plot:
            if self.plot_path is None or force:
                self.plot_path, self.gif_path = self.get_plot(self.subset, self.cluster_cols, self.id,
                                                              coms=(self.comA, self.comB),
                                                              stds=(self.stdA, self.stdB),
                                                              show=show, gif=gif)
        if snapshot:
            if self.snapshot_path is None or force:
                self.snapshot_path = self.show(snapshot=True, show_session=False)







    def get_com(self):
        from maths import find_com, distance
        anglesA = np.array(self.subset[["a0", "a1", "a2"]].values)
        anglesB = np.array(self.subset[["b0", "b1", "b2"]].values)
        self.comA = find_com(anglesA)
        self.comB = find_com(anglesB)
        distancesA = [distance(point, self.comA) for point in anglesA]
        distancesB = [distance(point, self.comB) for point in anglesB]
        self.stdA = np.std(distancesA)
        self.stdB = np.std(distancesB)
        #print(self.stdA, self.comA)
        #print(self.stdB, self.comB)


    def plot_matrix(self, show=False, **kwargs):
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        cvals = (0,0.5,1)
        colors = ("blue", "yellow", "red")
        norm = plt.Normalize(min(cvals), max(cvals))
        tuples = list(zip(map(norm, cvals), colors))
        cmap = mpl.colors.LinearSegmentedColormap.from_list("colormap", tuples)

        outer_ids_matrix = np.array([self.outer_ids_binary] * len(self.outer_ids_binary))
        colors2 = ((0, 0, 0, 1), (0.0, 0.0, 0.0, 0.0))
        tuples2 = list(zip([0, 1], colors2))
        cmap2 = mpl.colors.LinearSegmentedColormap.from_list("black0s", tuples2)

        fig, axes = plt.subplots(2, 2,
                                 # sharex="col",
                                 gridspec_kw={'height_ratios': [4, 2], "width_ratios": [2, 4]},
                                 figsize=(12, 9.6))

        ax = axes[0, 1]
        axLeft = axes[0, 0]
        axBottom = axes[1, 1]
        fig.subplots_adjust(right=0.8)
        main_fig1 = ax.imshow(self.matrix.T, cmap=cmap)
        cbar = plt.colorbar(main_fig1, cax=fig.add_axes([0.85, 0.15, 0.05, 0.7]))
        ax.imshow(outer_ids_matrix, cmap=cmap2)
        ax.imshow(outer_ids_matrix.T, cmap=cmap2)
        #print(self.oneDmatrix1)
        #print(self.oneDmatrix2)
        #print(self.outer_ids_complete)
        for n, (p1b, p2b, ob) in enumerate(zip(self.oneDmatrix1, self.oneDmatrix2, self.outer_ids_complete)):
            if n == 0:
                continue
            n2 = - n
            p1a, p2a, oa = self.oneDmatrix1[n-1], self.oneDmatrix2[n-1], self.outer_ids_complete[n-1]
            if oa is None:
                c1a = "black"
                c2a = "black"
            else:
                c1a = cmap(p1a)
                c2a = cmap(p2a)

            if ob is None:
                c1b = "black"
                c2b = "black"
            else:
                c1b =  cmap(p1b)
                c2b =  cmap(p2b)

            middle1 = (p1a + p1b) / 2
            middle2 = (p2a + p2b) / 2


            axBottom.plot((n - 1, n-0.5), (p1a, middle1), c=c1a, linestyle='--', linewidth=0.5)
            axBottom.plot((n, n - 0.5), (p1b, middle1), c=c1b, linestyle='--', linewidth=0.5)
            axLeft.plot( (p2a, middle2), (n2 + 1, n2 + 0.5), c=c2a, linestyle='--', linewidth=0.5)
            axLeft.plot( (p2b, middle2), (n2, n2 + 0.5), c=c2b, linestyle='--', linewidth=0.5)


        cbar_title = "Occurence"
        cbar.set_label(cbar_title)
        #ax.set_title(self.id + " N={}".format(self.ndimers))
        fig.suptitle(self.id + " N={}".format(self.ndimers))


        #fig.tight_layout()
        local["heatmaps"] = "images/heatmaps"
        fig_path = os.path.join(local.heatmaps, self.id + ".png")
        plt.savefig(fig_path)
        if show:
            plt.show(block=vars.block)
        return fig_path




    def get_matrix(self, threshold, plot = True, **kwargs):
        from imports import load_single_pdb
        print2("Generating cluster matrix")
        matrix = None
        self.subset.sort_values(by="id", inplace=True)
        progress = ProgressBar(len(self.subset), silent=True)
        for point in self.subset.itertuples():
            dimer = load_single_pdb(point.id, pickle_folder=local.dimers, first_only=True, quiet=True)
            is1to2 = point.is1to2
            if point.reversed:
                is1to2 = not is1to2

            if matrix is None:
                matrix = dimer.contact_surface.get_contact_map(threshold=threshold, transposed=not is1to2)
            else:
                matrix = np.add(matrix, dimer.contact_surface.get_contact_map(threshold=threshold, transposed=not is1to2))
            progress.add(info=point.id)

        oneDmatrix1 = [sum(i)/ len(self.outer_ids_complete) for i in matrix]
        oneDmatrix2 = [sum(i) / len(self.outer_ids_complete) for i in matrix.T]

        self.matrix, self.oneDmatrix1, self.oneDmatrix2 = matrix, oneDmatrix1, oneDmatrix2
        if plot:
            self.plot_matrix(**kwargs)

        return self.matrix, self.oneDmatrix1, self.oneDmatrix2

    @staticmethod
    def get_plot(subset, cluster_cols, cluster_id, coms=(None,None),stds=(None,None), gif=True, id_labels=False, save=True, show=False, show_outliers=False):
        subset.sort_values(by="id", inplace=True)
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(20,10))
        ax1 = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122, projection='3d')
        axes = ax1, ax2
        print2("Plotting cluster angles")
        progress = ProgressBar(len(subset), silent=True)
        for point in subset.itertuples():

            cl1 = point.__getattribute__(cluster_cols[0])
            cl2 = point.__getattribute__(cluster_cols[1])
            cols = []
            for cl in (cl1, cl2):
                if cl == -1:
                    cols.append("black")
                else:
                    cols.append("C" + str(cl))
            if not show_outliers and "black" in cols:
                continue
            ax1.scatter(point.a0, point.a1, point.a2, c=cols[1], edgecolors=cols[0], s=50, linewidths=2)
            ax2.scatter(point.b0, point.b1, point.b2, c=cols[1], edgecolors=cols[0], s=50, linewidths=2)
            if id_labels:
                ax1.text(point.a0, point.a1, point.a2, point.id)
                ax2.text(point.a0, point.a1, point.a2, point.id)
            progress.add(info=point.id)
        if None not in coms:
            ax1.scatter(coms[0][0], coms[0][1], coms[0][2], c=cols[0], s=stds[0]*180, linewidths=2, alpha=0.3)
            ax2.scatter(coms[1][0],coms[1][1], coms[1][2], c=cols[1], s=stds[1]*180, linewidths=2, alpha=0.3)
        ax_labels = ["0", "1", "2"]
        title = "CLUSTER_{}".format(cluster_id)
        for ax, l in zip(axes, ("a", "b")):
            ax.set_xlabel(l+ax_labels[0])
            ax.set_ylabel(l+ax_labels[1])
            ax.set_zlabel(l+ax_labels[2])
            ax.set_xlim(0, 180)
            ax.set_ylim(0, 180)
            ax.set_zlim(0, 180)
            ax.set_title(l+" "+title + " N={}".format(len(subset)))

        fig_savepath = None
        gif_savepath = None
        if save:
            local["dihedral_figs"] = "images/dihedral_figs"
            fig_savepath = os.path.join(local.dihedral_figs, title + ".png")
            plt.savefig(fig_savepath)
        if gif:
            local["dihedral_gifs"] = "images/dihedral_gifs"
            gif_savepath = mpl_to_gif(fig, axes, name=title, folder=local.dihedral_gifs)
        if show:
            plt.show(block = vars.block)
        return fig_savepath, gif_savepath



    def pickle(self):
        import pickle
        local["pickles"] = "pickles"
        pickle_folder = os.path.join(local.pickles, self.pickle_folder)
        os.makedirs(pickle_folder, exist_ok=True)
        local[self.pickle_folder] = "pickles/{}".format(self.pickle_folder)
        file_name = "{}{}".format(self.id, self.pickle_extension)
        self.pickle_path = os.path.join(pickle_folder, file_name)
        with open(self.pickle_path, 'wb') as f:
            pickle.dump(self, f)


    def merge(self, cluster2):
        self.merged.append(cluster2.id)
        self.merged.extend(cluster2.merged)
        sub2 = cluster2.subset
        sub2.rename(columns={'a0': 'b0', 'a1': 'b1', 'a2': 'b2',
                             'b0': 'a0', 'b1': 'a1', 'b2': 'a2' }, inplace=True)
        for row in sub2.itertuples():
            sub2.loc[row.Index, "reversed"] = not row.reversed
        self.subset = pd.concat([self.subset, sub2], axis= 0)
        self.ndimers = len(self.subset)
        self.subset.sort_values(by="id", inplace=True)
        self.remove_identical()
        self.ndimers = len(self.subset)
        cluster2.redundant = True
        cluster2.redundant_to = self.id
        cluster2.pickle()
        self.pickle()


    def reprocess_cluster(self, **kwargs):
        self.process_cluster(**kwargs)
        self.pickle()


    def delete(self):
        try:
            os.remove(self.pickle_path)
        except:
            print("Failed to delete {}:".format(self.id))
            print(self.pickle_path)
            os.remove(self.pickle_path)

    def show(self, snapshot =True, show_session=False, chainbows=False, cluster_colours=False, show_snapshot=False, regenerate_matrix=False):
        if self.is_all and not show_session:
            return None

        from imports import load_references, load_single_pdb
        from superpose import superpose_many_chains
        from Bio.PDB import PDBParser, PDBIO, Structure
        ref = load_references(identifier=self.ref_name)[0]
        self.subset.sort_values(by="id", inplace=True)
        print(self.subset)
        chains_to_align = {ref.name: (ref.path, ref.chain, True, False, 0)}
        for n, row in enumerate(self.subset.itertuples()):
            dimer = load_single_pdb(identifier=row.id, pickle_folder=local.dimers, quiet=True)[0]
            name = row.id + str(row.is1to2)
            is1to2 = row.is1to2
            if row.reversed:
                is1to2 = not is1to2
            if is1to2:
                chains_to_align[name] = (dimer.replaced_path, row.mon1, row.is1to2, row.reversed, n+1)
            else:
                chains_to_align[name] = (dimer.replaced_path, row.mon2, row.is1to2, row.reversed, n+1)
        self.chains_to_align = chains_to_align
        local["cluster_pdbs"] = "exports/cluster_pdbs"

        def alter_bfactors(chain, value_list):
            assert  len(chain) == len(value_list)
            for atom, value in zip(chain.get_atoms(), value_list):
                atom.bfactor = value

        super_data = superpose_many_chains(chains_to_align, file_name=self.id + ".pdb", save_folder=local.cluster_pdbs)
        monster_path = super_data["out_path"]

        if regenerate_matrix:
            self.reprocess_cluster(matrix=True,force=True)
        if self.oneDmatrix1 is not None and self.oneDmatrix2 is not None:
            structure = PDBParser(QUIET=True).get_structure(self.id, monster_path)
            assert len(list(structure.get_models())) == len(chains_to_align)
            for model, (key,value) in zip(structure.get_models(), chains_to_align.items()):
                model.id = key
                for chain in model.get_chains():
                    is_chain1 = chain.id == value[1]
                    if is_chain1:
                        alter_bfactors(chain, self.oneDmatrix1)
                    else:
                        alter_bfactors(chain, self.oneDmatrix2)


            new_paths = []
            for n, model in enumerate(structure.get_models()):
                new_structure = Structure.Structure(model.id)
                new_structure.add(model)
                exporter = PDBIO()
                exporter.set_structure(new_structure)
                new_folder = os.path.join(local.cluster_pdbs, monster_path.split(".")[0])
                os.makedirs(new_folder, exist_ok=True)
                new_path = model.id + "_{}.pdb".format(add_front_0(n, digits=4))
                new_path = os.path.join(new_folder, new_path)
                exporter.save(new_path)
                new_paths.append(new_path)


        snapshot_path = None
        if snapshot:
            local["snapshots"] = "snapshots"
            from pyMol import pymol_start, pymol_load_path, pymol_colour, pymol_list_to_bfactors, pymol_save_snapshot, \
                mpl_colours, mpl_ncolours, pymol_save_temp_session, pymol_open_session_terminal, pymol_split_states, \
                pymol_orient, pymol_reinitialize, pymol_get_all_objects
            pymol_start(show=False)
            pymol_reinitialize()
            for file in new_paths:
                pymol_load_path(file, os.path.basename(file))
            pymol_orient()
            if chainbows:
                pymol_colour("chainbow", "(all)")
                extra_id = "_chainbows"
            elif cluster_colours:
                pymol_colour(mpl_colours[self.c2 % mpl_ncolours], "(all)")
                extra_id = "_cluster_cols"
            else:
                pymol_colour("blue_yellow_red", "(all)", spectrum="b")
                extra_id = "_heatmap"
            if not self.is_all:
                snapshot_path = pymol_save_snapshot(self.id + extra_id, folder=local.snapshots)
                if show_snapshot:
                    open_file_from_system(snapshot_path)
            if show_session:
                session_path = pymol_save_temp_session()
                pymol_open_session_terminal(session_path)

        return snapshot_path


    def remove_identical(self):
        id_list = []
        for row in self.subset.itertuples():
            if not row.id in [i[0] for i in id_list]:
                #print(id_list)
                #print(row.id)
                #print(row.id in id_list)
                id_list.append((row.id, row.Index))
            else:
                #print("identical dimers found:", row.id)
                #print(id_list)
                #print("Removing:")
                if row.is1to2:
                    #print(self.subset.loc[id_list[[i[0] for i in id_list].index(row.id)][1]])
                    self.subset.drop(id_list[[i[0] for i in id_list].index(row.id)][1], inplace=True)
                else:
                    #print(self.subset.loc[row.Index])
                    self.subset.drop(row.Index, inplace=True)
        self.ndimers = len(self.subset)

    def show_mpl(self, save=True, gif = False, show=False, mergedMatrix = None, secondary=None, title=None):
        if self.matrix is None:
            self.process_cluster(matrix=True)
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Line3D
        from maths import points_to_line, get_middlepoint, normalize1D
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        if mergedMatrix is None:
            mergedMatrix = [m1+m2 for m1, m2 in zip(self.oneDmatrix1, self.oneDmatrix2)]
            mergedMatrix = normalize1D(mergedMatrix)
            #cvals = (min(mergedMatrix), (max(mergedMatrix)-min(mergedMatrix))/2, max(mergedMatrix))
            cvals = (0,0.5,1)
            colors = ("blue", "yellow", "red")
            norm = plt.Normalize(min(cvals), max(cvals))
            tuples = list(zip(map(norm, cvals), colors))
            cmap = mpl.colors.LinearSegmentedColormap.from_list("colormap", tuples)
        else:
            cmap = None

        atom_list = list(self.structure.get_atoms())
        def get_colour(value, cmap=None):
            if cmap is not None:
                colour = cmap(value)
            else:
                if type(value) in [int, float, np.int64]:
                    if value == -1:
                        colour = "black"
                    else:
                        colour = "C" + str(value)
                else:
                    colour = value
            return colour
        for n, (atom, value) in enumerate(zip(atom_list, mergedMatrix)):
            colour0 = get_colour(mergedMatrix[n-1], cmap=cmap)
            colour1 = get_colour(value, cmap=cmap)
            size = 100
            if secondary is not None:
                size=secondary[n]*200*len(secondary)
            ax.scatter(atom.coord[0], atom.coord[1], atom.coord[2], s=size, c=colour1)
            if n != 0:
                middle = get_middlepoint(atom_list[n-1].coord, atom.coord)
                line1 = Line3D(*points_to_line(atom_list[n-1].coord, middle), linewidth=5, color=colour0)
                line2 = Line3D(*points_to_line(middle, atom.coord), linewidth=5, color=colour1 )

                ax.add_line(line1)
                ax.add_line(line2)
        ax_labels = ["X", "Y", "Z"]

        #cbar = plt.colorbar(ax, cax=fig.add_axes([0.85, 0.15, 0.05, 0.7]))
        #cbar_title = "Occurence"
        #cbar.set_label(cbar_title)
        ax.set_aspect("equal")
        ax.view_init(elev=-55., azim=-80, roll= 80)
        if title is None:
            title = self.id + " N={}".format(self.ndimers)
        fig.suptitle(title)
        fig_savepath = None
        gif_savepath = None
        if save:
            local["res_coords"] = "images/res_coords"
            fig_savepath = os.path.join(local.res_coords, title + ".png")
            plt.savefig(fig_savepath)
        if gif:
            local["res_coords_gifs"] = "images/res_coords_gifs"
            gif_savepath = mpl_to_gif(fig, ax, name=title, folder=local.res_coords_gifs)
        if show:
            plt.show(block=vars.block)
        return fig_savepath, gif_savepath







def cluster_redundancy(**kwargs):
    from imports import load_clusters
    from maths import distance
    done_clusters = []
    for cluster1 in load_clusters(onebyone=True):
        if cluster1.is_all or cluster1.redundant or cluster1.outlier:
            done_clusters.append(cluster1.id)
            continue
        print2(cluster1.id)
        for cluster2 in load_clusters(onebyone=False, identifier=cluster1.ref_name):

            if cluster1.id == cluster2.id:
                continue
            if cluster2.id in done_clusters or cluster2.is_all or cluster2.redundant or cluster2.outlier:
                continue
            #print3(cluster2.id)
            d1 = distance(cluster1.comA, cluster2.comB)
            d2 = distance(cluster1.comB, cluster2.comA)
            v1 = 10+ cluster1.stdA + cluster2.stdB
            v2 = 10+ cluster1.stdB + cluster2.stdA

            if d1 <= v1:
                if d2 <= v2:
                    print3("Redundant cluster:", cluster2.id)
                    print4(d1, v1)
                    print4(d2, v2)
                    cluster1.merge(cluster2)
        done_clusters.append(cluster1.id)

    for cluster in load_clusters(onebyone=True):
        if cluster.redundant:
            print(cluster.id, "is redundant to", cluster.redundant_to)
            cluster.delete()
        elif len(cluster.merged) != 0:
            cluster.reprocess_cluster(force=True, **kwargs)



def cluster_dihedrals():
    pass



def get_faces():
    from imports import load_clusters
    from maths import normalize1D
    for cluster in load_clusters(identifier="all", onebyone=True):
        if not cluster.is_all:
            continue


        coord_array = np.array([atom.coord for atom in cluster.structure.get_atoms()])
        preference_array = [m1 + m2 for m1, m2 in zip(cluster.oneDmatrix1, cluster.oneDmatrix2)]
        from sklearn.cluster import AffinityPropagation, KMeans, BisectingKMeans, SpectralClustering, FeatureAgglomeration
        print(coord_array)

        print("###")
        preference_array = normalize1D(preference_array, add_to=1)
        for n, p in enumerate(preference_array):
            print(add_front_0(n,digits=3, zero=" ")+ "|"+"#"*round(p*100)+" "*(100-round(p*100))+"|")
        print("###")
        #print(preference_array)


        #model = AffinityPropagation(random_state=6, damping=0.95).fit(coord_array)
        #m = np.amin(model.affinity_matrix_)
        #print("#####", m)
        #median = np.median(model.affinity_matrix_)
        #print("#####", median)
        #maximum = np.amin(model.affinity_matrix_)
        def custom_metric(coord1, coord2):
            weight1 = coord1[3]
            weight2 = coord2[3]
            coord1 = coord1[:3]
            coord2 = coord2[:3]
            from maths import distance
            print((1-(weight1+weight2)/2))
            return abs(distance(coord1, coord2)) * (1-(weight1+weight2)/2)
        #affinity = [v*maximum*len(preference_array)/5 for v in preference_array]
        #print(affinity)
        weighted_array = []
        for a, p in zip(coord_array, preference_array):
            print([*a]+[p])
            weighted_array.append([*a]+[p])
        print(weighted_array)
        model = AffinityPropagation(random_state=6, convergence_iter=100, verbose=True).fit(weighted_array)
        print(model.labels_)
        #model2 = FeatureAgglomeration(n_clusters=4).fit(model.cluster_centers_)
        #print(model2.labels_)
        from scipy.cluster.hierarchy import linkage, cut_tree
        Z = linkage(weighted_array, method='average', metric=custom_metric)
        #print(Z)
        #y = cut_tree(Z, 3)
        #print(y)
        from scipy.cluster.hierarchy import dendrogram

        D = dendrogram(Z, p=2)
        cluster.show_mpl(show=True, save=False, title = cluster.id+" n_clusters = {}".format(len(set(D["leaves_color_list"]))), mergedMatrix = D["leaves_color_list"], secondary=preference_array)

        faces = []
        for c in set(model.labels_):
            face=dict(
            com = model.cluster_centers_[c],
            C = c,
            N = sum([1 for i in model.labels_ if i ==c]),
            M = sum([preference_array[n] for n, _ in enumerate(model.labels_) if model.labels_[n] == c])
            )
            print("C:{}, N:{}, M:{}".format(face["C"], face["N"], face["M"] ))
            faces.append(face)
        faces = sorted(faces, key=lambda face: face["N"], reverse=True)
        [print(face) for face in faces]
        labels = []
        for l in model.labels_:
            if l in [d["C"] for d in faces]:
                labels.append(l)
            else:
                labels.append(-1)
        cluster.show_mpl(show=True, save=False, title = cluster.id+" n_clusters = {}".format(len(faces)), mergedMatrix = labels, secondary=preference_array)

