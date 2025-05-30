import os, sys

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


def get_clusters(df, column, ref_name):
    #df = pd.read_csv(df_path)
    clusters = []
    for c in set(df[column]):
        subset = df[df[column] == c]
        clusters.append(Cluster(ref_name, subset, c, column))
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
                   gif = False, snapshot=False, first_matrix_only = True):
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
                if heatmap:
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
        if snapshot:
            if heatmap and False:
                cluster_snapshot(file=path,clusters=["all",subset], matrix1=oneDmatrix1, matrix2 = oneDmatrix2)
            else:
                cluster_snapshot(file=path,clusters=["all",subset], chainbows =False )
        if vars.block:
            plt.show(block = vars.block)
        plt.close()
        r.append([None, None, None])
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
                   split_by=None):
    sprint("Clustering angles")
    print1("Name:", cluster_name)
    print1("Folder:", folder)
    print1("Angles:", angles)
    dihedrals_df = pd.read_csv(dihedrals_path)
    if split_by is not None:
        clusters = set(dihedrals_df[split_by].values)
        for cluster in clusters:
            print3("Cluster:", cluster)
            subset_df = dihedrals_df[dihedrals_df[split_by] == cluster]
            subset_df[cluster_name] = quick_cluster(subset_df[angles], bandwidth=bandwidth)
            root[folder] = "dataframes/clustering2/{}".format(folder)
            print(subset_df)
            subset_df.to_csv(os.path.join(root[folder], os.path.basename(dihedrals_path).split(".")[0]+"-"+str(cluster)+".csv"))
    else:
        dihedrals_df[cluster_name] = quick_cluster(dihedrals_df[angles], bandwidth=bandwidth)
        root[folder] = "dataframes/clustering2/{}".format(folder)
        print(dihedrals_df)
        dihedrals_df.to_csv(os.path.join(root[folder], os.path.basename(dihedrals_path).split(".")[0]+"-all.csv"))
    return root[folder]



