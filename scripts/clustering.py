import os
from itertools import count

from numpy.ma.extras import average

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


def get_clusters(df_path, column, ref_name):
    df = pd.read_csv(df_path)
    clusters = []
    for c in set(df[column]):
        subset = df[df[column] == c]
        clusters.append(Cluster(ref_name, subset, c, column))
    return clusters







### DEPRECATED ###
def generate_sm_old(reference, force=False):
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
##################

def generate_sm(reference, force=False):
    print1("Generating SM for {}".format(reference.name))

    contacts_df = pd.read_csv(os.path.join(root.contacts, reference.name+".csv"), index_col=0)
    print(contacts_df)

    sm_ssd = pd.DataFrame(columns=["dimer1", "dimer2", "index1", "index2", "similarity"])
    n_dimers = len(contacts_df.columns) - 2
    if n_dimers <2:
        print1("Not enough dimers in {} dataframe".format(reference.name))
        return
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
    root["sms"] = "dataframes/clustering/sms"
    sm_path = os.path.join(root.sms, '{}.csv'.format(reference.name))
    sm_ssd.to_csv(sm_path,header=False, index=False)
    return sm_path

### DEPRECATED ###
def cc_analysis_old(reference, dimensions=3, force =False):
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
    return cc_out_path
##################

def cc_analysis(reference, dimensions=3, force =False):
    print("CC analysis for {}".format(reference.name))

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
    return cc_out_path

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


def cluster(reference, FORCE_ALL=False, DIMENSIONS = 3, score_id = ""):
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
    return
    cc_analysis(reference, force=FORCE_CC, dimensions=DIMENSIONS)
    clusterize_cc(reference, force=FORCE_CLUSTER, dimensions=DIMENSIONS)
    plot_cc(reference, labels=False, labels_centres=True, force=FORCE_PLOT, dimensions=DIMENSIONS)

    if reference.name == "GR":
        tprint("Comparing to Eva")
        df_cc = pd.read_csv(os.path.join(root.clustered, "GR_cc.csv"))
        scores = calculate_scores_GR(df_cc, score_id+str(DIMENSIONS))
        print1("Scores: cc: {}, eva: {}".format(scores[0], scores[1]))
        eprint("Compared successfully")




#### DEPRECATED ###
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
###################

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



def compare_contacts(reference):
    print1("Comparing contacts for GR")

    print2(reference)
    df_path = os.path.join(root.contacts, reference.name+".csv")
    contacts_df = pd.read_csv(df_path)
    print(contacts_df)
    if reference.name != "GR":
        return
    assert reference.name == "GR"
    if "GR_EVA.csv" in os.listdir(root.contacts):
        eva_df = pd.read_csv(os.path.join(root.contacts, "GR_EVA.csv"))
        print(eva_df)
    else:
        print("GR_EVA.csv not found (in contacts folder)")
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
        print("Best match for {}: {}, with {}% similarity, inverse: {}\n".format(dimer_id, best_match[0],
                                                                                 round(100 * best_match[1]),
                                                                                 best_match[2]))

        vars.classified_df.loc[len(vars.classified_df)] = [dimer_id, reference.name, best_match[0],
                                                 round(best_match[1] * 100), best_match[2]]
        progress.add()
    classified_path = os.path.join(root.dataframes, "classified_df.csv")
    vars.classified_df.to_csv(classified_path)
    return classified_path





if __name__ == "__main__":




    import setup
    from Globals import root, local, vars
    from imports import load_references, load_single_pdb
    from dataframes import save_dfs


    #### CONTACT LENGTH TESTING ########################################################################################
    vars["references"] = load_references()

    molecule_folder = local.many_pdbs
    molecule_list = sorted(os.listdir(molecule_folder))
    # print(len(vars.do_only), vars.do_only)
    if len(vars.do_only) > 0:
        molecule_list = [f for f in molecule_list if any([s in f for s in vars.do_only])]
    # [print(m) for m in sorted(molecule_list)]
    print1("Molecule list obtained:", len(molecule_list), "molecules")


    print(list(vars.clustering["contacts"].keys()))
    progress = ProgressBar(len(molecule_list))
    from surface import build_contact_arrays
    for dist in [3,4,5,6,7,8]:
        for m in molecule_list:
            if "lock" in m:
                sprint(".lock file detected:", m)
                continue
            filename = m.split(".")[0]
            sprint(filename)
            molecules = load_single_pdb(filename, local.molecules)
            for molecule in molecules:
                dimers = molecule.dimers
                for dimer in dimers:
                    print1(dimer)
                    if dimer.incomplete:
                        continue
                    build_contact_arrays(dimer, sasa=False, force=True, max_contact_length=dist)
                    #dimer.pickle()
                progress.add(info=molecule.id)
        save_dfs()

        for reference in vars.references:
            cluster(reference, score_id="under_{}_A_".format(dist))
    ####################################################################################################################

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

    # OLD
    '''clustering(FORCE_ALL=False,
               FORCE_SM=False,
               FORCE_CC=True,
               FORCE_CLUSTER=True,
               FORCE_PLOT=True,
               DIMENSIONS=5,
               )'''
    #from github import automatic_push_to_branch
    #automatic_push_to_branch(target="auto")
            