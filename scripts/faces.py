import os, sys

from clustering import quick_cluster
from imports import load_single_pdb
from utilities import *
from Globals import root, local, vars
import numpy as np
import pandas as pd
from maths import *



import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
try:
    matplotlib.use('QtAgg')
except:
    try:
        matplotlib.use('TkAgg')
    except:
        matplotlib.use('Agg')

GR_dict = {
    "front":     [526, 527, 528, 529, 530, 531, 532,
                  534, 535,
                  538, 539, 540, 541, 542, 543,
                  545, 546, 547, 548, 549, 550, 551,
                  575, 576, 579, 580, 582, 583, 586,
                  593, 597, 598, 600,
                  614, 615, 616, 617, 620, 622, 625, 626, 628, 633],
    "base":      [552, 553, 554, 555, 556, 558,
                  565, 566, 568, 569, 571, 571,
                  746, 747,
                  # Mixed
                  734, 736, 737, 740, 741, 742, 743, 744, 745
                  ],
    "back":      [637, 638, 641, 644, 645, 647,
                  730, 731, 733, 734, 735, 736, 737, 738,
                  740, 741, 742, 743, 744, 745, 748,
                  750, 751, 752, 755, 759, 761, 762, 764, 765, 766, 768, 769],
    "top":       [677, 678, 683, 684, 685,
                  687, 688, 689, 690, 691, 692, 693, 694, 695, 698, 672,
                  709, 710, 711, 712, 713, 715, 716, 717,
                  770, 771, 772, 773, 774, 775, 776, 777],
}

GR_groups = {1: ("front", "front"),
             2: ("front", "front"),
             3: ("front", "back"),
             4: ("front", "back"),
             5: ("front", "base"),
             6: ("front", "front"),
             7: ("front", "top"),
             8: ("front", "top"),
             9: ("front", "front"),
             10: ("front", "front"),
             11: ("front", "front"),
             12: ("base", "base"),
             13: ("base", "top"),
             14: ("base", "back"),
             15: ("back", "back"),
             16: ("back", "back"),
             17: ("back", "back"),
             18: ("back", "top"),
             19: ("back", "top"),
             20: ("top", "top"),
             }
GR_groups = {key:sorted(value) for key, value in GR_groups.items()}

GR_colours = {"front": "purple",
              "base": "salmon",
              "back": "cyan",
              "top": "pink",}


def define_faces_from_list(self, list):
    for contact in self.contacts:
        print(contact)
        pass

def get_face_coms(monomer):
    assert monomer.best_fit == "GR"
    assert monomer.replaced is not None
    coms = {face: [] for face in GR_dict.keys()}
    #print(coms)
    for atom in monomer.replaced.get_atoms():
        resn = atom.parent.id[1]
        for face in GR_dict:
            if resn in GR_dict[face]:
                coms[face].append(atom.coord)

    for face in coms.keys():
        x, y, z = 0, 0, 0
        l = len(coms[face])
        for coord in coms[face]:
            x += coord[0]
            y += coord[1]
            z += coord[2]

        coms[face] = x/l, y/l, z/l
    print(coms)
    return coms


def get_dimer_faces(dimer):
    shortest = (None, None, 999)
    ## Development
    if "face_coms" not in dimer.monomer1.__dict__:
        dimer.monomer1.face_coms = get_face_coms(dimer.monomer1)
        dimer.monomer1.pickle()
    if "face_coms" not in dimer.monomer2.__dict__:
        dimer.monomer2.face_coms = get_face_coms(dimer.monomer2)
        dimer.monomer2.pickle()
    ##
    for face1, com1 in dimer.monomer1.face_coms.items():
        for face2, com2 in dimer.monomer2.face_coms.items():
            d = distance(com1, com2)
            if d < shortest[2]:
                shortest = (face1, face2, d)
    return shortest





def get_pca(structure, n_components = 3, com = None, closer_to = None, solver = "covariance_eigh", clockwise =True):
    print6("Getting PCA")
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_components, random_state=6, svd_solver=solver)
    coords = [atom.coord for atom in structure.get_atoms()]
    if com is None:
        com = find_com(coords)
    if closer_to == "N":
        closer_to = get_terminals(structure)["N"]
    elif closer_to == "C":
        closer_to = get_terminals(structure)["C"]
    #print("COM:", com)
    coords = [c-com for c in coords]
    pca.fit(coords)
    #print(pca)
    #print(pca.components_)
    #print(pca.explained_variance_ratio_)
    #print(pca.singular_values_)

    if closer_to is not None:
        for n, component in enumerate(pca.components_):
            closer = closer_to - com
            #print(distance(component, closer),distance(component * -1, closer))
            if distance(component, closer) > distance(component * -1, closer):
                if vars.verbose:
                    print4("Reversed component", n, pca.components_[n], "-->", end=" ")
                pca.components_[n] = component * -1
                #print(pca.components_[n])


    else:
        pca.inverse = False
        if clockwise:
            dh1 = dihedral_angle(pca.components_[1], com, pca.components_[0], add(pca.components_[2], com))
            dh2 = dihedral_angle(pca.components_[2], com, pca.components_[0], add(pca.components_[1], com))
            print(dh1, dh2)
            if dh1 < dh2:
                pca.inverse = True


        for n, component in enumerate(pca.components_):
            if vars.verbose:
                print4("Component {}: {} / Value: {} --> Vector: {} / Inverse: {}".format(n,
                                                                         pca.components_[n],
                                                                         pca.explained_variance_[n],
                                                                         pca.components_[n] * pca.explained_variance_[n],
                                                                         pca.inverse))




    return pca

def pca_to_lines(pca, com, just_points = False):
    components = pca.components_
    variances = pca.explained_variance_
    points = []
    lines = []
    #print("COmponents:", components)
    for component, variance in zip(components,variances):
        c = add(com,component)
        v = scale(vector(com, c), variance)
        sc = add(com, v)
        #print( "C:", c)
        points.append((com, sc))
        lines.append(points_to_line(com, sc))
        #print(points)
    if just_points:
        return points
    else:
        return lines

def get_terminals(structure):
    atom_list = list(structure.get_atoms())
    terminals = dict(N = atom_list[0].coord,
                     C = atom_list[-1].coord)
    return terminals



def plot_atoms(structure, pca = None, block = True):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(structure.id)

    marker2c = dict(marker='o', linestyle=':',# markersize=15,
                           color='darkgrey',
                           markerfacecolor='tab:blue',
                           markerfacecoloralt='lightsteelblue',
                           markeredgecolor='brown',
                           fillstyle='right')

    coords = [atom.coord for atom in structure.get_atoms()]
    com = find_com(coords)
    #print("COM:", com)
    for atom in structure.get_atoms():
        coord = atom.coord - com
        colours = []
        #print(atom.parent.id[1])
        for face in GR_dict.keys():
            if atom.parent.id[1] in GR_dict[face]:
                colours.append(GR_colours[face])
        #print(colours)
        if len(colours) == 0:
            ax.scatter(*coord, c="gray", marker= "o", s=50)
        else:
            for col in colours:
                ax.scatter(*coord, c=col, marker = "o", s=50)

    if pca is not None:
        original_components = pca.components_ * pca.explained_variance_#*pca.singular_values_
        pca_lines = []
        pca_lines.append(points_to_line((0, 0, 0), original_components[0]))
        pca_lines.append(points_to_line((0, 0, 0), original_components[1]))
        pca_lines.append(points_to_line((0, 0, 0), original_components[2]))
        # print(pca_lines)
        for line in pca_lines:
            # print("line:", line)
            # print("line items:", line[0][0], line[0][1], line[1][0], line[1][1], line[2][0], line[2][1])
            ax.plot(line[0], line[1], line[2], color='orange')


    fig.tight_layout()
    ax.set_aspect('equal')
    plt.show(block=block)
    return fig, ax




def plot_pcas(pca_list, title="", dimensions = [0,1,2], mode="variance", comps=[0,1,2], cluster=None,bandwidth=0.2, method="MeanShift"):
    fig = plt.figure()
    print("N dimensions:", len(dimensions))
    if len(dimensions) == 3:
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.add_subplot(111)
    ax.set_title(title)
    #ax.scatter(0,0,0, marker= "o", c="red")
    points = []
    labels = []
    for pca in pca_list:
        if mode == "variance":
            #print(pca.explained_variance_)
            #print(pca.explained_variance_[[*dimensions]])
            #print([pca.explained_variance_[[*dimensions]]])
            #coords = pca.explained_variance_[[*dimensions]]
            coords = pca.explained_variance_ratio_[[*dimensions]]
            #print(coords)
            ax.scatter(*coords)
            points.append(coords)
            labels.append(pca.parent_id)
        elif mode == "components":
            print(pca.components_)
            c = 0
            for comp in pca.components_:
                comp = comp[dimensions]
                print("component:", comp)
                if c in comps:
                    ax.scatter(*comp, c="C{}".format(c))
                    points.append(comp)
                    labels.append(pca.parent_id)
                c+=1
        else:
            print("mode:", mode, "not valid, available: variance, components")
            quit()

    fig.tight_layout()
    ax.set_aspect('equal')
    if cluster is None:
        if "pymol" not in sys.argv or True:

            plt.show(block=vars.block)
        return None
    else:
        print(points)
        df = pd.DataFrame(points)
        print(df)
        df["cluster"] = quick_cluster(df, n_clusters = cluster, bandwidth= bandwidth, method=method)
        df["id"] = labels
        print(df)
        plot_points(df)
        return df


def plot_points(df, title = "", labels = True):
    fig = plt.figure()
    dimensions = []
    coord_cols = []
    print(df)
    for col in df.columns:
        if type(col) is int or col[0] == "_":
            int_col = int(str(col).replace("_", ""))
            dimensions.append(int_col)
            coord_cols.append(col)

    #print(dimensions)
    #print(coord_cols)
    coord_df = df
    #print(coord_df)
    axis = [a for a in ["x", "y", "z"][:len(dimensions)]]
    coord_df.rename(columns={col:str(a) for a, col in zip(axis, coord_cols)}, inplace=True)
    print("N dimensions:", len(dimensions))
    if len(dimensions) == 3:
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.add_subplot(111)
    ax.set_title(title)
    #print(coord_df)
    '''for point, cluster in zip(df[[*dimensions]].values, df[["cluster"]].values):
        ax.scatter(*point, c="C{}".format(cluster[0]) )'''

    ax.set_xlabel(dimensions[0])
    ax.set_ylabel(dimensions[1])
    if len(dimensions) > 2:
        ax.set_zlabel(dimensions[2])
    print(coord_df)
    #print(dimensions)
    for row in coord_df.itertuples():
        #print(row)

        coords = [row.__getattribute__(a) for a in axis]

        #print(coords)
        if "cluster" in  df.columns:
            if row.cluster == -1:
                c ="black"
            else:
                c = "C{}".format(row.cluster)
            ax.scatter(*coords, c=c )
        else:
            ax.scatter(*coords)
        if labels and "id" in df.columns:
            if "id" in df.columns:
                if "cluster" in df.columns:
                    ax.text(*coords, s=row.id, c=c)
                else:
                    ax.text(*coords, s=row.id)



    fig.tight_layout()
    ax.set_aspect('equal')
    plt.savefig(os.path.join(local.temp, "temp_plot.png"))
    if not "pymol" in sys.argv or True:
        plt.show(block=vars.block)



def get_component_angles(components):
    angles = []
    for n, component in enumerate(components):
        i = n+1
        if i == len(components):
            i = 0
        a = angle_between_vectors(component, components[i])
        angles.append(a)
    return angles



def get_pca_df(in_path, subfolder, only_pcas = False, force = False, splitted=True):
    print1("Generating PCA df")
    pcas = []
    pca_df = pd.DataFrame(columns = ["id", "P0", "P1", "P2", "variance_0", "variance_1", "variance_2", "singular_0", "singular_1", "singular_2", "a01", "a12", "a20"])
    root["pcas"] = "dataframes/clustering/pcas"
    contacts_df = pd.read_csv(in_path)
    name = os.path.basename(in_path).split(".")[0]

    if splitted:
        subfolder = subfolder.format("pcas")
        os.makedirs(os.path.join(root.pcas, subfolder), exist_ok=True)
        pca_path = os.path.join(root.pcas, subfolder, name + ".csv")
        if name + ".csv" in os.listdir(os.path.join(root.pcas, subfolder)) and not force:
            print2("Skipping pca load for {}".format(name))
            return pca_path
    else:
        pca_path = os.path.join(root.pcas, name + ".csv")
        if name + ".csv" in os.listdir(os.path.join(root.pcas)) and not force:
            print2("Skipping pca load for {}".format(name))
            return pca_path



    from imports import load_single_pdb
    progress = ProgressBar(len(contacts_df.columns))
    #print(contacts_df.columns)
    for d in contacts_df.columns:
        #print1(d)
        dimers = load_single_pdb(d, pickle_folder=local.dimers, quiet=True)
        print2(dimers)
        for dimer in dimers:
            #print(dimer)
            pcas.append(dimer.pca)
            #pca_df.loc[dimer.id] = [dimer.id, *dimer.pca.components_, *dimer.pca.explained_variance_]
            angles = get_component_angles(dimer.pca.components_)
            pca_df.loc[dimer.id] = [dimer.id, *dimer.pca.components_, *dimer.pca.explained_variance_ratio_, *dimer.pca.singular_values_, *angles]
            progress.add(info=dimer.id)
    print(pca_df)

    if only_pcas:
        return pcas
    else:
        print(pca_path)
        pca_df.to_csv(pca_path)
        return pca_path




class ResPair:
    def __init__(self, res1, res2):
        self.res1 = res1
        self.res2 = res2
        self.id1 = res1.id[1]
        self.id2 = res2.id[1]
        self.id = "{}:{}".format(self.id1, self.id2)
        self.coord1 = [atom.coord for atom in self.res1.get_atoms() if atom.id == "CA"][0]
        self.coord2 = [atom.coord for atom in self.res2.get_atoms() if atom.id == "CA"][0]
        self.vector12 = vector(self.coord1, self.coord2)
        self.vector21 = vector(self.coord2, self.coord1)
        self.dist = length(self.vector12)







class ContactSurface:
    def __init__(self, structure1, structure2, outer_ids=None):
        self.structure1 = structure1
        self.structure2 = structure2
        self.outer_ids = []

        self.residues = [res.id[1] for res in structure1.get_residues()]
        if outer_ids is not None:
            self.outer_ids = outer_ids
            self.residues = [resid for resid in self.residues if resid in self.outer_ids]
        self.res_pairs, self.s_matrix, self.d_s_matrix = self.get_atom_pairs()
        #self.d_s_matrix = np.tril(self.d_t_matrix.T, -1) + self.d_t_matrix





    def get_atom_pairs(self):
        assert len(list(self.structure1.get_residues())) == len(list(self.structure2.get_residues()))
        t_matrix = []
        d_t_matrix = []
        new_pairs = {}
        res_list1 = [res for res in self.structure1.get_residues()]# if res.id[1] in self.outer_ids]
        res_list2 = [res for res in self.structure2.get_residues()]# if res.id[1] in self.outer_ids]
        assert len(res_list1) == len(res_list2)
        for n, res1 in enumerate(res_list1):
            inner1 = False
            if res1.id[1] not in self.outer_ids:
                inner1 = True
            for res2 in list(res_list2):
                inner2 = False
                if res2.id[1] not in self.outer_ids:
                    inner2 = True
                if inner1 or inner2:
                    p = None
                    dist = 666
                else:
                    p = ResPair(res1, res2)
                    dist = p.dist
                    new_pairs[p.id] = p
                if len(t_matrix) == n:
                    t_matrix.append([p])
                    d_t_matrix.append([dist])
                else:
                    t_matrix[n].append(p)
                    d_t_matrix[n].append(dist)
        print6(len(new_pairs))
        return new_pairs, np.array(t_matrix), np.array(d_t_matrix)


    @staticmethod
    def normalize_matrix(matrix, n = None):
        if n == "max":
            n = amax(matrix)
        elif n is None:
            return matrix
        max_func = lambda x: x/n
        vec_max_func = np.vectorize(max_func)
        return vec_max_func(matrix)


    @staticmethod
    def get_heat_map(matrix, title="Heat map", normalize = None, plot = True, show=False,
                     colors = None, cvals=None,
                     folder = None, percentage = False,
                     outer_ids_complete = None):

        if normalize is not None:
            if percentage:
                normalize /= 100
            matrix = ContactSurface.normalize_matrix(matrix, n=normalize)

        if outer_ids_complete is None:
            oneDmatrix1 = [mean(i) for i in matrix]
            oneDmatrix2 = [mean(i) for i in matrix.T]
            outer_ids_complete = [True] * len(oneDmatrix1)
        else:
            oneDmatrix1 = [sum(i)/len(outer_ids_complete) for i in matrix]
            oneDmatrix2 = [sum(i) / len(outer_ids_complete) for i in matrix.T]
        if plot:
            fig, axes = plt.subplots(2,2,
                                     #sharex="col",
                                     gridspec_kw={'height_ratios': [4, 2], "width_ratios": [2,4]},
                                     figsize=(12,12))

            ax = axes[0,1]
            axLeft = axes[0,0]
            axBottom = axes[1,1]
            if colors is None:
                colors = ("blue", "yellow", "red")

            if cvals is None:
                if len(colors) == 3:
                    cvals = [0, 0.5, 1]
                elif len(colors) == 2:
                    cvals = [0, 1]
            norm = plt.Normalize(min(cvals), max(cvals))
            #print(norm)
            tuples = list(zip(map(norm, cvals), colors))
            #print(tuples)
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("colormap", tuples)
            #print(cmap)
            hm = ax.imshow(matrix.T, cmap=cmap)
            max_n = len(oneDmatrix1)
            for n, (p1,p2, o) in enumerate(zip(oneDmatrix1, oneDmatrix2, outer_ids_complete)):
                if n == 0:
                    continue
                n2 = - n
                if o is None:
                    axBottom.plot((n-1,n), (oneDmatrix1[n-1], p1), c="black",  linestyle='--', linewidth=0.5)
                    axLeft.plot((oneDmatrix2[n - 1], p2), (n2 + 1, n2),  c="black", linestyle='--', linewidth=0.5)
                else:
                    axBottom.plot((n - 1, n), (oneDmatrix1[n - 1], p1), c=cmap(p1), linestyle='--', linewidth=0.5)
                    axLeft.plot((oneDmatrix2[n - 1], p2), (n2 + 1, n2), c=cmap(p2), linestyle='--', linewidth=0.5)
            #fig.tight_layout()
            fig.subplots_adjust(right=0.8)
            cbar = plt.colorbar(hm, cax=fig.add_axes([0.85, 0.15, 0.05, 0.7]))
            cbar_title = "Occurence"
            if percentage:
                cbar_title = "% "+ cbar_title
            cbar.set_label(cbar_title)

            ax.set_title(title + " N={}".format(normalize))

            if folder is None:
                root["heatmaps"] = "images/heatmaps"
                folder = root.heatmaps
            fig_path = os.path.join(folder, title + ".png")
            plt.savefig(fig_path)
            if show:
                plt.show(block=vars.block)
        return matrix, oneDmatrix1, oneDmatrix2

    @staticmethod
    def heat_map_to_pdb(matrix, structure, inplace=False, outer_list=None):

        assert len(matrix) == len(list(structure.get_residues()))

        i = 0
        if not inplace:
            structure = structure.copy()
        for res in structure.get_residues():
            ca = [atom for atom in res.get_atoms() if atom.id == "CA"][0]
            ca.bfactor = mean(matrix[i])
            i +=1
        return structure

    @staticmethod
    def display_heatmap(matrix, title, structure, n_samples=None,  show_pymol=True,obj_name =None, show_heatmap=True,
                        colors = None, cvals=None, percentage=False):
        if n_samples is not None:
            norm = n_samples
            if percentage:
                norm /= 100
        else:
            norm = None
        matrix, oneDmatrix1, oneDmatrix2 = ContactSurface.get_heat_map(matrix,
                                                         title=title,
                                                         normalize= norm,
                                                         colors=colors,
                                                         cvals=cvals,
                                                         show=show_heatmap,
                                                         percentage=percentage)

        if show_pymol:
            from pyMol import pymol_colour, pymol_temp_show
            structure = ContactSurface.heat_map_to_pdb(matrix, structure)
            name = pymol_temp_show(structure, name=obj_name)
            pymol_colour(colour="blue_yellow_red", sele=name, spectrum="b", minimum=0, maximum=amax(oneDmatrix1) / 2)


    @staticmethod
    def is_above_threshold(value, threshold, equal=True, inverse=False, as_bool=True):
        r = None
        if value == 666:
            r = False
        else:
            if equal:
                r = value >= threshold
            else:
                r = value > threshold
            if inverse:
                r= not r
        if as_bool:
            return r
        else:
            if r:
                return 1
            else:
                return 0

    def get_contact_map(self, threshold=10, as_bool=False, transposed=False):
        vec_fun = np.vectorize(self.is_above_threshold)
        contact_matrix = vec_fun(self.d_s_matrix, threshold=threshold,as_bool=as_bool, inverse=True)
        #self.get_heat_map(contact_matrix)
        if transposed:
            contact_matrix = contact_matrix.T
        return contact_matrix




def plot_cluster_heatmap(df_path):
    df = pd.read_csv(df_path)





























