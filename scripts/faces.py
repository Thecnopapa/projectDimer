import os
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
    matplotlib.use('TkAgg')
except:
    matplotlib.use('QtAgg')


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

GR_colours = {"front": "purple",
              "base": "salmon",
              "back": "cyan",
              "top": "pink",}


def define_faces_from_list(self, list):
    for contact in self.contacts:
        print(contact)
        pass


def get_pca(structure, n_components = 3, com = None, closer_to = None):
    print3("Getting PCA")
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_components, random_state=6)
    coords = [atom.coord for atom in structure.get_atoms()]
    if com is None:
        com = find_com(coords)
    if closer_to is None:
        closer_to = get_terminals(structure)["N"]
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
                print4("Reverse component", n, pca.components_[n], "-->", end=" ")
                pca.components_[n] = component * -1
                print(pca.components_[n])



    return pca

def pca_to_lines(pca, com, just_points = False):
    components = pca.components_ * pca.explained_variance_ratio_ * pca.singular_values_
    points = []
    lines = []
    print("COmponents:", components)
    for component in components:
        c = [component[i] + com[i] for i in range(len(component))]
        print( "C:", c)
        points.append((com, c))
        lines.append(points_to_line(com, c))
        print(points)
    if just_points:
        return points
    else:
        return lines

def get_terminals(structure):
    atom_list = list(structure.get_atoms())
    terminals = dict(N = atom_list[0].coord,
                     C = atom_list[-1].coord)
    return terminals



def plot_atoms(structure, pca = None):
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
    print("COM:", com)
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
        original_components = pca.components_ * pca.explained_variance_ratio_*pca.singular_values_
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
    plt.show(block=True)


























