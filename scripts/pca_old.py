from utilities import *
import numpy as np
import pandas as pd
from maths import *












def atoms_to_eigenvectors(atoms, parent, kwargs):
    # print("atoms_to_eigenvectors")
    atoms_df = atoms_to_df(atoms)
    #print(atoms_df)
    #print1(kwargs['name'])
    com = find_com(atoms)
    vector_df = atoms_to_vectors(atoms_df, com)
    #print(vector_df)

    pca = run_pca(vector_df, kwargs = kwargs)
    kwargs["invert"] = True
    points, components, rotated, is_inverted = process_pca(pca, vector_df, kwargs=kwargs)



    kwargs["pca"] = pca
    kwargs["original_components"] = components
    kwargs["components"] = rotated
    kwargs["points"] = points

    if kwargs["save"]:
        plot_atoms(parent, kwargs=kwargs)
    from Pairs import Axis
    original = {  # Scaled by 1/5 for representation
        "c0": Axis(com, add(com, components[0] / 5)),
        "c1": Axis(com, add(com, components[1] / 5)),
        "c2": Axis(com, add(com, components[2] / 5)),
        "centre": com
    }
    if is_inverted:
        c0o = (components[0] * -1)
        # print(parent.name, c0o, components[0])
        original["c0o"] = Axis(com, add(com, c0o / 5))
    origin = (0, 0, 0)
    rotated = {
        "Xaxis": Axis(origin, rotated[0]),
        "Yaxis": Axis(origin, rotated[1]),
        "Zaxis": Axis(origin, rotated[2]),
        "centre": origin
    }

    return pca, original, rotated


def process_pca(pca, df, kwargs):
    original_components = pca.components_
    # print("original components:", original_components)
    #print("explained variance:", pca.explained_variance_)
    scaled_components = original_components * pca.explained_variance_
    # print(df)
    is_inverted = False
    if kwargs["invert"]:
        if not is_closer_to_c(scaled_components[0], kwargs["com_c"], kwargs["com_n"])[0]:
            # print("Inverting c0", scaled_components[0])
            scaled_components[0] *= -1
            is_inverted = True
            #print("Inverted c0", scaled_components[0])

    # rotate component 0 -> Z
    rotation, points, rotated_components = rotate_all(points=df, components=scaled_components, axis=2, component=0)
    # print("rotated components:", components)
    # print(points)

    # orient component 1 -> Y
    rotation, points, rotated_components = rotate_all(points=points, components=rotated_components, axis=1, component=1,
                                                      fixed_axis=[2])
    #print("rotated components:", components)

    '''# if component 2 <0 invert all
    if False:  # components[2][0] < 0:
        rotation, points, components = rotate_all(points=points, components=components, axis=2, component=0,
                                                  inverted=True)
        rotation, points, components = rotate_all(points=points, components=components, axis=1, component=1,
                                                  fixed_axis=[2], inverted=False)
    #print("rotated components:", components)'''

    return points, scaled_components, rotated_components, is_inverted


def is_closer_to_c(component, com_c, com_n):
    dist_c = distance(component, com_c)
    dist_n = distance(component, com_n)
    if dist_c < dist_n:
        return True, dist_c, dist_n
    else:
        return False, dist_c, dist_n

def atoms_to_df(atoms):
    #print("atoms_to_df")
    df = pd.DataFrame(columns=['atom', 'X', 'Y', 'Z'])
    for atom in atoms:
        df.loc[len(df)] = [atom.id, atom.coord[0], atom.coord[1], atom.coord[2]]
    return df


def atoms_to_vectors(atoms_df, com):
    #print("atoms_to_vectors")
    vector_df = pd.DataFrame(columns=['atom', "X", "Y", "Z", "x", "y", "z"])
    for index, row in atoms_df.iterrows():
        vector_df.loc[index] = [row.atom, row.X - com[0], row.Y - com[1], row.Z - com[2],None,None,None]
    return vector_df


def run_pca(df, kwargs):
    #print("run_pca")
    pca = 1

    from sklearn.decomposition import PCA
    if "n_components" in kwargs:
        pca = PCA(n_components=kwargs['n_components'])
    else:
        pca = PCA()
    pca.fit_transform(df[["X", "Y", "Z"]])
    #print(pca)
    #print(pca.explained_variance_ratio_)
    #print(pca.singular_values_)


    return pca


def plot_atoms(parent, kwargs):
    #print("plot_atoms")
    import matplotlib.pyplot as plt
    from Globals import dirs
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #print2(kwargs['name'])
    title = "_".join([parent.name_ns, kwargs["chain"]])
    ax.set_title(title)

    ### THIS CODE RELIES ON COMPONENTS BEING IN 3 DIFFERENT AXIS! aka xyz / 90 degree rotations
    ### Now not so much but still

    points = kwargs["points"]
    original_components = kwargs["original_components"]
    components = kwargs["components"]

    # original:
    ax.scatter(points['X'], points['Y'], points['Z'], c="r")
    # rotated:
    ax.scatter(points['x'], points['y'], points['z'], c="b")

    ax.scatter(0,0,0, c="black", marker='o', s=100)
    # original:
    pca_lines = []
    pca_lines.append(points_to_line((0, 0, 0), original_components[0]))
    pca_lines.append(points_to_line((0, 0, 0), original_components[1]))
    pca_lines.append(points_to_line((0, 0, 0), original_components[2]))
    # print(pca_lines)
    for line in pca_lines:
        # print("line:", line)
        #print("line items:", line[0][0], line[0][1], line[1][0], line[1][1], line[2][0], line[2][1])
        ax.plot(line[0], line[1], line[2], color='orange')

    # rotated:
    pca_lines = []
    pca_lines.append(points_to_line((0, 0, 0), components[0]))
    pca_lines.append(points_to_line((0, 0, 0), components[1]))
    pca_lines.append(points_to_line((0, 0, 0), components[2]))
    #print(pca_lines)
    for line in pca_lines:
        # print("line:", line)
        #print("line items:", line[0][0], line[0][1], line[1][0], line[1][1], line[2][0], line[2][1])
        ax.plot(line[0], line[1], line[2], color='g')

    ax.text(x=components[0][0], y=components[0][1], z=components[0][2], s="0")
    ax.text(x=components[1][0], y=components[1][1], z=components[1][2], s="1")
    ax.text(x=components[2][0], y=components[2][1], z=components[2][2], s="2")

    #plot_ellipsoid(ax, components*pca.explained_variance_, color= "b")
    #plot_ellipsoid(ax, components*pca.explained_variance_ratio_, color="r")
    #plot_ellipsoid(ax, components*pca.singular_values_, color="g")
    #plt.show(block=True)
    fname = "{}.png".format(title)
    fdir = os.path.join(dirs.figures,"ellipsoids")
    os.makedirs(fdir, exist_ok=True)
    fig.savefig(os.path.join(fdir, fname))
    #print("Fig:", fname, "saved at", fdir)
    #save_plot(title, fig, subfolder="ellipsoids")
    #fig.savefig()
    plt.axis("equal")
    if kwargs["show"] == True:
        plt.show(block=True)
    plt.close()





def rotate_all(**kwargs):
    #print("rotate_all")
    from Maths import length
    components = kwargs['components']
    points = kwargs['points']

    which_axis = kwargs['axis']
    which_component = kwargs['component']

    component = components[which_component]
    #print("component:", component)
    if "fixed_axis" in kwargs:
        #print("fixed axis:", kwargs["fixed_axis"])
        for a in kwargs["fixed_axis"]:
            # print("a:", a)
            #print(component[a])
            component[a] = 0
    #print("component:", component)

    axis = [0, 0, 0]
    axis[which_axis] = length(component)
    if "inverted" in kwargs:
        if kwargs["inverted"]:
            axis[which_axis] = -axis[which_axis]
    # y_axis = (0, 0, length(components[0]))
    #print("axis:", axis)
    rotation = rotation_matrix_from_vectors(axis, component)  # Vector, Target (although function says opposite?
    # If wrong rotation: rotation.T could solve issues
    #print("rotation matrix:", rotation)
    points = rotate_points(points, rotation)
    vectors = rotate_vectors(components, rotation)
    return rotation, points, vectors


def rotate_points(points, rotation_matrix):
    #print("rotate_points")
    #print(points)
    for index, point in points.iterrows():
        #print(index, ":", point.X, point.Y, point.Z)
        new_point = dot(point[["X", "Y", "Z"]], rotation_matrix)
        points.loc[index, ["x", "y", "z"]] = new_point
    return points

def rotate_vectors(vectors, rotation_matrix):
    #print("rotate_vectors")
    new_vectors = []
    for vector in vectors:
        new_vectors.append(dot(vector, rotation_matrix))
    return new_vectors

def plot_ellipsoid_v2(ax, components, color = "black"):
    print("plot_ellipsoid_v2")

    pass


def plot_ellipsoid(ax, components, color = "black"):
    print("plot_ellipsoid")
    import numpy as np
    import numpy.linalg as linalg
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # your ellipsoid and center in matrix form
    A = np.array([components[0], components[1], components[2]])
    center = [0, 0, 0]

    # find the rotation matrix and radii of the axes
    U, s, rotation = linalg.svd(A)
    radii = 1.0 / np.sqrt(s)
    #radii = np.sqrt(s)
    #radii = s

    # now carry on with EOL's answer
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
    for i in range(len(x)):
        for j in range(len(x)):
            [x[i, j], y[i, j], z[i, j]] = np.dot([x[i, j], y[i, j], z[i, j]], rotation) + center
    ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color=color, alpha=0.2)