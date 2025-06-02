
from scipy.linalg import expm, norm
import random as rnd
from numpy import *
import numpy as np
import math


def dihedral_angle(p0, p1, p2, p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    #p0 = p[0]
    #p1 = p[1]
    #p2 = p[2]
    #p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def dihedral_angle2(p0, p1, p2, p3):
    b0 = scale(vector(p0, p1), -1.0)
    b1 = vector(p1, p2)
    b2 = vector(p2, p3)
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))



def pnt2line(pnt, start, end):
    line_vec = vector(start, end)
    pnt_vec = vector(start, pnt)
    line_len = length(line_vec)
    line_unitvec = unit(line_vec)
    pnt_vec_scaled = scale(pnt_vec, 1.0 / line_len)
    t = dot(line_unitvec, pnt_vec_scaled)
    '''if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0'''
    nearest = scale(line_vec, t)
    dist = distance(nearest, pnt_vec)
    nearest = add(nearest, start)
    return np.array([nearest[0], nearest[1], nearest[2]])


def M(axis, theta):
    return expm(cross(eye(3), axis / norm(axis) * theta))


def rotate_around_axis(axis, theta, point):
    m = M(axis, theta)
    return dot(m, point)


def get_closest(start, end, point):
    start_end = end - start
    start_point = point - start
    return pnt2line(point, start, end)


def get_middlepoint(v1, v2):
    return (v1 + v2) / 2.0


def get_vector(start, end):
    return end - start


def test():
    start = np.array([0, 0, 0])
    end = np.array([1, 1, 1])
    point = np.array([rnd.random(), rnd.random(), rnd.random()])
    get_closest(start, end, point)


def angle_3_points(a, b, c):
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)

def angle_between_vectors(u, v):
    dot_product = sum(i * j for i, j in zip(u, v))
    norm_u = math.sqrt(sum(i ** 2 for i in u))
    norm_v = math.sqrt(sum(i ** 2 for i in v))
    cos_theta = dot_product / (norm_u * norm_v)
    angle_rad = math.acos(cos_theta)
    angle_deg = math.degrees(angle_rad)
    return angle_deg


def dot(v, w):
    x, y, z = v
    X, Y, Z = w
    return x * X + y * Y + z * Z


def length(v):
    if len(v) == 3:
        x, y, z = v
        return math.sqrt(x * x + y * y + z * z)
    elif len(v) == 2:
        x, y = v
        return math.sqrt(x * x + y * y)



def vector(b, e):
    if len(b) == 3:
        x, y, z = b
        X, Y, Z = e
        return (X - x, Y - y, Z - z)
    elif len(b) == 2:
        x, y = b
        X, Y = e
        return (X - x, Y - y)


def unit(v):
    x, y, z = v
    mag = length(v)
    return (x / mag, y / mag, z / mag)


def distance(p0, p1, **kwargs):
    return length(vector(p0, p1))

def d2(p0, p1, root = False, **kwargs):
    #print(p0, p1)
    if root:
        return sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2 + (p0[2] - p1[2])**2)
    else:
        return (p0[0] - p1[0])**2 + (p0[1] - p1[1])**2 + (p0[2] - p1[2])**2

def scale(v, sc):
    x, y, z = v
    return (x * sc, y * sc, z * sc)

def add_multiple(vectors):
    total = [0] * len(vectors[0])
    for v in vectors:
        total = add(total, v)
    return total


def add(v, w):
    #TODO: make this take any dimensions
    x, y, z = v
    X, Y, Z = w
    return (x + X, y + Y, z + Z)

def get_line(start, vector):
    end = start + vector
    return [start[0], end[0]], [start[1], end[1]], [start[2], end[2]]


def angle_dif_3d(angles1, angles2):
    if len(angles1) != len(angles2):
        print("input must be same length")
        quit()
    diffs = []
    input_len = len(angles1)
    for i in range(input_len):
        pos1 = abs(angles1[i] - angles2[i])
        if are_different_sign(angles1[i], angles2[i]):
            pos2 = 180 - abs(angles1[i]) + 180 - abs(angles2[i])
            if pos1 > pos2:
                diffs.append(pos2)
            else:
                diffs.append(pos1)
        else:
            diffs.append(pos1)
    n_diffs = []
    for n in diffs:
        # n_diffs.append(abs((n / 180) - 0.5) * 2)
        n_diffs.append(abs(n / 180))
        # n_diffs.append(n / 180)
    print(angles1, "\n", angles2)
    print(diffs)
    print(1 - (sum(n_diffs) / input_len))
    return 1 - (sum(n_diffs) / input_len)


def are_different_sign(angle1, angle2):
    if angle1 < 0 and angle2 >= 0:
        return True
    elif angle1 >= 0 and angle2 < 0:
        return True
    else:
        return False

def get_closest_point(point, points):
    shortest_dist = 999
    shortest_coords = None
    for p in points:
        #print("shorter?",point,p)
        dist = distance(p, point)
        #print(dist, shortest_dist)
        if dist < shortest_dist:
            #print("shorter:" ,dist, shortest_dist)
            shortest_dist = dist
            shortest_coords = p
    #print("shortest:", shortest_coords)
    return shortest_coords

def sub_vectors(start,end):
    new_vector = []
    for dim in range(len(start)):
        new_vector.append(start[dim] - end[dim])
    return tuple(new_vector)

def string_numbers():
    return ["1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]

def points_to_line(start,end):
    return list(zip(start,end))

def find_com(atoms):
    import types
    if isinstance(atoms, types.GeneratorType):
        atoms = list(atoms)
    x = 0
    y = 0
    z = 0
    for atom in atoms:
        if isinstance(atom, np.ndarray):
            x += atom[0]
            y += atom[1]
            z += atom[2]
        else:
            #print(com, atom.coord)
            x += atom.coord[0]
            y += atom.coord[1]
            z += atom.coord[2]
    x /= len(atoms)
    y /= len(atoms)
    z /= len(atoms)
    #print(len(atoms), end = " -> ")
    return x, y, z




def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    if any(v): #if not all zeros then
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

    else:
        return np.eye(3) #cross of all zeros only occurs on identical directions



def difference_between_boolean_pairs(A1, A2, B1, B2):
    diffX = [0,0] # Matches / Total
    diffx = [0,0]

    if A1 or B1:
        diffX[1] += 0.5
        if A1 == B1:
            diffX[0] += 0.5
    if A2 or B2:
        diffX[1] += 0.5
        if A2 == B2:
            diffX[0] += 0.5

    if A1 or B2:
        diffx[1] += 0.5
        if A1 == B2:
            diffx[0] += 0.5
    if A2 or B1:
        diffx[1] += 0.5
        if A2 == B1:
            diffx[0] += 0.5

    return diffX, diffx



