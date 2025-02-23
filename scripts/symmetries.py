import os, sys
from utilities import *
from Globals import root, local, vars
import numpy as np
import pandas as pd


import setup

# Arcimboldo Path
#/cri4/iain/borges-arcimboldo/ARCIMBOLDO_FULL
#from cri4.iain.borges-arcimboldo.ARCIMBOLDO_FULL import SELSLIB2



def get_cell_dim(card):
    return [card["a"], card["b"], card["c"], card["alpha"], card["beta"], card["gamma"]]

def calculate_parameters(card):
    cell_dim = get_cell_dim(card)
    parameters = {}
    parameters["A"] = A = float(cell_dim[0])
    parameters["B"] = B = float(cell_dim[1])
    parameters["C"] = C = float(cell_dim[2])
    parameters["alphaDeg"] = alphaDeg = float(cell_dim[3])
    parameters["betaDeg"] = betaDeg = float(cell_dim[4])
    parameters["gammaDeg"] = gammaDeg = float(cell_dim[5])
    parameters["alpha"] = alpha = (alphaDeg * 2 * np.pi) / 360
    parameters["beta"] = beta = (betaDeg * 2 * np.pi) / 360
    parameters["gamma"] = gamma = (gammaDeg * 2 * np.pi) / 360
    parameters["c_a"] = c_a = np.cos(alpha)
    parameters["c_b"] = c_b = np.cos(beta)
    parameters["c_g"] = c_g = np.cos(gamma)
    parameters["s_g"] = s_g = np.sin(gamma)
    parameters["q"] = q = np.sqrt(1 + 2 * c_a * c_b * c_g - c_a ** 2 - c_b ** 2 - c_g ** 2)
    parameters["uu"] = uu = s_g / (q * C)
    parameters["vv"] = vv = (c_b * c_g - c_a) / (q * B * s_g)
    parameters["uuy"] = uuy = 1 / (B * s_g)
    parameters["vvz"] = vvz = -1 * (c_g / (A * s_g))
    parameters["uuz"] = uuz = (c_a * c_g - c_b) / (q * A * s_g)
    parameters["vvy"] = vvy = 1 / A

    return parameters


def get_crystal(file_path):

    print1("Parsing crystalcard from", file_path)

    def line_contents_no_blank(line):
        l = line.split(" ")
        for c in l:
            if c == "" or c == "\n":
                l.remove(c)
        return l

    cryst1 = ""
    origx1 = ""
    origx2 = ""
    origx3 = ""
    scale1 = ""
    scale2 = ""
    scale3 = ""

    with open(file_path, "r") as file:  # Open the file in question
        for line in file:  # Check every line for the keywords
            if "CRYST1" in line: cryst1 = line
            if "ORIGX1" in line: origx1 = line_contents_no_blank(line)
            if "ORIGX2" in line: origx2 = line_contents_no_blank(line)
            if "ORIGX3" in line: origx3 = line_contents_no_blank(line)
            if "SCALE1" in line: scale1 = line_contents_no_blank(line)
            if "SCALE2" in line: scale2 = line_contents_no_blank(line)
            if "SCALE3" in line: scale3 = line_contents_no_blank(line)
    if cryst1 == "":
        return None
    crystal = {
        "a": float(clean_string(cryst1[7:16])),
        "b": float(clean_string(cryst1[16:25])),
        "c": float(clean_string(cryst1[25:34])),
        "alpha": float(clean_string(cryst1[34:41])),
        "beta": float(clean_string(cryst1[41:48])),
        "gamma": float(clean_string(cryst1[48:55])),
        "group": line_contents_no_blank(cryst1[55:67])[0:4],
        "Z": clean_string(cryst1[67:70]),
        "ori1": origx1[1:5],
        "ori2": origx2[1:5],
        "ori3": origx3[1:5],
        "scale1": scale1[1:5],
        "scale2": scale2[1:5],
        "scale3": scale3[1:5],
    }
    return crystal



def convertFromOrthToFrac(orth_coords, parameters):
    x, y, z = orth_coords

    nx = (x * parameters["vvy"]) + (y * parameters["vvz"]) + (z * parameters["uuz"])
    ny = (y * parameters["uuy"]) + (z * parameters["vv"])
    nz = z * parameters["uu"]

    return nx, ny, nz



def convertFromFracToOrth(frac_coords, parameters):
    t1, t2, t3 = frac_coords

    tz = t3 / parameters["uu"]
    ty = (t2 - tz * parameters["vv"]) / parameters["uuy"]
    tx = (t1 - ty * parameters["vvz"] - tz * parameters["uuz"]) / parameters["vvy"]

    return tx, ty, tz



def generate_displaced_copy(original, distance = 99.5):
    if original is None:
        return None
    from copy import deepcopy
    displaced = deepcopy(original)
    for atom in displaced.get_atoms():
        if atom.is_disordered():
            for d_atom in atom:
                d_atom.coord = [x+distance for x in d_atom.coord]
        else:
            atom.coord = [x+distance for x in atom.coord]
    return displaced

def find_nearest_neighbour(original,  params):
    
    displaced = generate_displaced_copy(original)

    o_atoms = original.get_atoms()
    d_atoms = displaced.get_atoms()


    assert sum(1 for _ in o_atoms) == sum(1 for _ in d_atoms)
    
    a = params["A"]
    b = params["B"]
    c = params["C"]
    c_a = params["c_a"]
    c_b = params["c_b"]
    c_g = params["c_g"]

    
    for o_atom, d_atom in zip(o_atoms, d_atoms):
        deltaX = ((d_atom[0] - d_atom[0]) % 1) - 0.5
        deltaY = ((d_atom[1] - d_atom[1]) % 1) - 0.5
        deltaZ = ((d_atom[2] - d_atom[2]) % 1) - 0.5
        d2 = a^2*deltaX^2 +b^2*deltaY^2 + c^2*deltaZ^2 +2*b*c*c_a*deltaY*deltaZ +2*a*c*c_b*deltaX*deltaZ +2*a*b*c_g*deltaX*deltaY

        if d_atom.is_disordered():
            # might not be 100% accurate but placeholder for now
            assert o_atom.is_disordered() and d_atom.is_disordered()
            for dis_d in d_atom:
                dis_d.coord[0] = atom.coord[0] + deltaX
                dis_d.coord[1] = atom.coord[1] + deltaY
                dis_d.coord[2] = atom.coord[2] + deltaZ

        else:
            d_atom.coord[0] = atom.coord[0] + deltaX
            d_atom.coord[1] = atom.coord[1] + deltaY
            d_atom.coord[2] = atom.coord[2] + deltaZ

    return displaced






if __name__ == "__main__":

    import setup
    from Globals import *


    from imports import load_from_files
    molecules = load_from_files(local.many_pdbs)
    for molecule in molecules:
        print(molecule.id)
        molecule.read_card()
        molecule.generate_fractional()
        molecule.export_fractional()
        print(molecule.fractional)

        molecule.get_neighbour()
        molecule.export_neighbour()
        print(molecule.neighbour)
        molecule.pickle()
