import os, sys

from utilities import *
from Globals import root, local, vars
import numpy as np
import pandas as pd


import setup

# Arcimboldo Path
#/cri4/iain/borges-arcimboldo/ARCIMBOLDO_FULL
#from cri4.iain.borges-arcimboldo.ARCIMBOLDO_FULL import SELSLIB2


def get_space_group(card):
    from spaceGroups import dictio_space_groups
    new_group = None
    raw_group = " ".join(card["group"]).strip()
    #print(raw_group, end=" -> ")
    for key, group in dictio_space_groups.items():
        if raw_group == group["symbol"]:
            new_group = group["symbol"]
            new_key = key
            break
    if new_group is None:
        print("Missing crystal group")
        quit()
    else:
        #print(new_group)
        pass
    return new_group, new_key

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

    #print1("Parsing crystalcard from", file_path)

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

def entity_to_frac(entity, params):
    for atom in entity.get_atoms():
        if atom.is_disordered():
            for d_atom in atom:
                d_atom.coord = convertFromOrthToFrac(d_atom.coord, params)
        else:
            atom.coord = convertFromOrthToFrac(atom.coord, params)
    return entity


def convertFromFracToOrth(frac_coords, parameters):
    t1, t2, t3 = frac_coords

    tz = t3 / parameters["uu"]
    ty = (t2 - tz * parameters["vv"]) / parameters["uuy"]
    tx = (t1 - ty * parameters["vvz"] - tz * parameters["uuz"]) / parameters["vvy"]

    return tx, ty, tz
def entity_to_orth(entity, params):
    for atom in entity.get_atoms():
        if atom.is_disordered():
            for d_atom in atom:
                d_atom.coord = convertFromFracToOrth(d_atom.coord, params)
        else:
            atom.coord = convertFromFracToOrth(atom.coord, params)
    return entity


def generate_displaced_copy(original, distance = 99.5, rotation = None, reverse=True):
    print3("Generating displaced copy")
    if original is None:
        return None
    displaced = original.copy()

    try:
        print(x for x in distance)
    except:
        distance = [distance] * 3
    print(distance)

    if rotation is None:
        for atom in displaced.get_atoms():
            if atom.is_disordered():
                for d_atom in atom:
                    d_atom.coord = [x+d for x, d in zip(d_atom.coord, distance)]
            else:
                atom.coord = [x+d for x, d in zip(atom.coord, distance)]
    else:
        if reverse:

            tra = [-r for r in rotation["tra"]]
            #print(rotation["rot"])
            rot = np.reshape(rotation["rot"], (1,9))
            rot = [-r for r in rot]
            rot = np.reshape(rot, (3,3))
            #print(rot)

        else:
            rot = rotation["rot"]
            tra = rotation["tra"]
        for atom in displaced.get_atoms():
            if atom.is_disordered():
                for d_atom in atom:
                    x, y, z = d_atom.coord
                    nx = (rot[0][0] * x) + (rot[0][1] * y) + (rot[0][2] * z) + tra[0]+distance[0]
                    ny = (rot[1][0] * x) + (rot[1][1] * y) + (rot[1][2] * z) + tra[1]+distance[1]
                    nz = (rot[2][0] * x) + (rot[2][1] * y) + (rot[2][2] * z) + tra[2]+distance[2]
                    d_atom.coord = [nx, ny, nz]

            else:
                x, y, z = atom.coord
                nx = (rot[0][0] * x) + (rot[0][1] * y) + (rot[0][2] * z) + tra[0]+distance[0]
                ny = (rot[1][0] * x) + (rot[1][1] * y) + (rot[1][2] * z) + tra[1]+distance[1]
                nz = (rot[2][0] * x) + (rot[2][1] * y) + (rot[2][2] * z) + tra[2]+distance[2]
                atom.coord = [nx,ny,nz]

    return displaced

def print_all_coords(entity):
    for atom in entity.get_atoms():
        print(atom.coord)

'''
def find_nearest_neighbour(original, params, key):
    print1("Finding nearest neighbour")
    print2("Space group:", key)
    from spaceGroups import dictio_space_groups
    rotation_set = dictio_space_groups[key]
    min_d2 = 0
    if len(rotation_set["symops"].items()) <= 1:
        return None
    for op_number, operation in rotation_set["symops"].items():
        print3("Symmetry operation:", op_number)
        if op_number == 1:
            continue
        displaced = generate_displaced_copy(original, distance= 99.5, rotation=operation)

        assert sum(1 for _ in original.get_atoms()) == sum(1 for _ in displaced.get_atoms())

        a = params["A"]
        b = params["B"]
        c = params["C"]
        c_a = params["c_a"]
        c_b = params["c_b"]
        c_g = params["c_g"]

        from maths import find_com
        # This shouldn't work, coms in fractional are distorted
        com = find_com(original.get_atoms())
        #print("COM:", com)
        for o_atom, d_atom in zip(original.get_atoms(), displaced.get_atoms()):
            #print(d_atom, o_atom)
            #deltaX = ((d_atom.coord[0] - o_atom.coord[0]) % 1) - 0.5
            #deltaY = ((d_atom.coord[1] - o_atom.coord[1]) % 1) - 0.5
            #deltaZ = ((d_atom.coord[2] - o_atom.coord[2]) % 1) - 0.5
            deltaX = ((d_atom.coord[0] - com[0]) % 1) - 0.5
            deltaY = ((d_atom.coord[1] - com[1]) % 1) - 0.5
            deltaZ = ((d_atom.coord[2] - com[2]) % 1) - 0.5
            #print4("Distances:",deltaX, deltaY, deltaZ)

            d2 = (a**2)*(deltaX**2) + (b**2)*(deltaY**2) + (c**2)*(deltaZ**2) +2*b*c*c_a*deltaY*deltaZ +2*a*c*c_b*deltaX*deltaZ +2*a*b*c_g*deltaX*deltaY
            #print3(d2)
            if d2 >= min_d2:
                try:
                    if d2 < o_atom.d2[0]:
                        o_atom.d2 = (d2, op_number, (deltaX,deltaY,deltaZ))
                        #print4("D2:", d2)
                    else:
                        #print4("D2:", d2, end="\r")
                        pass
                except:
                    o_atom.d2 = (d2, op_number, (deltaX,deltaY,deltaZ))
                    #print4("D2:", d2)
            else:
                #print4("D2: identity", end="\r")
                pass

    neighbour = original.copy()
    for atom in neighbour.get_atoms():
        try:
            atom.d2
        except:
            atom.d2 = (0, 1)
        atom.bfactor = atom.d2[1]
        if atom.is_disordered():
            for d_atom in atom:
                d_atom.coord = coord_add(com, atom.d2[2])

        else:
            atom.coord = coord_add(com, atom.d2[2])

    return neighbour
'''

def coord_operation_entity(entity, key, op_n, reverse = False):
    for atom in entity.get_atoms():
        if atom.is_disordered():
            for d_atom in atom:
                d_atom.coord = coord_operation(d_atom.coord, key, op_n, reverse=reverse)
        else:
            atom.coord = coord_operation(atom.coord, key, op_n, reverse=reverse)
    return entity

def coord_operation(coord, key, op_n, reverse = False):
    from spaceGroups import dictio_space_groups
    rotation = dictio_space_groups[key]
    #print(rotation)
    operation = rotation["symops"][op_n]
    rot = operation["rot"]
    tra = operation["tra"]

    x, y, z = coord
    if reverse:
        nx = (-rot[0][0] * x) + (-rot[0][1] * y) + (-rot[0][2] * z) - tra[0]
        ny = (-rot[1][0] * x) + (-rot[1][1] * y) + (-rot[1][2] * z) - tra[1]
        nz = (-rot[2][0] * x) + (-rot[2][1] * y) + (-rot[2][2] * z) - tra[2]
    else:
        nx = (rot[0][0] * x) + (rot[0][1] * y) + (rot[0][2] * z) + tra[0]
        ny = (rot[1][0] * x) + (rot[1][1] * y) + (rot[1][2] * z) + tra[1]
        nz = (rot[2][0] * x) + (rot[2][1] * y) + (rot[2][2] * z) + tra[2]

    return nx, ny, nz

def coord_add(coord, deltas, subtract = False):
    x, y, z = coord
    if subtract:
        nx = x - deltas[0]
        ny = y - deltas[1]
        nz = z - deltas[2]
    else:
        nx = x + deltas[0]
        ny = y + deltas[1]
        nz = z + deltas[2]
    return [nx, ny, nz]

def get_operation(key,op_n):
    from spaceGroups import dictio_space_groups
    return dictio_space_groups[key]["symops"][op_n]

def find_nearest_neighbour_by_chain(fractional, params, key, orth_struct):
    print1("Finding nearest neighbour")
    print2("Space group:", key)
    from spaceGroups import dictio_space_groups
    rotation_set = dictio_space_groups[key]
    min_d2 = 0
    if len(rotation_set["symops"].items()) <= 1:
        return None

    from maths import find_com
    chains = []
    for chain, orth_chain in zip(fractional.get_chains(), orth_struct.get_chains()):
        if sum([1 for _ in chain.get_residues()]) >= 100:

            orth_chain.com = find_com(orth_chain.get_atoms())
            chain.orth_com = find_com(orth_chain.get_atoms())
            chain.com = convertFromOrthToFrac(chain.orth_com, params)
            print2("COM:", [round(x, 2) for x in chain.com], [round(x, 2) for x in chain.orth_com])
            chains.append(chain)

    print2(chains)

    for op_number, operation in rotation_set["symops"].items():
        print2("Symmetry operation:", op_number)

        displaced = generate_displaced_copy(fractional, distance= 99.5, rotation=operation)

        assert sum(1 for _ in fractional.get_atoms()) == sum(1 for _ in displaced.get_atoms())

        a = params["A"]
        b = params["B"]
        c = params["C"]
        c_a = params["c_a"]
        c_b = params["c_b"]
        c_g = params["c_g"]

        for chain in chains:


            for o_atom, d_atom in zip(fractional.get_atoms(), displaced.get_atoms()):
                if op_number == 1 and o_atom.get_full_id()[2] == chain.id:
                    continue

                rel_deltaX = ((d_atom.coord[0] - o_atom.coord[0]) % 1) - 0.5
                rel_deltaY = ((d_atom.coord[1] - o_atom.coord[1]) % 1) - 0.5
                rel_deltaZ = ((d_atom.coord[2] - o_atom.coord[2]) % 1) - 0.5
                deltaX = ((d_atom.coord[0] - chain.com[0]) % 1) - 0.5
                deltaY = ((d_atom.coord[1] - chain.com[1]) % 1) - 0.5
                deltaZ = ((d_atom.coord[2] - chain.com[2]) % 1) - 0.5
                #print4("Distances:",deltaX, deltaY, deltaZ)

                d2 = (a**2)*(deltaX**2) + (b**2)*(deltaY**2) + (c**2)*(deltaZ**2) +2*b*c*c_a*deltaY*deltaZ +2*a*c*c_b*deltaX*deltaZ +2*a*b*c_g*deltaX*deltaY
                #print3(d2)
                try:
                    o_atom.d2.keys()
                except:
                    o_atom.d2 = {}
                #print(o_atom.d2.keys(), chain.id, )
                if chain.id in o_atom.d2.keys():
                    #print3(d2, o_atom.d2[chain.id][0])
                    if d2 < o_atom.d2[chain.id][0]:
                            o_atom.d2[chain.id] =  [d2, op_number, (deltaX,deltaY,deltaZ), (rel_deltaX, rel_deltaY, rel_deltaZ)]
                            #print4("D2:", d2)
                else:
                    o_atom.d2[chain.id] = [d2, op_number, (deltaX,deltaY,deltaZ), (rel_deltaX, rel_deltaY, rel_deltaZ)]
                    #print4("D2:", d2)



    print2("Generating neighbour")
    neighbour = fractional.copy()
    neighbour = entity_to_orth(neighbour, params)
    lines = []
    i = 1
    for chain in chains:
        new_model = neighbour[0].copy()
        new_model.id += i
        neighbour.add(new_model)
        i +=1
        print3(chain, new_model.id)

        for atom in new_model.get_atoms():

            ### Development
            from maths import add, sub_vectors
            # lines.append([atom.coord, sub_vectors(atom.coord,convertFromFracToOrth(atom.d2[chain.id][2], params))])
            print(atom.coord)
            #lines.append((convertFromFracToOrth(atom.coord, params),
            #              sub_vectors(convertFromFracToOrth(atom.coord, params), convertFromFracToOrth(atom.d2[chain.id][3], params))))
            lines.append((atom.coord, sub_vectors(atom.coord, atom.d2[chain.id][3])))
            # print(lines[-1])

            ###

            #atom.d2[chain.id]
            '''try:
                atom.d2[chain.id]
            except:
                #print("Atom", atom.id,"has no d2:")#, atom.__dict__)
                atom.d2[chain.id] = [0, 1, (0,0,0)]'''

            #print(atom.d2[chain.id])
            atom.bfactor = atom.d2[chain.id][1]
            orth_delta = convertFromFracToOrth(atom.d2[chain.id][2], params)
            delta = atom.d2[chain.id][2]
            rel_delta = atom.d2[chain.id][3]
            if atom.is_disordered():
                for d_atom in atom:
                    d_atom.coord = coord_add(chain.orth_com, orth_delta)
                    #d_atom.coord = coord_operation(d_atom.coord, key, atom.d2[chain.id][1])
                    #d_atom.coord = coord_add(d_atom.coord, rel_delta)
            else:
                atom.coord = coord_add(chain.orth_com, orth_delta)
                #atom.coord = coord_operation(atom.coord, key, atom.d2[chain.id][1])
                #atom.coord = coord_add(atom.coord, rel_delta)

    #neighbour = entity_to_orth(neighbour, params)
    #print_all_coords(neighbour)
    return neighbour, lines


def check_matching_coords(subset, target):
    from maths import distance
    bool_list = []
    for atom1 in subset:
        for atom2 in target:
            if atom1.id == atom2.id:
                print([round(x) for x in atom1.coord], [round(x) for x in atom2.coord], distance(atom1.coord, atom2.coord),distance(atom1.coord, atom2.coord) <= 1 )
                if distance(atom1.coord, atom2.coord) <= 1:
                    bool_list.append(True)
                else:
                    bool_list.append(False)
                break
    matching = not (False in bool_list) and len(bool_list) > 0
    return matching, bool_list


def reconstruct_relevant_neighbours(self, neighbour, params, key):
    print1("Reconstructing relevant neighbours")
    if neighbour is None:
        return None
    dimers = []
    mates = []

    original = neighbour[0]

    chains = []
    for chain in original.get_chains():
        if sum([1 for _ in chain.get_residues()]) >= 100:
            chains.append(chain)

    for chain in chains:
        model = neighbour[chains.index(chain)+1]
        print2(chain, model)
        operations = list(set([atom.d2[chain.id][1] for atom in model.get_atoms()]))
        for i, op in enumerate(operations):
            num_atoms = sum([1 for atom in model.get_atoms() if atom.d2[chain.id][1] == op])
            print(num_atoms)
            num_contacts = {}
            for c in list(set([atom.get_full_id()[2] for atom in model.get_atoms() if atom.d2[chain.id][1] == op])):
                n=0
                for atom in model.get_atoms():
                    if atom.d2[chain.id][1] != op or atom.get_full_id()[2] != c:
                        continue
                    atom.is_contact = False

                    for c_atom in chain.get_atoms():
                        from maths import distance
                        if distance(atom.coord, c_atom.coord) <= 8:
                            atom.is_contact = True
                            n += 1
                            break
                num_contacts[c] = n
            operations[i] = (op, num_atoms, num_contacts)
        print3("Operations:", operations)



        for op in operations:
            print3(op)
            for id, n_contacts in op[2].items():
                print4(id, n_contacts)
                new_chain = None
                if n_contacts < 4:
                    continue
                for c in chains:
                    if c.id == id:
                        new_chain = c.copy()
                        break
                print5(new_chain)

                if new_chain is not None:
                    model_atoms = [atom for atom in model.get_atoms() if atom.d2[chain.id][0] == op[0] and atom.is_contact]
                    '''deltas = list([atom.d2[new_chain.id][3] for atom in model_atoms])
                    [print(d) for d in deltas]'''
                    fractional = entity_to_frac(new_chain, params)
                    new_chain = generate_displaced_copy(fractional, distance = 0, rotation = get_operation(key, op[0]), reverse = False)
                    new_chain = entity_to_orth(new_chain, params)

                    print4("N atoms in model:", len(model_atoms))
                    #print([atom.is_contact for atom in model_atoms])
                    matching = check_matching_coords(model_atoms, new_chain.get_atoms())
                    print4(matching)
                    from molecules import BioObject
                    mates.append(BioObject.export(self, "mates", new_chain, "_mate_{}_{}_{}".format(chain.id, op[0], id)))


    return mates








if __name__ == "__main__":

    import setup
    from Globals import *


    if len (sys.argv) > 2:
        print(sys.argv[2:])
        if sys.argv[2] == "all":
            vars["do_only"] = []
        else:
            vars["do_only"] = sys.argv[2:]
    else:
        vars["do_only"] = ["1M2Z"]

    from imports import load_from_files
    molecules = load_from_files(local.many_pdbs, force_reload =False)
    tprint("Generating symmetries...")
    progress = ProgressBar(len(molecules))
    for molecule in molecules:
        sprint(molecule.id)
        molecule.get_all_dimers()
        molecule.pickle()
        progress.add(info=molecule.id)


    eprint("Symmetries generated")
