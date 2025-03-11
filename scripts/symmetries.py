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
        if atom.is_disordered() > 0:
            for d_atom in atom:
                d_atom.coord = convertFromOrthToFrac(d_atom.coord, params)
        #else:
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
        if atom.is_disordered() > 0:
            for d_atom in atom:
                d_atom.coord = convertFromFracToOrth(d_atom.coord, params)
        #else:
        atom.coord = convertFromFracToOrth(atom.coord, params)
    return entity


def generate_displaced_copy(original, distance = 99.5, key = None, op_n = None):
    print6("Generating displaced copy")
    if original is None:
        return None
    displaced = original.copy()

    try:
        _ = [x for x in distance]
    except:
        distance = [distance] * 3
    print6(distance)

    if key is None and op_n in None:
        for atom in displaced.get_atoms():
            if atom.is_disordered() > 0:
                for d_atom in atom:
                    d_atom.coord = [x+d for x, d in zip(d_atom.coord, distance)]
            else:
                atom.coord = [x+d for x, d in zip(atom.coord, distance)]
    else:
        coord_operation_entity(displaced, key=key, op_n=op_n, distance =distance)
        for atom in displaced.get_atoms():
            if atom.is_disordered() > 0:
                print(atom.coord)
                [print(d_atom.coord) for d_atom in atom]

    return displaced

def print_all_coords(entity, head = 5):
    n = 0
    for atom in entity.get_atoms():
        print(atom.coord)
        n += 1
        if n >= head:
            break



def coord_operation_entity(entity, key, op_n, distance = (0,0,0)):
    for atom in entity.get_atoms():
        if atom.is_disordered() > 0:
            #print(atom.is_disordered())
            for d_atom in atom:
                #print(d_atom.get_full_id())
                d_atom.coord = coord_operation(d_atom.coord, key, op_n, distance =distance)
        #else:
        atom.coord = coord_operation(atom.coord, key, op_n, distance=distance)
    return entity

def coord_operation(coord, key, op_n, distance = (0,0,0)):
    from spaceGroups import dictio_space_groups
    rotation = dictio_space_groups[key]
    #print(rotation)
    operation = rotation["symops"][op_n]
    rot = operation["rot"]
    tra = operation["tra"]

    x, y, z = coord

    nx = (rot[0][0] * x) + (rot[0][1] * y) + (rot[0][2] * z) + tra[0]+distance[0]
    ny = (rot[1][0] * x) + (rot[1][1] * y) + (rot[1][2] * z) + tra[1]+distance[1]
    nz = (rot[2][0] * x) + (rot[2][1] * y) + (rot[2][2] * z) + tra[2]+distance[2]

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

def get_neigh_from_coord(coord, target_atoms, max_distance, count_to = 1, params = None):
    if params is not None:
        distance_fun = get_fractional_distance
        max_distance = max_distance**2
    else:
        from maths import distance as distance_fun
    num_contacts = 0
    contacts = []
    is_contact = False
    for atom in target_atoms:
        if distance_fun(coord, atom.coord, params = params) <= max_distance:
            num_contacts += 1
            contacts.append([[c for c in coord], [c for c in atom.coord]])
            if num_contacts >= count_to:
                is_contact = True
                break
    return is_contact, num_contacts, contacts
                    
    
def get_fractional_distance(coord1, coord2, params):
    from maths import vector
    deltaX, deltaY, deltaZ = vector(coord1, coord2)

    a = params["A"]
    b = params["B"]
    c = params["C"]
    c_a = params["c_a"]
    c_b = params["c_b"]
    c_g = params["c_g"]

    d2 = (a**2)*(deltaX**2) + (b**2)*(deltaY**2) + (c**2)*(deltaZ**2) +2*b*c*c_a*deltaY*deltaZ +2*a*c*c_b*deltaX*deltaZ +2*a*b*c_g*deltaX*deltaY

    return d2



def find_relevant_mates(orth_struct, params, key):
    print1("Finding relevant mates")
    print2("Space group:", key)

    fractional = entity_to_frac(orth_struct.copy(), params)


    from spaceGroups import dictio_space_groups
    rotation_set = dictio_space_groups[key]
    min_d2 = 0
    if len(rotation_set["symops"].items()) <= 1:
        return None, None

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

    mates = []
    remaining_ids = [chain.id for chain in chains]
    from molecules import Mate

    for fixed_chain in chains:
        print2("Fixed chain:", fixed_chain)
        fixed_atoms = [atom for atom in fixed_chain.get_atoms()]
        com = fixed_chain.com
          
        for op_number, operation in rotation_set["symops"].items():
            print3("Operation:", op_number)
            displaced = generate_displaced_copy(fractional, distance= 99.5, key =key, op_n=op_number)
            #print_all_coords(displaced)

            for moving_chain in chains:
                if moving_chain.id not in remaining_ids:
                    continue
                if moving_chain.id == fixed_chain.id and op_number == 1:
                    continue

                print4("Moving chain:", moving_chain)
                #print_all_coords(moving_chain)
                
                moving_atoms = [atom for atom in displaced.get_atoms() if atom.get_full_id()[-3] == moving_chain.id]
                #print(moving_atoms)
                
                mate = None

                for atom in moving_atoms:
                    #print(atom, atom.get_full_id(), atom.coord, atom.is_disordered()>0)
                    deltaX = ((atom.coord[0] - com[0]) % 1) - 0.5
                    deltaY = ((atom.coord[1] - com[1]) % 1) - 0.5
                    deltaZ = ((atom.coord[2] - com[2]) % 1) - 0.5

                    new_coordX = com[0] + deltaX
                    new_coordY = com[1] + deltaY
                    new_coordZ = com[2] + deltaZ
                    new_coord = [new_coordX, new_coordY, new_coordZ]

                    position= [(n_coord - d_coord + 99.5) for n_coord, d_coord in zip(new_coord, atom.coord)]
                    for p in position:
                        assert p % 1 == 0
                    #print(position)

                    position = tuple([int(p) for p in position])
                    if any([p >= 10 for p in position]):
                        print(atom.get_full_id())
                        print("Position:", position)
                        print("Original:", atom.coord)
                        #print("Deltas:", deltaX, deltaY, deltaZ)
                        print("New coord:", new_coord)
                        quit()
                    
                    is_contact, _, contacts =  get_neigh_from_coord(new_coord, fixed_atoms, max_distance = 8, params = params)
                    if is_contact:
                        if mate is None:
                            mate = Mate(op_number, operation, params, fixed_chain, moving_chain)
                        else:
                            if position in mate.positions.keys():
                                mate.positions[position]["n_contacts"] += 1
                                mate.positions[position]["contacts"].extend(contacts)
                            else:
                                mate.positions[position] = {"position": position,
                                                            "n_contacts": 1,
                                                            "contacts": contacts}
                if mate is None:
                    print5("No relevant contacts")
                    continue
                print5(mate, mate.positions.keys())
                mates.append(mate)

        remaining_ids.remove(fixed_chain.id)
    return mates







    
def deprecated():  ##########
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

            

            mate = Mate(op_number, operation, params, chain) 
            for o_atom, d_atom in zip(fractional.get_atoms(), displaced.get_atoms()):
                if op_number == 1 and o_atom.get_full_id()[2] == chain.id:
                    continue

                rel_deltaX = ((d_atom.coord[0] - o_atom.coord[0]) % 1) - 0.5
                rel_deltaY = ((d_atom.coord[1] - o_atom.coord[1]) % 1) - 0.5
                rel_deltaZ = ((d_atom.coord[2] - o_atom.coord[2]) % 1) - 0.5
                deltaX = ((d_atom.coord[0] - chain.com[0]) % 1) - 0.5
                deltaY = ((d_atom.coord[1] - chain.com[1]) % 1) - 0.5
                deltaZ = ((d_atom.coord[2] - chain.com[2]) % 1) - 0.5
                new_coordX = o_atom.coord[0] + rel_deltaX
                new_coordY = o_atom.coord[1] + rel_deltaY
                new_coordZ = o_atom.coord[2] + rel_deltaZ
                new_coord = [new_coordX, new_coordY, new_coordZ]
                if op_number == 99:
                    print4("Distances (op: {}):". format(op_number),new_coordX, new_coordY, new_coordZ, end = " -> ")

                position = [int(n_coord - d_coord + 99.5) for n_coord, d_coord in zip(new_coord, d_atom.coord)]
                position = tuple(position)

                if op_number == 99:
                    print(position)

                d2 = (a**2)*(deltaX**2) + (b**2)*(deltaY**2) + (c**2)*(deltaZ**2) +2*b*c*c_a*deltaY*deltaZ +2*a*c*c_b*deltaX*deltaZ +2*a*b*c_g*deltaX*deltaY
                #print3(d2)
                try:
                    o_atom.d2.keys()
                except:
                    o_atom.d2 = {}
                #print(o_atom.d2.keys(), chain.id, )
                if chain.id in o_atom.d2.keys():
                    #print3(d2, o_atom.d2[chain.id][0])
                    if d2 < o_atom.d2[chain.id]["distance"]:
                            o_atom.d2[chain.id] =  {"distance": d2,
                                                    "operation":op_number,
                                                    "deltas": (deltaX,deltaY,deltaZ),
                                                    "rel_deltas": (rel_deltaX, rel_deltaY, rel_deltaZ),
                                                    "new_coords": (new_coordX, new_coordY, new_coordZ),
                                                    "position": position}
                            #print4("D2:", d2)
                else:
                    o_atom.d2[chain.id] = {"distance": d2,
                                           "operation":op_number,
                                           "deltas": (deltaX,deltaY,deltaZ),
                                           "rel_deltas": (rel_deltaX, rel_deltaY, rel_deltaZ),
                                           "new_coords": (new_coordX, new_coordY, new_coordZ),
                                           "position": position}
                    #print4("D2:", d2)


    
    print2("Generating neighbour")
    neighbour = fractional.copy()
    #neighbour = entity_to_orth(neighbour, params)
    all_lines = []
    i = 1
    for chain in chains:
        lines = []
        new_model = neighbour[0].copy()
        new_model.id += i
        neighbour.add(new_model)
        i +=1
        print3(chain, new_model.id)

        for atom in new_model.get_atoms():

            #print(atom.d2[chain.id])
            atom.bfactor = atom.d2[chain.id]["operation"]
            orth_delta = convertFromFracToOrth(atom.d2[chain.id]["deltas"], params)
            delta = atom.d2[chain.id]["deltas"]
            rel_delta = atom.d2[chain.id]["rel_deltas"]
            #print(atom.d2[chain.id][4])
            lines.append([convertFromFracToOrth(atom.coord, params)])
            if atom.is_disordered() > 0:
                for d_atom in atom:
                    #d_atom.coord = coord_add(chain.orth_com, orth_delta)
                    #d_atom.coord = coord_operation(d_atom.coord, key, atom.d2[chain.id]["operation"])
                    #d_atom.coord = coord_add(d_atom.coord, atom.d2[chain.id]["position"])
                    d_atom.coord = coord_add(d_atom.coord, rel_delta)
                    #d_atom.coord = convertFromFracToOrth(atom.d2[chain.id][4], params)
            else:
                #atom.coord = coord_add(chain.orth_com, orth_delta)
                #atom.coord = coord_operation(atom.coord, key, atom.d2[chain.id]["operation"])
                #atom.coord = coord_add(atom.coord, atom.d2[chain.id]["position"])
                atom.coord = coord_add(atom.coord, rel_delta)
                #atom.coord = convertFromFracToOrth(atom.d2[chain.id][4], params)
            lines[-1].append(convertFromFracToOrth(atom.coord, params))
            #print(atom.coord)
        all_lines.append(lines)
    neighbour = entity_to_orth(neighbour, params)
    #print_all_coords(neighbour)
    return neighbour, all_lines


def check_matching_coords(subset, target):
    from maths import distance
    bool_list = []
    for atom1 in subset:
        for atom2 in target:
            if atom1.get_full_id()[-2] == atom2.get_full_id()[-2] and atom1.get_full_id()[-3] == atom2.get_full_id()[-3]:
                
                if distance(atom1.coord, atom2.coord) <= 1:
                    bool_list.append(True)
                else:
                    bool_list.append(False)
                    print(atom1.get_full_id(), atom2.get_full_id())
                    print([round(x) for x in atom1.coord], [round(x) for x in atom2.coord], distance(atom1.coord, atom2.coord),distance(atom1.coord, atom2.coord) <= 1 )
                break
    matching = not (False in bool_list) and len(bool_list) > 0
    return matching, bool_list


def reconstruct_relevant_neighbours(neighbour, params, key):
    print1("Reconstructing relevant neighbours")
    from maths import distance
    from molecules import BioObject, Monomer, Dimer

    if neighbour is None:
        return None
    dimers = []
    monomers = []
    all_mates = []
    all_contacts = []

    original = neighbour[0]

    chains = []
    for chain in original.get_chains():
        if sum([1 for _ in chain.get_residues()]) >= 100:
            chains.append(chain)
    remaining_ids = [chain.id for chain in chains]
    for chain in chains:
        m_n = chains.index(chain)+1
        model = neighbour[m_n]
        print2("Chain:", chain, model)
        operations = list(set([atom.d2[chain.id]["operation"] for atom in model.get_atoms()]))
        mates = []
        contacts = []
        
        for i, op in enumerate(operations):
            print3( "Operation:", i, op)
            op_atoms = [atom for atom in model.get_atoms() if atom.d2[chain.id]["operation"] == op]
            num_atoms = len(op_atoms)
            positions = set([atom.d2[chain.id]["position"] for atom in op_atoms])

            for position in positions:
                op_atoms_by_pos = [atom for atom in op_atoms if atom.d2[chain.id]["position"] == position]
                print4("Position: {} (N:{})".format(position, len(op_atoms_by_pos)))

                new_chain_ids_by_pos = set([atom.get_full_id()[2] for atom in op_atoms_by_pos if atom.get_full_id()[2] in [chain.id for chain in chains]])
                #print5("New chain ids", new_chain_ids_by_pos) 
                for new_chain_id in new_chain_ids_by_pos:
                    if new_chain_id not in remaining_ids:
                        continue
                    print5("New chain id:", new_chain_id)
                    op_atoms_by_pos_by_chain = [atom for atom in op_atoms_by_pos if atom.get_full_id()[2] == new_chain_id]
                    

                    num_contacts = 0
                
                    for atom in op_atoms_by_pos_by_chain:
                        atom.is_contact = False
                        for o_atom in chain.get_atoms():
                            if distance(atom.coord, o_atom.coord) <= 8:
                                atom.is_contact = True
                                num_contacts += 1
                                contacts.append([[np.float64(c) for c in o_atom.coord] , [np.float64(c) for c in atom.coord]])
                                
                                break

                    
                    
                    if num_contacts < 5:
                        print6("Not enough contacts ({})".format(num_contacts))
                        continue
                    else:
                        print6("N: {}, Contacts: {}".format(len(op_atoms_by_pos_by_chain), num_contacts))

                    mate_chain = None
                    
                    for c in chains:
                        if c.id == new_chain_id:
                            mate_chain = c.copy()

                    if mate_chain is None:
                        print6("Error finding reference chain")
                        continue
                    mate_chain = entity_to_frac(mate_chain, params)
                    mate_chain = generate_displaced_copy(mate_chain, distance = position, rotation = get_operation(key, op), reverse = False)
                    mate_chain = entity_to_orth(mate_chain, params)
                    matching = check_matching_coords(op_atoms_by_pos_by_chain, mate_chain.get_atoms())
                    if not matching[0]:
                        print6("Mate does not match closest atoms")
                        print(matching)
                    mates.append(BioObject.export(self, "mates", mate_chain, "_mate_{}_{}_{}_{}".format(chain.id, i, new_chain_id, position)))
                    print6("Mate generated succesfully:", mate_chain.get_full_id()) 
                    dimer = chain.copy()
    


        remaining_ids.remove(chain.id)

        all_mates.append(mates)




        all_contacts.append(contacts)
    #self.contacts = all_contacts
    return {"mates": all_mates,
            "contacts": all_contacts,
            "monomers": monomers,
            "dimers": dimers}



def symmetries(molecule, progress):
    sprint(molecule.id)
    molecule.get_all_dimers()
    molecule.pickle()
    progress.add(info=molecule.id)


if __name__ == "__main__":

    import setup
    from Globals import *


    if len (sys.argv) > 2:
        print(sys.argv)
        if "all" in sys.argv:
            vars["do_only"] = []
        else:
            vars["do_only"] = sys.argv[2:]
    else:
        vars["do_only"] = ["1M2Z"]

    from imports import load_from_files
    
    molecules = load_from_files(local.many_pdbs, force_reload = False)
    tprint("Generating symmetries...")
    progress = ProgressBar(len(molecules))

    MAX_WORKERS = 0
    if MAX_WORKERS <= 1:
        for molecule in molecules:
            symmetries(molecule, progress)
    else:

        import concurrent.futures
        pool = concurrent.futures.ThreadPoolExecutor(max_workers=32)
        for molecule in molecules:
            pool.submit(symmetries, molecule, progress)
        pool.shutdown(wait = True)

    eprint("Symmetries generated")
