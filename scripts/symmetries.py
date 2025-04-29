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
    #print6("Generating displaced copy")
    if original is None:
        return None
    displaced = original.copy()

    try:
        _ = [x for x in distance]
    except:
        distance = [distance] * 3
    #print6(distance)

    if key is None and op_n in None:
        for atom in displaced.get_atoms():
            if atom.is_disordered() > 0:
                for d_atom in atom:
                    d_atom.coord = [x+d for x, d in zip(d_atom.coord, distance)]
            else:
                atom.coord = [x+d for x, d in zip(atom.coord, distance)]
    else:
        #print(key,op_n)
        coord_operation_entity(displaced, key=key, op_n=op_n, distance =distance)
        for atom in displaced.get_atoms():
            if atom.is_disordered() > 0:
                #print(atom.coord)
                #[print(d_atom.coord) for d_atom in atom]
                pass

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


class Contact:
    def __init__(self, atom, position,  target_list,
                 max_distance = None, min_contacts=0, params= None, coord = None, count_to_min = False,
                 ref_dict = None):
        #print2(atom.parent.id[1], len(target_list), type(target_list))
        if coord is None:
            self.coord = atom.coord
        else:
            self.coord = coord
        self.atom = atom
        self.max_distance = max_distance
        if self.max_distance is None:
            self.backup_distance = None
        else:
            self.backup_distance = max_distance-2
        self.min_contacts = min_contacts
        self.count_to_min = count_to_min
        self.params = params
        self.position = position["position"]
        self.target_list = target_list
        self.ref_dict = ref_dict
        self.shortest_contact = None
        self.is_contact = False
        self.is_backup = False
        self.face = None
        self.backup_face = None
        self.face_opposite = None
        self.backup_face_opposite = None
        self.num_backup_contacts = 0
        self.num_contacts = 0
        self.all_contacts = []

        self.reprocess_contacts(self.target_list, self.ref_dict)


    def __repr__(self):
        import math
        #return "Contact in c. {} r. {} ({}-{}) d: {} b: {} N:[{}/{}]".format(self.atom.get_full_id()[-2][1], self.atom.get_full_id()[-3], self.face, self. face_opposite, round(self.shortest_contact["distance"]), self.is_backup, self.num_contacts,self.num_backup_contacts)

        return "Contacts ({}) in res {} (chain {}, face: {}), shortest: {} A (res: {}, face: {})".format(self.num_contacts,
                                                                                     self.atom.get_full_id()[-2][1],
                                                                                     self.atom.get_full_id()[-3],
                                                                                     self.face,
                                                                                     round(math.sqrt(self.shortest_contact["distance"]), 2),
                                                                                     self.shortest_contact["target_atom"].parent.id[1],
                                                                                     self.shortest_contact["face"])

    def reprocess_contacts(self, target_list= None, ref_dict=None):
        if target_list is not None:
            self.get_contact_contacts(target_list)
        if ref_dict is not None:
            self.get_contact_faces(ref_dict)


    def get_contact_faces(self, ref_dict):
        face_dict = {}
        backup_face_dict = {}
        op_face_dict = {}
        op_backup_face_dict = {}
        for face, res_list in ref_dict.items():
            #print(self.atom.parent.id[1], res_list,self.atom.parent.id[1] in res_list )
            if self.atom.parent.id[1] in res_list:
                if face in face_dict.keys():
                    face_dict[face] += 1
                else:
                    face_dict[face] = 1
                    backup_face_dict[face] = 0
                if self.is_backup:
                    backup_face_dict[face] += 1

            for contact in self.all_contacts:
                if contact["target_atom"].parent.id[1] in res_list:
                    contact["face"] = face
                    if face in op_face_dict.keys():
                        op_face_dict[face] += 1
                    else:
                        op_face_dict[face] = 1
                        backup_face_dict[face] = 0
                    if contact["backup"]:
                        backup_face_dict[face] += 1

        if len(face_dict) > 0:
            self.face = sort_dict(face_dict, as_list=True)[0][0]
            if len(backup_face_dict) > 0:
                self.backup_face = sort_dict(backup_face_dict, as_list=True)[0][0]

        if len(op_face_dict) > 0:
            self.face_opposite = sort_dict(op_face_dict, as_list=True)[0][0]
            if len(op_backup_face_dict) > 0:
                self.backup_face_opposite = sort_dict(op_backup_face_dict, as_list=True)[0][0]




    def get_contact_contacts(self, target_atoms):

        if self.params is not None:
            distance_fun = get_fractional_distance
            if self.max_distance is not None:
                max_distance = self.max_distance ** 2
        else:
            from maths import d2
            distance_fun = d2
            if self.max_distance is not None:
                max_distance = self.max_distance ** 2
        if self.backup_distance is None:
            self.backup_distance = self.max_distance/2
        backup_distance = self.backup_distance ** 2
        self.is_contact = False
        self.is_backup = False
        self.num_contacts = 0
        self.all_contacts = []
        self.num_backup_contacts = 0
        for atom in target_atoms:
            #print("  ",atom.parent.id[1])
            dist = distance_fun(self.coord, atom.coord, params=self.params)
            #print(dist)
            if self.max_distance is None or dist <= max_distance:
                backup = False
                self.num_contacts += 1
                if dist <= backup_distance:
                    backup = True
                    self.is_backup = True
                    self.num_backup_contacts += 1
                line = [[c for c in self.coord], [c for c in atom.coord]]
                if self.params is not None:
                    line = [convertFromFracToOrth(coord, self.params) for coord in line]

                new_contact = {
                    "atom": self.atom,
                    "target_atom": atom,
                    "line":  line,
                    "distance": dist,
                    "backup_distance": backup_distance,
                    "maximum_distance": max_distance,
                    "face": None,
                    "backup": backup,
                }
                self.all_contacts.append(new_contact)
                if self.shortest_contact is None or dist < self.shortest_contact["distance"]:
                    self.shortest_contact = new_contact
                #print(self.num_contacts, self.min_contacts,self.num_contacts >= self.min_contacts )
                self.is_contact = True
                if self.num_contacts >= self.min_contacts:
                    #self.is_contact = True
                    if self.count_to_min:
                        #break
                        pass



def get_neigh_from_coord(coord, target_atoms, max_distance, count_to = 1, params = None): #Deprecated
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



def find_relevant_mates(self, orth_struct, params, key, minimum_chain_length = 100, contact_distance = 8, min_contacts = 0):
    print1("Finding relevant mates")
    print2("Space group:", key)

    fractional = entity_to_frac(orth_struct.copy(), params)


    from spaceGroups import dictio_space_groups
    rotation_set = dictio_space_groups[key]
    min_d2 = 0
    if len(rotation_set["symops"].items()) <= 1:
        return None
    from molecules import Mate, Monomer
    from maths import find_com
    monomers = []
    for chain, orth_chain in zip(fractional.get_chains(), orth_struct.get_chains()):
        if sum([1 for _ in chain.get_residues()]) >= minimum_chain_length:

            orth_chain.com = find_com(orth_chain.get_atoms())
            chain.orth_com = find_com(orth_chain.get_atoms())
            chain.com = convertFromOrthToFrac(chain.orth_com, params)
            #print2("COM:", [round(x, 2) for x in chain.com], [round(x, 2) for x in chain.orth_com])
            monomers.append(Monomer(self.name, chain.id, chain, self))

    print2(monomers)

    mates = []
    remaining_ids = [monomer.id for monomer in monomers]


    for fixed_monomer in monomers:
        fixed_chain = fixed_monomer.fractional_structure
        print2("Fixed monomer:", fixed_monomer.id)
        fixed_atoms = [atom for atom in fixed_chain.get_atoms()]
        com = fixed_chain.com

        for op_number, operation in rotation_set["symops"].items():
            print3("Operation:", op_number)
            displaced = generate_displaced_copy(fractional, distance= 99.5, key =key, op_n=op_number)
            #print_all_coords(displaced)
            for moving_monomer in monomers:
                if moving_monomer.id not in remaining_ids:
                    continue
                if moving_monomer.id == fixed_monomer.id and op_number == 1:
                    continue
                moving_chain = moving_monomer.fractional_structure

                #print4("Moving chain:", moving_chain)
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
                    
                    contact =  Contact(atom, position, fixed_atoms, coord = new_coord, max_distance = contact_distance, min_contacts = min_contacts, params = params)
                    if contact.is_contact:
                        if mate is None:
                            mate = Mate(op_number, operation, params, fixed_monomer, moving_monomer)
                        else:
                            if position in mate.positions.keys():
                                mate.positions[position]["contacts"].append(contact)
                                mate.positions[position]["n_contacts"] += 1
                            else:
                                mate.positions[position] = {"position": position,
                                                            "value": abs(sum(position)),
                                                            "n_contacts": 1,
                                                            "contacts": [contact]}
                if mate is None:
                    #print5("No relevant contacts")
                    continue
                print5(mate, mate.positions.keys())
                mates.append(mate)

        remaining_ids.remove(fixed_monomer.id)
    return mates







def symmetries(molecule):
    sprint(molecule.id)
    #molecule.get_all_dimers()
    molecule.pickle()



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
