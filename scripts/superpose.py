import os
import Globals
from Globals import *
from utilities import *
import platform
if os.name == 'nt':
    Globals.vars["gesamt"] = "C:/Program Files/CCP4-9/CCP4/bin/gesamt.exe"
else:
    Globals.vars["gesamt"] = "gesamt"

def superpose_single(id, fixed, moving):
    import subprocess
    local["super_raw"] = "superposed/super_raw"
    os.makedirs(local.super_raw, exist_ok=True)
    #print(id)
    out_path =  os.path.join(local["super_raw"], id + ".pdb")
    super_line = [Globals.vars.gesamt, fixed, moving, "-o", out_path ]
    #print(super_line)
    gesamt_out = subprocess.run(super_line, capture_output=True, text=True)
    if "show_gesamt" in vars:
        if not(len(vars.do_only) == 0 or vars.do_only is None):
            print(gesamt_out.stdout)
            print(gesamt_out.stderr)

    data = {"out_path": out_path}
    t_matrix_lines = 0
    centroid_lines = 0
    data["map"] =[]
    for line in gesamt_out.stdout.splitlines():
        #print(centroid_lines, line)
        if "Q-score" in line:
            data["q_score"] = float(line.split(":")[1])
        if "RMSD" in line:
            data["rmsd"] = float(line.split(":")[1])
        if "Aligned residues" in line:
            data["aligned_residues"] = int(line.split(":")[1])
        if "Sequence Id:" in line:
            data["identity"] = float(line.split(":")[2])
        if t_matrix_lines > 0:
            l = line.split()
            data["t_matrix"]["Rx"].append(float(l[0]))
            data["t_matrix"]["Ry"].append(float(l[1]))
            data["t_matrix"]["Rz"].append(float(l[2]))
            data["t_matrix"]["T"].append(float(l[3]))
            t_matrix_lines -= 1
        if "Rx" in line:
            data["t_matrix"] = {"Rx": [],"Ry": [],"Rz": [], "T": []}
            t_matrix_lines = 3
        if "FIXED |" in line:
            data["nres"] = int(line.split("|")[1])
        if "MOVING |" in line:
            data["ref_nres"] = int(line.split("|")[1])

        if centroid_lines > 0 :
            if "FIXED" in line:
                l = line.split()
                data["centroids"]["self"] = [float(i) for i in l[1:4]]
            elif "MOVING" in line:
                l = line.split()
                data["centroids"]["reference"] = [float(i) for i in l[1:4]]
            elif "Distance" in line:
                data["centroids"]["distance"] = float(line.split(":")[1])

            elif "cosines" in line:
                l = line.split(":")[1]
                l = l.split()
                data["centroids"]["cosines"] = [float(i) for i in l]
            elif "Angle between" in line:
                data["centroids"]["angle"] = float(line.split(":")[1])
        if "CENTROIDS" in line:
            centroid_lines = 8
            data["centroids"] = {}

        if "Polar angles" in line:
            l = line.split(":")
            data["ccp4_angles"] = {"polar": [float(i) for i in l[1].split()]}
        if "Euler angles" in line:
            l = line.split(":")
            data["ccp4_angles"]["euler"] = [float(i) for i in l[1].split()]
        if "Orthogonal translation" in line:
            l = line.split(":")
            data["ccp4_angles"]["translation"] = [float(i) for i in l[1].split()]

        if line.endswith("|") and line.startswith("|"):
            l = line.split("|")
            if len(l) != 5:
                continue
            res1 = l[1].split(" ")[-2]
            dist = l[2][4:8]
            res2 = l[3].split(" ")[-2]
            if res1 == "":
                res1 = None
            else:
                res1 = get_digits(res1, integer= True)
            if res2 == "":
                res2 = None
            else:
                res2 = get_digits(res2, integer=True)
            if dist == "    ":
                dist = None
            #print(res1, dist, res2)
            data["map"].append({"res1": res1, "res2": res2, "distance": dist})

        centroid_lines -=1

    return data

def create_maps(fixed, moving, raw_map):
    fixed_to_moving = {}
    moving_to_fixed = {}
    fixed_residues = {res.id[1]:res for res in fixed.get_residues()}
    moving_residues = {res.id[1]:res for res in moving.get_residues()}
    for row in raw_map:
        if row["res1"] is not None:
            if row["res2"] is not None:
                fixed_to_moving[row["res1"]] = moving_residues[row["res2"]]
            else:
                fixed_to_moving[row["res1"]] = None

        if row["res2"] is not None:
            if row["res1"] is not None:
                moving_to_fixed[row["res2"]] = fixed_residues[row["res1"]]
            else:
                moving_to_fixed[row["res1"]] = None



    return fixed_to_moving, moving_to_fixed














