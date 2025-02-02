import os
import globals
from globals import *
from utilities import *

if os.name == 'nt':
    globals.vars["gesamt"] = "C:\Program Files\CCP4-9\CCP4/bin\gesamt.exe"

def superpose_single(id, fixed, moving):
    import subprocess
    local["superposed_monomers"] = "superposed/superposed_monomers"
    os.makedirs(local.superposed_monomers, exist_ok=True)
    out_path =  os.path.join(local["superposed_monomers"], id + ".pdb")
    super_line = [globals.vars.gesamt, fixed, moving, "-o", out_path ]
    #print(super_line)
    gesamt_out = subprocess.run(super_line, capture_output=True, text=True)
    #print(gesamt_out)
    data = {}
    t_matrix_lines = 0
    for line in gesamt_out.stdout.splitlines():
        if "Q-score" in line:
            data["q_score"] = float(line.split(":")[1])
        if "RMSD" in line:
            data["rmsd"] = float(line.split(":")[1])
        if "Aligned residues" in line:
            data["aligned_residues"] = int(line.split(":")[1])
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
    return data












