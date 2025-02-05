import os

import scipy.stats

from utilities import *
from globals import root, local, vars


import pymol


def pymol_start(show=False, quiet=True):
    # Start PyMol
    sprint("Starting PyMol...")
    if show:  # Debug stuff / -c: no interface, -q no startup message, -Q completely quiet
        pymol.finish_launching(["pymol"])  #
        if quiet:
            pymol.finish_launching(["pymol", "-q"])
    else:
        pymol.finish_launching(["pymol", "-cqQ"])

def pymol_reset():
    pymol.cmd.delete("All")

def pymol_load_name(file_name, folder):
    pymol.cmd.load(os.path.join(local[folder], file_name),file_name)

def pymol_load_path(path):
    print(path, os.path.basename(path),"\n")
    pymol.cmd.load(path,os.path.basename(path), state=0)


def pymol_save_small(file_name, folder, dpi=300, height=100, width=100):
    image_path = os.path.join(folder, file_name+".png")
    pymol.cmd.png(image_path, width=width, height=height, dpi=dpi)
    return image_path

def generate_preview(path, folder=""):
    try:
        root[folder] = "previews/{}".format(folder)
        name = os.path.basename(path).split(".")[0]
        pymol_reset()
        pymol_load_path(path)
        preview_path = pymol_save_small(name, root[folder], dpi=50, height=150, width=150)
        return preview_path
    except:
        print("Failed to generate preview for {} ({})".format(name,path))
        return None