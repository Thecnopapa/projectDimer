import os

from pkg_resources import non_empty_lines

from utilities import *
from Globals import root, local, vars


import pymol


def pymol_colour_everything(start_at=0):
    # List of available colours
    colours = ['green', 'cyan', 'red', 'yellow', 'violet','blue',
               'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine',
               'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',
               'wheat', 'white', 'grey']
    ncolours = len(colours)
    num_states = pymol.cmd.count_states("(all)")
    # Loop over objects
    i = 0
    for obj in get_all_obj():
        for state in range(1,num_states):
            pymol.cmd.color(colours[i+start_at], obj+" and state "+str(state))
            i += 1
            if i >= ncolours:
                i = 0

def pymol_start(show=False, quiet=True):
    # Start PyMol
    sprint("Starting PyMol...")
    if show:  # Debug stuff / -c: no interface, -q no startup message, -Q completely quiet
        pymol.finish_launching(["pymol"])  #
        if quiet:
            pymol.finish_launching(["pymol", "-q"])
    else:
        pymol.finish_launching(["pymol", "-cqQ"])

def pymol_close():
    sprint("Closing PyMol...")
    pymol.cmd.quit()

def pymol_reset():
    pymol.cmd.delete("All")

def pymol_load_name(file_name, folder):
    pymol.cmd.load(os.path.join(folder,file_name), file_name)
    return file_name

def pymol_load_path(path,  name=None, state = -1,):
    if name is None:
        name = os.path.basename(path)
    pymol.cmd.load(path,name)
    pymol.cmd.set("state", state)
    return name

def pymol_hide(sele, hide = None):
    all_reprs = ["lines","spheres","mesh","ribbon","cartoon","sticks","dots","surface","labels","nonbonded","nb_spheres"]
    if hide is not None:
        if hide == "all":
            hide = " ".join(all_reprs)
        pymol.cmd.hide(representation=hide, selection=sele)

def pymol_format(representation,identifier="", hide="all", colour = None, spectrum=None):
    for obj in get_all_obj():
        if identifier in obj:
            pymol_hide(obj, hide=hide)
            pymol.cmd.show(representation=representation, selection=obj)
            if colour is not None:
                pymol_colour(colour, obj=obj, spectrum=spectrum)


def pymol_colour(colour, obj = "(all)", sele = None, spectrum=None):
    if sele is not None:
        sele_str = "({} and {})".format(obj,sele)
    else:
        sele_str = "({})".format(obj)
    if spectrum is not None:
        pymol.cmd.spectrum(spectrum, colour, sele_str)
    else:
        print("colouring", sele_str, colour)
        pymol.cmd.colour(colour, sele_str)


def pymol_save_small(file_name, folder, dpi=300, height=100, width=100, save_session=None):
    image_path = os.path.join(folder, file_name+".png")
    pymol.cmd.png(image_path, width=width, height=height, dpi=dpi)
    return image_path

def pymol_save_session(file_name, folder):
    path = os.path.join(folder, file_name + ".pse")
    pymol.cmd.save(path)
    return path

def get_all_obj():
    return pymol.cmd.get_names(type='objects')

def generate_preview(path, folder="", state = 0,save_session=True):
    if save_session:
        local[folder+"_sessions"] = "sessions/{}_sessions".format(folder)
    local[folder] = "previews/{}".format(folder)
    name = os.path.basename(path).split(".")[0]
    pymol_reset()
    pymol_load_path(path, state=state)
    pymol.cmd.orient("(all)")
    if state == 0:
        pymol_colour_everything()
    else:
        pymol_colour_everything(start_at=state-1)
    if save_session:
        pymol_save_session(name, local[folder+"_sessions"])
    preview_path = pymol_save_small(name.upper(), local[folder], dpi=50, height=150, width=150)
    return preview_path

def pymol_align_all():
    all_obj = pymol.cmd.get_names(type='objects')
    print(all_obj)
    obj1 = all_obj[0]
    for obj2 in all_obj:
        if obj2 != obj1:
            pymol_align__obj(obj1, obj2)
    pymol.cmd.orient("(all)")


def pymol_align__obj(obj1, obj2):
    pymol.cmd.align(obj2, obj1)
    print3("aligned:", obj1, obj2)
    pymol.cmd.orient("(all)")


def pymol_align_chains(chains_to_align):
    all_obj = pymol.cmd.get_names(type='objects')
    print1(all_obj)
    obj1, chain1 = all_obj[0], chains_to_align[0]
    sele1 = "{} and c. {}".format(obj1, chain1)
    for obj2, chain2 in zip(all_obj, chains_to_align):
        if obj2 != obj1:
            sele2 = "{} and c. {}".format(obj2, chain2)
            pymol_align__obj(sele1, sele2)


def pymol_symmetries(obj = None):
    if obj is None:
        obj = pymol.cmd.get_names(type='objects')[0]
    pymol.cmd.symexp("sym", obj, obj, 6)

def pymol_set_state(state):
    pymol.cmd.set("state", state)

def pymol_colour_all(colour):
    pymol.cmd.color(colour, "(all)")

def pymol_orient():
    pymol.cmd.orient("(all)")

def pymol_group(identifier = "sym", name = "symmetries"):
    group = []
    for obj in pymol.cmd.get_names(type='objects'):
        if identifier in obj:
            group.append(obj)

    pymol.cmd.group(name, " ".join(group))

def pymol_draw_line(coord1, coord2, name = "d", state = -1):
    #print("Distance between:", coord1, "and", coord2)
    pymol.cmd.pseudoatom("tmp1", pos=coord1, state=state)
    pymol.cmd.pseudoatom("tmp2", pos=coord2, state=state)
    pymol.cmd.distance(name, "tmp1","tmp2", state=state)
    pymol.cmd.delete("tmp1")
    pymol.cmd.delete("tmp2")




















