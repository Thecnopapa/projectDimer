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
    if vars.pymol_started:
        if not quiet:
            print("(PyMol) PyMol already started")
        return
    sprint("(PyMol) Starting PyMol...")
    if show:  # Debug stuff / -c: no interface, -q no startup message, -Q completely quiet
        pymol.finish_launching(["pymol"])  #
        if quiet:
            pymol.finish_launching(["pymol", "-q"])
    else:
        pymol.finish_launching(["pymol", "-cqQ"])
    vars["pymol_started"] = True

def pymol_close():
    sprint("(PyMol) Closing PyMol...")
    pymol.cmd.quit()

def pymol_reset():
    pymol.cmd.delete("All")

def pymol_load_name(file_name, folder):
    pymol.cmd.load(os.path.join(folder,file_name), file_name)
    return file_name

def pymol_load_path(path,  name=None, state = -1,):
    if name is None:
        name = os.path.basename(path)
    print("(PyMol) Loading:",name)
    pymol.cmd.load(path,name)
    #pymol.cmd.set("state", state)
    return name

def pymol_disable(sele, silent=True):
    if not silent:
        print("(PyMol) Disabling:", sele)
    pymol.cmd.disable(sele)

def pymol_hide(sele, hide = None, silent = True):
    all_reprs = ["lines","spheres","mesh","ribbon","cartoon","sticks","dots","surface","labels","nonbonded","nb_spheres"]
    if hide is not None and len(sele) >0:
        if hide == "all":
            hide = " ".join(all_reprs)
        if not silent:
            print("(PyMol) Hiding:", sele)
        pymol.cmd.hide(representation=hide, selection=sele)

def pymol_format(representation,identifier="", hide="all", colour = None, spectrum=None, quiet=False):
    for obj in get_all_obj():
        if identifier in obj:
            pymol_hide(obj, hide=hide)
            if not quiet:
                print("(PyMol) Formatting:", obj, representation)
            pymol.cmd.show(representation=representation, selection=obj)
            if colour is not None:
                pymol_colour(colour, obj=obj, spectrum=spectrum)


def pymol_colour(colour, obj = "(all)", sele = None, spectrum=None, silent =False):
    if sele is not None:
        sele_str = "({} and {})".format(obj,sele)
    else:
        sele_str = "({})".format(obj)
    if not silent:
        print("(PyMol) Colouring:", sele_str, colour, spectrum)
    if spectrum is not None:
        pymol.cmd.spectrum(spectrum, colour, sele_str)
    else:
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
    #print(all_obj)
    obj1 = all_obj[0]
    for obj2 in all_obj:
        if obj2 != obj1:
            pymol_align__obj(obj1, obj2)
    pymol.cmd.orient("(all)")


def pymol_align__obj(obj1, obj2, orient = True, quiet=False):
    r = pymol.cmd.align(obj2, obj1)
    if not quiet:
        print3("aligned:", obj1, obj2)
    if orient:
        pymol.cmd.orient("(all)")
    return r




def pymol_align_chains(chains_to_align):
    #all_obj = pymol.cmd.get_names(type='objects')
    #print1(all_obj)
    obj1, chain1 = chains_to_align[0]
    sele1 = "{} and c. {}".format(obj1, chain1)
    for obj2, chain2 in chains_to_align[1:]:
        if obj2 != obj1:
            sele2 = "{} and c. {}".format(obj2, chain2)
            #print(sele1)
            #print(sele2)
            pymol_align__obj(sele1, sele2)

def pymol_align_chains_best(chains_to_align, double_best = False):
    #all_obj = pymol.cmd.get_names(type='objects')
    #print1(all_obj)
    obj1, chain1a, chain1b = chains_to_align[0]
    sele1a = "{} and c. {}".format(obj1, chain1a)
    sele1b = "{} and c. {}".format(obj1, chain1b)


    for obj2, chain2a, chain2b in chains_to_align[1:]:
        sele2a = "{} and c. {}".format(obj2, chain2a)
        sele2b = "{} and c. {}".format(obj2, chain2b)
        if obj2 != obj1:
            ali2a = pymol_align__obj(sele1a, sele2a, orient=False)
            print4(ali2a)
            ali2b = pymol_align__obj(sele1a, sele2b, orient=False)
            print4(ali2b)
            if ali2a[0] > ali2b[0]:
                ali1a = pymol_align__obj(sele1a, sele2a, orient=False, quiet=True), sele2a
            else:
                ali1a = ali2b, sele2b

            if double_best:
                ali2a = pymol_align__obj(sele1b, sele2a, orient=False)
                print4(ali2a)
                ali2b = pymol_align__obj(sele1b, sele2b, orient=False)
                print4(ali2b)
                if ali2a[0] < ali2b[0]:
                    ali1b = pymol_align__obj(sele1b, sele2a, orient=False, quiet=True), sele2a
                else:
                    ali1b = ali2a, sele2b

                if ali1a[0][0] < ali1b[0][0]:
                    ali = pymol_align__obj(sele1a, ali1a[1], orient=False,quiet=True)




def pymol_symmetries(obj = None):
    if obj is None:
        obj = pymol.cmd.get_names(type='objects')[0]
    print("(PyMol) Generating symmetries for: ", obj, )
    pymol.cmd.symexp("sym", obj, obj, 6)

def pymol_set_state(state):
    print("(PyMol) Setting state to: ", state, )
    pymol.cmd.set("state", state)


def pymol_orient():
    print("(PyMol) Orienting all")
    pymol.cmd.orient("(all)")

def pymol_show_cell():
    print("(PyMol) Displaying cell")
    pymol.cmd.show("cell")

def pymol_group(identifier = "sym", name = None):
    group = []
    if name is None:
        name = identifier
    print("(PyMol) Grouping:", identifier, "in", name)
    for obj in pymol.cmd.get_names(type='objects'):
        if type(identifier) is str:
            if identifier in obj:
                group.append(obj)
        if type(identifier) is list:
            for i in identifier:
                if i in obj:
                    group.append(obj)

    pymol.cmd.group(name, " ".join(group))

def pymol_draw_line(coord1, coord2, name = "d", state = -1, quiet= True):
    if not quiet:
        print("(PyMol) ({}) Line between:".format(name), coord1, "and", coord2, end="\r")
    pymol.cmd.pseudoatom("tmp1", pos=coord1, state=state)
    pymol.cmd.pseudoatom("tmp2", pos=coord2, state=state)
    pymol.cmd.distance(name, "tmp1","tmp2", state=state)
    pymol.cmd.delete("tmp1")
    pymol.cmd.delete("tmp2")

def pymol_paint_contacts(obj, contact_list, colour ="yellow"):
    print("(PyMol) Colouring contacts in {}".format(obj))
    for chain, resn, in contact_list:
        sele = "c. {} and i. {}".format(chain, resn)
        pymol_colour(colour, obj, sele, silent = True)


def pymol_temp_show(structure, disable = False):
    from Bio.PDB import PDBIO
    local["temp"] = "temp"
    exporting = PDBIO()
    exporting.set_structure(structure)
    all_obj = [n.upper() for n in pymol_get_all_objects()]
    #print(all_obj)
    name = "temp_{}.pdb".format(structure.id)
    if name.upper() in all_obj:
        n = 1
        while name.upper() in all_obj:
            name = "temp_{}_{}.pdb".format(structure.id, n)
            n += 1
    path = os.path.join(local.temp, name)
    exporting.save(path)
    pymol_start(show=True)
    if disable:
        pymol_disable("(all)")
    loaded = pymol_load_path(path)
    os.remove(path)

def pymol_get_all_objects():
    return pymol.cmd.get_names(type='objects')


def pymol_move(sele, distance = [50, 0, 0], camera = 0, state = 0):
    pymol.cmd.translate(distance, sele , state,camera)


def pymol_sphere(coords, name = None, colour="white", state = -1, scale = 8):
    if name is None:
        n_spheres = sum(1 for a in ("tmp_sphere" in obj for obj in pymol_get_all_objects()) if a)
        name = "tmp_sphere_{}".format(n_spheres)
    print(name)
    #pymol.cmd.pseudoatom(pos = coords, object = name)
    pymol.cmd.pseudoatom(object=name, pos=coords, color = colour, elem="Ca", state = state)
    pymol.cmd.set("sphere_scale", scale)#, selection="({})".format(name))
    pymol.cmd.set("sphere_transparency", 0.25)
    pymol_format("sphere", name, quiet= True)


def pymol_disable(sele = "all"):
    pymol.cmd.disable(name = sele )





def pymol_paint_all_faces(obj):
    print("(PyMOL) paining all faces in {}".format(obj.id))
    from faces import GR_colours, GR_dict
    assert obj.best_fit == "GR"
    for o in pymol_get_all_objects():
        if obj.id in o:
            for face, ress in GR_dict.items():
                sele = "({})".format(" or ".join( "i. {}".format(res) for res in ress))
                #print(sele)
                pymol_colour(GR_colours[face], o, sele, silent=True)











