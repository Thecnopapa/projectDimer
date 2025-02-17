import os



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
    pymol.cmd.load(os.path.join(local[folder], file_name),file_name)

def pymol_load_path(path, state = 0):
    pymol.cmd.load(path,os.path.basename(path))
    pymol.cmd.set("state", state)

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








