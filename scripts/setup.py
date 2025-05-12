import os, sys, platform

from scipy.ndimage import black_tophat


def setup(local_path=None, deepness_of_script =2):

    import Globals
    root = __file__
    for i in range(deepness_of_script):
        root = os.path.dirname(root)
    os.chdir(root)
    Globals.set_root(os.getcwd())
    if local_path is None:
        local_path = os.path.abspath("../localdata/projectB")
    Globals.set_local(local_path)



# Change console tab name
try:
    import inspect
    for frame in inspect.stack()[1:]:
        if frame.filename[0] != '<':
            importing_file = frame.filename
            break
    print('\33]0;{}\a'.format(os.path.basename(importing_file)), end='', flush=True)
except:
    pass


if os.name == "nt":
    setup("C:/Users/iainv/localdata/projectB")
elif "cri4" in __file__:
    setup("/localdata/iain/_local/projectB")
elif "EMERALD" in platform.node():
    setup()
elif "GARNET" in platform.node():
    setup("/mnt/d/localdata/projectB")
else:
    setup()

#### Set up essential folders ###
from Globals import root, local, vars
from utilities import *
local["molecules"] = "pickles/molecules"
local["monomers"] = "pickles/monomers"
local["dimers"] = "pickles/dimers"
local["temp"] = "temp"
vars["pymol_started"] = False
try:
    vars["tab_name"] = importing_file
except:
    pass

vars["blacklist"] = []
if "blacklist" in os.listdir(root.pdb_lists):
    with open(os.path.join(root.pdb_lists, "blacklist"), "r") as f:
        for line in f:
            vars["blacklist"].append(line)

if "force" in sys.argv or "-f" in sys.argv:
    vars["force"] = True
    supress(sys.argv.remove, "force")
    supress(sys.argv.remove, "-f")
else:
    vars["force"] = False

if "verbose" in sys.argv or "-v" in sys.argv:
    vars["verbose"] = True
    supress(sys.argv.remove, "verbose")
    supress(sys.argv.remove, "-v")
else:
    vars["verbose"] = False

if "quiet" in sys.argv or "-q" in sys.argv:
    supress(sys.argv.remove, "quiet")
    supress(sys.argv.remove, "-q")
    vars["quiet"] = True
else:
    vars["quiet"] = False


if "block" in sys.argv or "-q" in sys.argv:
    supress(sys.argv.remove, "block")
    supress(sys.argv.remove, "-b")
    vars["block"] = True
else:
    vars["block"] = False


def get_sys_vars():
    for n, arg in enumerate(sys.argv):
        if arg.startswith("--") or arg.startswith("-"):
            arg = arg.strip("-")
            value = sys.argv[n + 1]
            try:
                value = int(value)
            except:
                pass





