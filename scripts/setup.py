import os, platform

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
else:
    setup()

#### Set up essential folders ###
from Globals import root, local, vars
local["molecules"] = "pickles/molecules"
local["temp"] = "temp"
vars["pymol_started"] = False
try:
    vars["tab_name"] = importing_file
except:
    pass



