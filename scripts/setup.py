import os, platform

def setup(local_path=None):

    import Globals

    root = os.path.dirname(os.path.dirname(__file__))
    os.chdir(root)
    Globals.set_root(os.getcwd())
    if local_path is None:
        local_path = os.path.abspath("../localdata/projectB")
    Globals.set_local(local_path)






if os.name == "nt":
    setup("C:/Users/iainv/localdata/projectB")
elif "cri4" in __file__:
    setup("/localdata/iain/_local/projectB")
else:
    setup()
