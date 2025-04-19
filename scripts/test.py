
import pandas as pd
import os
import setup
import time
import sys
# Imports that need globals initialised:
from Globals import root, local, vars
from utilities import *
from imports import *

sprint("Loading References")
vars["references"] = load_references(force_reload=False)
print1("References loaded")


from faces import *

objects = load_single_pdb("all", pickle_folder=local.dimers)
for dimer in objects:
    #print(dimer.face1, dimer.face2)
    print(get_dimer_faces(dimer))
    dimer.pickle()


'''objects = load_single_pdb("1M2Z", pickle_folder=local.monomers)
pcas = []
for obj in objects:
    obj.face_coms = get_face_coms(obj)
    fig, ax = plot_atoms(obj.replaced, block = False)

    for face, com in obj.face_coms.items():
        fig.tight_layout()
        com = [c-o for c, o in zip(com, obj.com)]


        ax.scatter(*com, s=500, c=GR_colours[face])
    ax.set_aspect('equal')
    plt.show(block=False)
    plt.savefig("test.png")
    quit()'''


