
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

for self in vars["references"]:
    self.best_fit = self.name
    self.replaced = self.structure
    self.com = find_com(self.structure.get_atoms())
    if self.name == "GR":
        self.face_coms = get_face_coms(self)
    self.pca = dict(pca=get_pca(self.structure, com=self.com),
                    com=self.com,
                    chain=self.chain)
    self.pickle()

quit()
dimer_list = os.listdir(local.dimers)
progress = ProgressBar(len(dimer_list))
for dimer_name in dimer_list:
    sprint(dimer_name)
    dimers = load_single_pdb(dimer_name, pickle_folder=local.dimers)
    for dimer in dimers:
        #print(dimer.face1, dimer.face2)
        if dimer.best_fit != "GR":
            progress.add(info=dimer.id)
            continue
        dimer.face1, dimer.face2, dimer.interface_distance = get_dimer_faces(dimer)
        dimer.pickle()
        progress.add(info=dimer.id)

def plot_face_coms():
    objects = load_single_pdb("1M2Z", pickle_folder=local.monomers)
    pcas = []
    for obj in objects:
        obj.face_coms = get_face_coms(obj)
        fig, ax = plot_atoms(obj.replaced, block = False)

        for face, com in obj.face_coms.items():
            fig.tight_layout()
            com = [c-o for c, o in zip(com, obj.com)]


            ax.scatter(*com, s=500, c=GR_colours[face])
        ax.set_aspect('equal')
        plt.show(block=True)
        plt.savefig("test.png")


