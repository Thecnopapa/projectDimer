
import pandas as pd
import os
import setup
import time
import sys
# Imports that need globals initialized:
from Globals import root, local, vars
from utilities import *
from imports import *

sprint("Loading References")
vars["references"] = load_references(force_reload=False)
print1("References loaded")




from faces import *

'''dimer_list = sorted(os.listdir(local.dimers))
progress = ProgressBar(len(dimer_list))
face_df = pd.DataFrame(columns=["ID", "face1", "face2", "contact_face1", "contact_face2"])
print(face_df)'''

'''for dimer_name in dimer_list:
    sprint(dimer_name)
    dimers = load_single_pdb(dimer_name, pickle_folder=local.dimers)
    for dimer in dimers:
        if dimer.best_fit != "GR":
            progress.add(info=dimer.id)
            continue
        print(dimer.face1, dimer.face2, dimer.contact_face1, dimer.contact_face2)
        face_df.loc[dimer.id] = [dimer.id, dimer.face1, dimer.face2, dimer.contact_face1, dimer.contact_face2]
        progress.add(info=dimer.id)
face_df.to_csv("/cri4/iain/projectB/dataframes/clustering/faces/GR.csv")'''

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

def find_face_discrpancies(path = "/cri4/iain/projectB/dataframes/clustering/faces/GR.csv"):
    face_df = pd.read_csv(path)
    discrepancies = {
        "double": {

    }}
    for row in face_df.itertuples():
        #print(row.face1, row.face2)
        if row.face1 == row.face2 and row.contact_face1 == row.contact_face2:
            if row.face1 != row.contact_face1:
                d =  row.face1+ " --> "+ row.contact_face1
                print("Discrepancy (DOUBLE):", d)
                if d in discrepancies["double"]:
                    discrepancies["double"][d] += 1
                else:
                    discrepancies["double"][d] = 1

    print(discrepancies)

find_face_discrpancies()