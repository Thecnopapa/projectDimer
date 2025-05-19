
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

dimer_list = sorted(os.listdir(local.dimers))
progress = ProgressBar(len(dimer_list))
face_df = pd.DataFrame(columns=["ID", "face1", "face2", "contact_face1", "contact_face2"])
#print(dimer_list)

n_dimers = 0
thresholds = range(1,11)
ref_name = "GR"
contact_maps = {threshold:None for threshold in thresholds}
for dimer_name in dimer_list:
    sprint(dimer_name)
    #print(local.dimers)
    dimers = load_single_pdb(dimer_name, pickle_folder=local.dimers)
    for dimer in dimers:
        if dimer.best_fit == ref_name and dimer.best_fit != "Missmatch":
            if "contact_surface" not in dimer.__dict__ or False:
                dimer.contact_surface = ContactSurface(dimer.monomer1.replaced, dimer.monomer2.replaced)
                dimer.pickle()
            if dimer.contact_surface is not None:
                for threshold, contact_map in contact_maps.items():
                    if contact_map is None:
                        contact_maps[threshold] = dimer.contact_surface.get_contact_map(threshold=threshold)
                    else:
                        contact_maps[threshold] = np.add(contact_map, dimer.contact_surface.get_contact_map(threshold=threshold))
                n_dimers += 1
        print(contact_maps)
        progress.add(info=dimer.id)

for threshold, contact_map in contact_maps.items():
    ContactSurface.get_heat_map(contact_map, title = "{} heat-map, threshold = {}, N = {}".format(ref_name, threshold, n_dimers), normalize = n_dimers)








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

#find_face_discrpancies()