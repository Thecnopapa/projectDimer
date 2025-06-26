
import pandas as pd
import os
import setup
import time
import sys
# Imports that need globals initialized:
from Globals import root, local, vars
from utilities import *
from imports import *
from faces import *
from maths import *


sprint("Loading References")
vars["references"] = load_references(force_reload=True)
print1("References loaded")
for ref in vars.references:
    ref.pickle()

quit()
from clustering import cluster_redundancy
cluster_redundancy()




quit()

df1 = pd.DataFrame([[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]])
df2 = pd.DataFrame([[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]])
print(df1)
print(df2)
print("#")
print(pd.concat([df1,df2], axis=0))


quit()

from ARDB import *

parse_ardb()
ardb_sequence = parse_ardb_sequence()
print(ardb_sequence)
for ref in vars["references"]:
    print(ref)
    print(ref.sequence)
    if ref.name == "AR":
        from alignments import *
        al = get_alignment_map(ref.sequence,ardb_sequence)
        [print("{} : {}".format(k,v)) for k,v in al.items()]









quit()








from clustering import get_faces, compare_all_with_eva, cluster_redundancy, get_space_groups, generate_cluster_grids


for cluster in load_clusters("ALL", onebyone=True):
    tprint(cluster.id)
    cluster.cluster_dihedrals(method="HDBSCAN")
    cluster.pickle()

quit()
generate_cluster_grids(identifier="ER")
get_space_groups(identifier="ER")

quit()


list_path = os.path.join(root.pdb_lists, "rcsb_pdb_ids_20250524034011.txt")

l = []
with open(list_path, "r") as f:
    for line in f:
        l.extend(line.split(","))
print(l)
list_to_table(l, ncols = 10, path=os.path.join(root.pdb_lists, "rcsb_pdb_ids_20250524034011-70x10.csv"))






quit()
from clustering import get_faces, compare_all_with_eva, cluster_redundancy, get_space_groups
#get_faces(force=False)
#compare_all_with_eva()
#quit()

from clustering import generate_cluster_grids, get_space_groups

#face_combinations = generate_cluster_grids(identifier="GR", use_faces="generated")
get_space_groups(identifier="GR", use_faces="generated")






quit()
for ref in vars.references:
    sprint(ref.id)
    print1(ref.header)



dimer_list = sorted(os.listdir(local.dimers))
progress = ProgressBar(len(dimer_list))
face_df = pd.DataFrame(columns=["ID", "face1", "face2", "contact_face1", "contact_face2"])
#print(dimer_list)

n_dimers = 0
thresholds = range(1,11)
ref_name = "GR"
ref = [r for r in vars.references if r.name == ref_name][0]
contact_maps = {threshold:None for threshold in thresholds}
for dimer_name in dimer_list:
    sprint(dimer_name)
    #print(local.dimers)
    dimers = load_single_pdb(dimer_name, pickle_folder=local.dimers, quiet=True)
    for dimer in dimers:
        if dimer.best_fit == ref_name and dimer.best_fit != "Missmatch":
            if "contact_surface" not in dimer.__dict__ or True:
                dimer.contact_surface = ContactSurface(dimer.monomer1.replaced, dimer.monomer2.replaced, outer_ids=ref.outer_ids)
                dimer.pickle()
            if dimer.contact_surface is not None:# and len(list(dimer.monomer1.replaced.get_residues())) == 229:
                for threshold, contact_map in contact_maps.items():
                    if contact_map is None:
                        contact_maps[threshold] = dimer.contact_surface.get_contact_map(threshold=threshold)
                    else:
                        contact_maps[threshold] = np.add(contact_map, dimer.contact_surface.get_contact_map(threshold=threshold))
                n_dimers += 1
        print(contact_maps)
        progress.add(info=dimer.id)




for threshold, contact_map in contact_maps.items():
    if len(contact_map) == 0:
        continue
    from pyMol import pymol_start, pymol_save_temp_session, pymol_open_session_terminal

    pymol_start(show=False)
    title = "{} heat-map, threshold = {}, N = {}".format(ref_name, threshold, n_dimers)

    ContactSurface.display_heatmap(matrix=contact_map, title=title,
                                   structure=[ref.structure for ref in vars.references if ref.name == ref_name][0],
                                   n_samples=n_dimers,
                                   show_pymol=True,
                                   obj_name = "{}-T{}-N{}".format(ref_name, threshold, n_dimers),
                                   show_heatmap=False,
                                   #colors=["blue", "yellow", "red", "red"],
                                   #cvals=[0, 0.25, 0.5, 1],
                                   percentage=True
                                   )
session = pymol_save_temp_session(name="heatmap_plot_{}-N{}".format(ref_name,n_dimers) + ".pse")
pymol_open_session_terminal(session)














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