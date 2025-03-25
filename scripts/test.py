
from utilities import *
import pandas as pd
import os
import setup

# Imports that need globals initialised:
from Globals import root, local, vars


from imports import pickle, export, download_pdbs, load_references, load_single_pdb


# Load/Import molecule references
sprint("Loading References")
vars["references"] = load_references()
print1("References loaded")

for reference in vars["references"]:
    if reference.name != "GR":
        continue
    print1("Generating SM for {}".format(reference.name))

    contacts_df = pd.read_csv(os.path.join(root.contacts, reference.name + ".csv"), index_col=0)
    print(contacts_df)

    sm_ssd = pd.DataFrame(columns=["dimer1", "dimer2", "index1", "index2", "similarity","diffX", "diffx"])
    n_dimers = len(contacts_df.columns) - 2
    if n_dimers < 2:
        print1("Not enough dimers in {} dataframe".format(reference.name))
        continue
    n_res = len(contacts_df)

    progress = ProgressBar(n_dimers)
    index1 = 0
    from maths import difference_between_boolean_pairs

    for id1, contacts1 in zip(contacts_df.columns, contacts_df._iter_column_arrays()):
        if id1 in ["ResNum", "ResName"]:
            continue
        progress.add(info="{}".format(id1), show_time=True)
        index1 += 1
        index2 = 0
        for id2, contacts2 in zip(contacts_df.columns, contacts_df._iter_column_arrays()):
            if id2 in ["ResNum", "ResName"]:
                continue
            index2 += 1
            #if index2 <= index1:
            #    continue
            diffX = [0, 0]
            diffx = [0, 0]
            for res in range(n_res):
                c1a, c1b = clean_list([contacts1[res]], delimiter=",", format="bool")
                c2a, c2b = clean_list([contacts2[res]], delimiter=",", format="bool")

                resX, resx = difference_between_boolean_pairs(c1a, c1b, c2a, c2b)
                diffX[0] += resX[0]
                diffX[1] += resX[1]
                diffx[0] += resx[0]
                diffx[1] += resx[1]

            dX = diffX
            dx = diffx

            if diffX[1] != 0:
                diffX = diffX[0] / diffX[1]
            else:
                diffX = 0

            if diffx[1] != 0:
                diffx = diffx[0] / diffx[1]
            else:
                diffx = 0

            similarity = max(diffX, diffx)
            # similarity = similarity / n_res
            # print(id1, id2, round(similarity,2))
            sm_ssd.loc[len(sm_ssd)] = id1, id2, index1, index2, int(similarity*100), dX, dx
        print(sm_ssd.sort_values("similarity", ascending=False).to_string())
        quit()
    print(sm_ssd)


