import os
from utilities import *
from globals import root, local, vars
import numpy as np
import pandas as pd

import Bio.PDB
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, StructureBuilder, Structure, SASA


def get_contact_res(self=None, monomer1=None, monomer2=None, dimer_structure = None, radius=1.4, n_points=100):

    sr = SASA.ShrakeRupley(n_points=n_points, probe_radius=radius)

    if self is not None:
        monomer1 = self.monomer1
        monomer2 = self.monomer2
    if dimer_structure is None:
        dimer_structure = self.replaced_structure

    structure1 = monomer1.structure
    structure2 = monomer2.structure

    sr.compute(structure1,level="R")
    sr.compute(structure2,level="R")
    dimer_sasa = sr.compute(dimer_structure,level="R")

    sasa_df = pd.DataFrame(columns=["ID", "name", "sasa1", "sasa2", "sasa1_dimer", "sasa2_dimer", "is_contact"])

    if monomer1.id == monomer2.id:
        for i, res in enumerate(monomer1.get_residues()):
            sasa_df[len(sasa_df)] = [res.id, res.resname, structure1[i].sasa, structure2[i].sasa, dimer_structure[0][i].sasa, dimer_structure[1][i].sasa]
    if self is not None:
        self.sasa_df = sasa_df
    return sasa_df


if __name__ == "__main__":

    import globals

    globals.set_root("../")
    if os.name == "nt":
        globals.set_local("C:/Users/iainv/localdata/_local/projectB")
    elif os.name == "posix":
        globals.set_local("/localdata/iain/_local/projectB")
    from globals import root, local, vars

    tprint("Loading dimers")
    from imports import load_dimers
    dimers = load_dimers()
    eprint("Dimers loaded")

    tprint("Calculating SASAs")
    progress = ProgressBar(len(dimers))
    for dimer in dimers:
        get_contact_res(dimer)
        progress.add(info = dimer.id)
    eprint("SASAs calculated")



    from imports import pickle
    print("Pickling...")
    pickle(dimers)
    print("Done")