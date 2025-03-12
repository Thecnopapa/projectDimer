# Project_B by Iain Visa @ IBMB-CSIC / UB

# Essential imports
import os, sys
import pandas as pd
from utilities import *
import platform




def main(PROCESS_ALL = False,
         LARGE_DATASET = True,
         DO_ONLY = [],
         GENERATE_SYMMETRIES = True,
         MAX_THREADS = 1,
         ):


    # Set up
    vars["do_only"] = DO_ONLY
    print(local.list())

    # Enable garbage collection
    enable_garbage_collector()

    from dataframes import save_dfs, create_dfs
    from imports import pickle, export


    # Download large dataset
    tprint("Downloading large data")
    from imports import download_pdbs
    if "many_pdbs" not in local.list() and LARGE_DATASET:
        download_pdbs(os.path.join(root.pdb_lists,"list_1500"), "many_pdbs", terminal = True)
    eprint("Large data downloaded")


    # Load/Import molecule references
    tprint("Loading files")
    from imports import load_references, load_single_pdb
    from molecules import Reference, PDB
    vars["references"] = load_references()


    # Create dataframes
    sprint("Creating dataframes")
    create_dfs(vars.references)
    eprint("Dataframes created")


    if LARGE_DATASET:
        molecule_folder = local.many_pdbs
    else:
        molecule_folder = local.experimental
    molecule_list = os.listdir(molecule_folder)
    #print(len(vars.do_only), vars.do_only)
    if len(vars.do_only) > 0:
        molecule_list = [f for f in molecule_list if any([s in f for s in vars.do_only])]
    #[print(m) for m in sorted(molecule_list)]
    progress = ProgressBar(len(molecule_list))
    for m in sorted(molecule_list):
        if "lock" in m:
            sprint(".lock file detected:", m)
            continue

        filename = m.split(".")[0]
        sprint(filename)
        molecules = load_single_pdb(filename, local.molecules, local.many_pdbs, force_reload=PROCESS_ALL)
        for molecule in molecules:
            if GENERATE_SYMMETRIES:
                molecule.get_all_dimers()
            else:
                molecule.get_dimers()
            molecule.pickle()
            progress.add(info=molecule.id)
        save_dfs()



    # Save and exit
    tprint("Saving data")
    save_dfs() # Save dataframes and generate figures



if __name__ == "__main__":
    # Setup paths and globals
    print(sys.argv)
    DO_ONLY = []

    if "force" in sys.argv:
        PROCESS_ALL = True
        sys.argv.remove("force")
    else:
        PROCESS_ALL = False

    if len(sys.argv) > 2:
        DO_ONLY = sys.argv[2:]

    import setup

    # Imports that need globals initialised:
    from Globals import root, local, vars

    main(PROCESS_ALL=PROCESS_ALL, # Ignore saved pickles and generate everything from scratch
         LARGE_DATASET = True, # Use a large dataset (delete all local data previously to avoid errors)
         DO_ONLY = DO_ONLY, # ( list of strings / string) Names of pdbs to be processed (CAPS sensitive, separated by space) e.g ["5N10", "1M2Z"] or "5N10 1M2Z"
         )
    quit()

    from surface import surface
    surface(FORCE_SASA=True,
            FORCE_SIMILARITY=True
            )

    from clustering import clustering
    clustering(FORCE_ALL=False,
               )



