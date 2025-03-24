# Project_B by Iain Visa @ IBMB-CSIC / UB

# Essential imports
import os, sys
import pandas as pd

from utilities import *
import platform
import numpy as np




def main(PROCESS_ALL = False,
         SKIP_SYMMETRY = False,
         SKIP_DIMERS = False,
         FORCE_CONTACTS = False,
         LARGE_DATASET = True,
         DO_ONLY = [],
         GENERATE_SYMMETRIES = True,
         MAX_THREADS = 1,
         MINIMUM_CHAIN_LENGTH = 100,
         CONTACT_DISTANCE = 8,
         MINIMUM_CONTACTS = 0,
         FORCE_SASA = False,
         BALL_SIZE = 1.6,
         VERBOSE = False,
         ):


    ###### SET UP ######################################################################################################
    tprint("SET UP")

    vars["do_only"] = DO_ONLY
    vars["verbose"] = VERBOSE

    # Enable garbage collection
    enable_garbage_collector() # Idk if it works

    from dataframes import save_dfs, create_dfs, create_clustering_dfs
    from imports import pickle, export, download_pdbs, load_references, load_single_pdb


    # Download large dataset
    sprint("Downloading large dataset")
    if "many_pdbs" not in local.list() and LARGE_DATASET:
        download_pdbs(os.path.join(root.pdb_lists,"list_1500"), "many_pdbs", terminal = True)
    print1("Large dataset downloaded")


    # Load/Import molecule references
    sprint("Loading References")
    vars["references"] = load_references()
    print1("References loaded")


    # Create dataframes
    sprint("Creating dataframes")
    create_dfs(vars.references)
    create_clustering_dfs(vars.references)
    print1("Dataframes created")


    # Getting molecule list
    sprint("Obtaining molecule list")
    if LARGE_DATASET:
        molecule_folder = local.many_pdbs
    else:
        molecule_folder = local.experimental
    molecule_list = sorted(os.listdir(molecule_folder))
    #print(len(vars.do_only), vars.do_only)
    if len(vars.do_only) > 0:
        molecule_list = [f for f in molecule_list if any([s in f for s in vars.do_only])]
    #[print(m) for m in sorted(molecule_list)]
    print1("Molecule list obtained:", len(molecule_list), "molecules")

    eprint("SET UP")
    ###### SYMMETRY & DIMER GENERATION #################################################################################
    tprint("SYMMETRY & DIMER GENERATION")

    if not SKIP_SYMMETRY and not PROCESS_ALL:
        progress = ProgressBar(len(molecule_list))
        for m in molecule_list:
            if "lock" in m:
                sprint(".lock file detected:", m)
                continue
            filename = m.split(".")[0]
            sprint(filename)
            molecules = load_single_pdb(filename, local.molecules, molecule_folder, force_reload=PROCESS_ALL)
            for molecule in molecules:
                if GENERATE_SYMMETRIES:
                    molecule.get_all_dimers(force=PROCESS_ALL,
                                            minimum_chain_length=MINIMUM_CHAIN_LENGTH,
                                            contact_distance=CONTACT_DISTANCE,
                                            min_contacts=MINIMUM_CONTACTS,
                                            )
                else:
                    molecule.get_dimers()
                molecule.pickle()
                progress.add(info=molecule.id)
        save_dfs(clustering=False)


    eprint("SYMMETRY & DIMER GENERATION")
    ###### DIMER ANALYSIS ##############################################################################################
    tprint("DIMER ANALYSIS")

    if not SKIP_DIMERS and not PROCESS_ALL:
        print(list(vars.clustering["contacts"].keys()))
        progress = ProgressBar(len(molecule_list))
        from surface import build_contact_arrays
        for m in molecule_list:
            if "lock" in m:
                sprint(".lock file detected:", m)
                continue
            filename = m.split(".")[0]
            sprint(filename)
            molecules = load_single_pdb(filename, local.molecules)
            for molecule in molecules:
                dimers = molecule.dimers
                for dimer in dimers:
                    print1(dimer)
                    if dimer.incomplete:
                        continue
                    build_contact_arrays(dimer, sasa=False, force= FORCE_CONTACTS)
                    dimer.pickle()
                progress.add(info=molecule.id)
        for reference in vars.references:
            if reference.name+".csv" in root.contacts:
                reference.contacts_df = pd.read_csv(root.contacts[reference.name+".csv"])
        save_dfs(general=False)

    eprint("DIMER ANALYSIS")
    ###### CLUSTERING ##################################################################################################
    tprint("CLUSTERING")


    from clustering import compare_contacts, cluster, cc_analysis
    for reference in vars.references:
        sprint(reference.name)
        if reference.name == "GR":
            compare_contacts(reference)
        cc_path = cc_analysis(reference)
        #cluster(reference)
        save_dfs(general=False)



    eprint("CLUSTERING")
    ###### SAVE & EXIT #################################################################################################
    tprint("SAVE & EXIT")

    # Save and exit
    #save_dfs() # Save dataframes and generate figures

    eprint("DONE")
########################################################################################################################


if __name__ == "__main__":
    # Setup paths and globals
    print(sys.argv)
    DO_ONLY = []
    if "force" in sys.argv or "-f" in sys.argv:
        PROCESS_ALL = True
        supress(sys.argv.remove, "force")
        supress(sys.argv.remove, "-f")
    else:
        PROCESS_ALL = False
    if "verbose" in sys.argv or "-v" in sys.argv:
        VERBOSE = True
        supress(sys.argv.remove, "verbose")
        supress(sys.argv.remove, "-v")
    else:
        VERBOSE = False

    print(sys.argv)
    if len(sys.argv) >= 2 and not "all" in sys.argv:
        DO_ONLY = [arg.upper() for arg in sys.argv[1:]]
    import setup

    # Imports that need globals initialised:
    from Globals import root, local, vars

    main(PROCESS_ALL=PROCESS_ALL, # Ignore saved pickles and generate everything from scratch
         SKIP_SYMMETRY = True,
         SKIP_DIMERS = False,
         FORCE_CONTACTS = False,
         LARGE_DATASET = True, # Use a large dataset (delete all local data previously to avoid errors)
         DO_ONLY = DO_ONLY, # ( list of strings / string) Names of pdbs to be processed (CAPS sensitive, separated by space) e.g ["5N10", "1M2Z"] or "5N10 1M2Z"
         GENERATE_SYMMETRIES=True,
         MAX_THREADS=1, # Number of threads, might not be implemented yet, (0 or 1 deactivate threading)
         MINIMUM_CHAIN_LENGTH=100, # Minimum number of residues to consider a chain for dimerization (to ignore ligands and small molecules)
         CONTACT_DISTANCE = 8, # Minimum (less or equal than) distance in Angstroms to consider a contact between atoms
         MINIMUM_CONTACTS = 0, # Minimum number of contacts to consider a dimer interface
         FORCE_SASA = True,
         BALL_SIZE = 1.6,
         VERBOSE = VERBOSE
         )

    quit()

    from surface import surface
    surface(FORCE_SASA=True,
            FORCE_SIMILARITY=True
            )

    from clustering import clustering
    clustering(FORCE_ALL=False,
               )



