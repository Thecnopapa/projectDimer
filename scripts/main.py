# Project_B by Iain Visa @ IBMB-CSIC / UB

# Essential imports
import os, sys
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

from utilities import *
import platform
import numpy as np




def main(PROCESS_ALL = False,
         SKIP_SYMMETRY = False,
         SKIP_DIMERS = False,
         SKIP_CLUSTERING = False,
         FORCE_CONTACTS = False,
         FORCE_COMPARE = True,
         ONLY_GR = False,
         FORCE_CLUSTERING = False,
         LARGE_DATASET = True,
         DO_ONLY = [],
         GENERATE_SYMMETRIES = True,
         MAX_THREADS = 1,
         MINIMUM_CHAIN_LENGTH = 100,
         CONTACT_DISTANCE_SYMMETRY = 8,
         CONTACT_DISTANCE_CLUSTERING = 16,
         MINIMUM_CONTACTS = 0,
         SASA = False,
         FORCE_SASA = False,
         BALL_SIZE = 1.6,
         VERBOSE = False,
         QUIET = False,
         REPROCESS_DIMERS = False,
         SPLIT_FACES = True,
         FORCE_SPLIT = False,
         CLUSTER_BY_PCA = True,
         N_CLUSTERS = 4,
         DIMENSIONS_PCA = [0,1,2],
         FACES_BY_COM = True,
         MINIMUM_SCORE = 0,
         COMPARE = True,
         SPLIT_FACES_ANYWAY = False,
         REMOVE_REDUNDANCY = True,
         CLUSTERING_METHOD = "KMeans",
         N_SAMPLE_MULTIPLIER = 0.5,
         QUANTILE = 0.1,
         BANDWIDTH = 0.1,
         HEATMAPS =True,
         GIFS = True,
         SNAPSHOTS = True,
         CHAINBOWS = False,
         GENERATE_CLUSTERS = True,
         DELETE_PREVIOUS = True,
         REFRESH_PLOTS = True,
         REPROCESS_CLUSTERS = False,
         USE_FACES = "eva",
         PLOT_DIHEDRALS = False,
         DIHEDRAL_ALGORITHM  = "MeanShift",
         FACE_ALGORITHM = "affinity",
         ):


    ###### SET UP ######################################################################################################
    tprint("SET UP")

    vars["do_only"] = DO_ONLY
    #vars["verbose"] = VERBOSE
    #vars["quiet"] = QUIET

    # Enable garbage collection
    enable_garbage_collector() # Idk if it works

    from dataframes import save_dfs, create_dfs, create_clustering_dfs
    from imports import pickle, export, download_pdbs, load_references, load_single_pdb, load_list_1by1

    sprint("Blacklist:", vars.blacklist)

    # Download large dataset
    sprint("Downloading large dataset")
    if "many_pdbs" not in local.list() and LARGE_DATASET:
        pdb_file = "rcsb_pdb_ids_20250524034011.txt"
        download_pdbs(os.path.join(root.pdb_lists, pdb_file), "many_pdbs", terminal = True)
    print1("Large dataset downloaded")


    # Load/Import molecule references
    sprint("Loading References")
    vars["references"] = load_references(force_reload=PROCESS_ALL)
    vars.log("References loaded:", *[vars.references])
    print1("References loaded")


    # Create dataframes
    sprint("Creating dataframes")
    create_dfs(vars.references)
    create_clustering_dfs(vars.references)
    for reference in vars.references:
        reference.restore_reference_dfs(reset=PROCESS_ALL)
        reference.pickle()

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
    print1("Molecule list obtained:", len(molecule_list), "molecules, DO_ONLY = ", DO_ONLY, log=True)

    eprint("SET UP")
    ###### SYMMETRY & DIMER GENERATION #################################################################################
    tprint("SYMMETRY & DIMER GENERATION")

    if not SKIP_SYMMETRY or PROCESS_ALL or len(os.listdir(local.molecules)) == 0:
        progress = ProgressBar(len(molecule_list))
        for m in molecule_list:
            if "lock" in m:
                sprint(".lock file detected:", m)
                progress.add(info=m)
                continue
            filename = m.split(".")[0]
            sprint(filename)
            molecules = load_single_pdb(filename, local.molecules, molecule_folder, force_reload=PROCESS_ALL)
            for molecule in molecules:
                molecule.get_all_dimers(force=PROCESS_ALL,
                                        minimum_chain_length=MINIMUM_CHAIN_LENGTH,
                                        contact_distance=CONTACT_DISTANCE_SYMMETRY,
                                        min_contacts=MINIMUM_CONTACTS,
                                        )
                molecule.pickle()
            progress.add(info=m)
        save_dfs()


    eprint("SYMMETRY & DIMER GENERATION")
    ###### DIMER ANALYSIS ##############################################################################################
    tprint("DIMER ANALYSIS v2")

    if not SKIP_DIMERS and not PROCESS_ALL:
            for dimer in load_list_1by1(pickle_folder=local.dimers):
                dimer.reprocess(contacts=False, faces=False)



    eprint("DIMER ANALYSIS v2")



    """tprint("DIMER ANALYSIS")

    if not SKIP_DIMERS or REPROCESS_DIMERS or PROCESS_ALL or (SPLIT_FACES_ANYWAY and FORCE_SPLIT):
        #print(list(vars.clustering["contacts"].keys()))
        progress = ProgressBar(len(molecule_list))
        from surface import build_contact_arrays
        c_arrays = {ref.name: [] for ref in vars.references}
        for m in molecule_list:
            if "lock" in m:
                sprint(".lock file detected:", m)
                progress.add(info=m)
                continue
            filename = m.split(".")[0]
            sprint(filename)
            molecules = load_single_pdb(filename, local.molecules)
            for molecule in molecules:
                dimers = molecule.dimers
                for dimer in dimers:
                    print1(dimer)
                    if REPROCESS_DIMERS:
                        dimer.process()
                    if dimer.incomplete:
                        continue
                    dimer.get_contacts(max_distance=CONTACT_DISTANCE_CLUSTERING, force= FORCE_CONTACTS)
                    if dimer.best_fit == "GR":
                        dimer.get_faces(by_com = FACES_BY_COM)
                        face_df = vars["clustering"]["faces"][dimer.best_fit]
                        face_df.loc[dimer.id] = [dimer.id, dimer.face1, dimer.face2, dimer.contact_face1, dimer.contact_face2]
                    build_contact_arrays(dimer, c_arrays, sasa=SASA, force=FORCE_CONTACTS or PROCESS_ALL, max_contact_length=CONTACT_DISTANCE_CLUSTERING)
                    dimer.pickle()

            progress.add(info=m)

        for key in c_arrays.keys():
            vars.clustering["contacts"][key] = pd.concat([vars.clustering["contacts"][key].iloc[:,0:2], *c_arrays[key]], axis=1)
            print(vars.clustering["contacts"][key])
        save_dfs(general=False, clustering=True)

        for reference in vars.references:
            sprint("Faces and contact dfs for {}".format(reference.name))
            reference.faces_df = vars.clustering["faces"][reference.name]
            reference.contacts_df = vars.clustering["contacts"][reference.name]
            reference.pickle()
            print(reference.faces_df)
            print(reference.contacts_df)



    eprint("DIMER ANALYSIS")"""
    ###### CLUSTERING ##################################################################################################
    tprint("CLUSTERING v2")

    if not SKIP_CLUSTERING or PROCESS_ALL and False:
        from clustering import generate_dihedrals_df, cluster_angles, create_clusters, cluster_redundancy, cluster_dihedrals, \
            get_faces, compare_all_with_eva
        from imports import load_clusters


        if GENERATE_CLUSTERS or len(os.listdir(local.cluster_pickles))<=3:
            if DELETE_PREVIOUS:
                if os.path.exists(local.cluster_pickles):
                    print("Deleting files in: {}".format(local.cluster_pickles))
                    for file in sorted(os.listdir(local.cluster_pickles)):
                        if ONLY_GR and "GR" not in file:
                            continue
                        print(file, end="\r")
                        os.remove(os.path.join(local.cluster_pickles, file))
            generate_dihedrals_df(force=False or PROCESS_ALL)



            for file in sorted(os.listdir(root.dihedrals)):
                sprint("Clustering1")
                if ONLY_GR and "GR" not in file:
                    continue
                dihedrals_path = os.path.join(root.dihedrals, file)
                cluster1_folder = cluster_angles(dihedrals_path,
                                                 bandwidth=40,
                                                 angles=["a0", "a1", "a2"],
                                                 cluster_name="angle_cluster1",
                                                 folder="angle_clusters1")



            for file in sorted(os.listdir(cluster1_folder)):
                sprint("Clustering2")
                if ONLY_GR and "GR" not in file:
                    continue
                dihedrals_path = os.path.join(cluster1_folder, file)
                cluster2_folder = cluster_angles(dihedrals_path,
                                                 bandwidth=30,
                                                 angles=["b0", "b1", "b2"],
                                                 cluster_name="angle_cluster2",
                                                 folder="angle_clusters2",
                                                 split_by="angle_cluster1",
                                                 save_together=True,
                                                 )



            if SNAPSHOTS:
                from pyMol import pymol_start
                pymol_start(show=False)


            for n, file in enumerate(sorted(os.listdir(cluster2_folder))):
                if "--1" in file:
                    continue
                if ONLY_GR and "GR" not in file:
                    continue
                ref_name = file.split(".")[0]
                ref = [ref for ref in vars.references if ref.name == ref_name][0]
                sprint(ref_name+ "({}/{})".format(n, len(os.listdir(cluster2_folder))))
                dihedrals_path = os.path.join(cluster2_folder, file)
                create_clusters(dihedrals_path, ref )

            cluster_redundancy()
            cluster_dihedrals()

        if ONLY_GR:
            identifier = "GR"
        else:
            identifier = "ALL"

        # ONLY FOR ALL- ALL
        get_faces(algorithm = FACE_ALGORITHM, force=True, gif=GIFS)
        quit()
        compare_all_with_eva()

        for cluster in load_clusters(identifier=identifier, onebyone=True):
                if cluster.is_all:
                    continue
                sprint(cluster.id)
                cluster.process_cluster(force=REPROCESS_CLUSTERS or USE_FACES not in cluster.faces.keys(),
                                        faces = True,
                                        use_face=USE_FACES,
                                        dihedral_algorithm=DIHEDRAL_ALGORITHM)
                cluster.pickle()

        cluster_colours = None
        if PLOT_DIHEDRALS:
            cluster_colours = "dihedrals"


        for cluster in load_clusters(identifier=identifier, onebyone=True):
            sprint(cluster.id)
            cluster.plot_cluster(force = REFRESH_PLOTS, plot=True,
                                 snapshot = SNAPSHOTS,
                                 face_colours = USE_FACES,
                                 cluster_colours=cluster_colours)
            cluster.pickle()


        from clustering import generate_cluster_grids, get_space_groups
        for ref in vars.references:
            if ref.name != "GR" and ONLY_GR:
                continue
            generate_cluster_grids(identifier=ref.name, use_faces=USE_FACES, face_algorithm=FACE_ALGORITHM)
            get_space_groups(identifier=ref.name, use_faces=USE_FACES)

        if GIFS:
            for cluster in load_clusters(identifier=identifier, onebyone=True):
                sprint(cluster.id)
                cluster.plot_cluster(plot=True, gif =GIFS, cluster_colours=cluster_colours)
                cluster.pickle()


        """matrix, oneDmatrix1, oneDmatrix2 = plot_dihedrals(dihedrals_path,
                                                            clusters="angle_cluster2",
                                                            subset_col="angle_cluster2",
                                                            subset = None,
                                                            heatmap = HEATMAPS, hm_threshold=10,
                                                            outer_ids_complete=ref.get_outer_res_list(complete_list=True),
                                                            gif=GIFS,
                                                            snapshot=SNAPSHOTS,
                                                            chainbows = CHAINBOWS,
                                                            include_all=True,)"""




        '''from clustering import compare_contacts, get_clusters, cluster, split_by_faces, cluster_by_face
        for reference in vars.references:
            sprint(reference.name)
            if reference.name != "GR" and ONLY_GR:
                reference.pickle()
                continue
            if SPLIT_FACES or SPLIT_FACES_ANYWAY:
                split_by_faces(reference, force=FORCE_SPLIT, by_com= FACES_BY_COM)
            if reference.name == "GR" and (COMPARE or PROCESS_ALL):
                if not "classified_df" in reference.__dict__ or (FORCE_COMPARE or PROCESS_ALL):
                    reference.classified_df = compare_contacts(reference, force = FORCE_COMPARE or PROCESS_ALL)
                from clustering import add_info_to_classified
                #save_dfs(general=False, clustering=True)
                add_info_to_classified(reference)
                from visualisation import classified_chart
                classified_chart()
                #reference.clusters_eva = get_clusters(reference.classified_df, column = "Best_Match", ref_name=reference.name)
    
    
            cluster_by_face(reference, FORCE_ALL= FORCE_CLUSTERING or PROCESS_ALL, minimum_score=MINIMUM_SCORE,
                            n_clusters=N_CLUSTERS, pca=CLUSTER_BY_PCA, pca_dimensions = DIMENSIONS_PCA,
                            splitted=SPLIT_FACES, rem_red=REMOVE_REDUNDANCY, method = CLUSTERING_METHOD,
                            quantile=QUANTILE, n_sample_multiplier=N_SAMPLE_MULTIPLIER, bandwidth = BANDWIDTH)
            reference.pickle()
        #save_dfs(general=False, clustering=True)'''



    eprint("CLUSTERING v2")
    ###### SAVE & EXIT #################################################################################################
    tprint("SAVE & EXIT")

    # Save and exit
    #save_dfs() # Save dataframes and generate figures
    ring_bell(times=3)

    eprint("DONE")
########################################################################################################################


if __name__ == "__main__":

    import setup

    # Setup paths and globals
    print(sys.argv)
    DO_ONLY = []


    print(sys.argv)
    if len(sys.argv) > 1 and not "all" in sys.argv:
        DO_ONLY = [arg.upper() for arg in sys.argv[2:]]


    # Imports that need globals initialized:
    from Globals import root, local, vars
    from utilities import *

    main(
        PROCESS_ALL=vars.force, # Master switch

        # Setup and data import
        VERBOSE=vars.verbose,
        QUIET=vars.quiet,
        DO_ONLY=DO_ONLY, # ( list of strings / string) Names of PDBs to be processed (CAPS sensitive?, separated by space) e.g ["5N10", "1M2Z"] or "5N10 1M2Z"
        LARGE_DATASET=True,  # Use a large dataset (delete all local data previously to avoid errors)
        MAX_THREADS=1,  # Number of threads, might not be implemented yet, (0 or 1 deactivate threading)

        # Symmetry calculations, and generation of Monomers + Dimers
        SKIP_SYMMETRY = True, # Skip the entire block (overridden by PROCESS_ALL)
        MINIMUM_CHAIN_LENGTH=100,# Minimum number of residues to consider a chain for dimerization (to ignore ligands and small molecules)
        CONTACT_DISTANCE_SYMMETRY=8,  # Minimum (less or equal than) distance in Angstroms to consider a contact between atoms
        MINIMUM_CONTACTS=0,  # Minimum number of contacts to consider a dimer interface

        # Dimer processing, includes contact calculation and face identification, generates contact dataframes
        SKIP_DIMERS = True, # Skip the entire block (overridden by PROCESS_ALL and REPROCESS_DIMERS)
        REPROCESS_DIMERS = True,
        FORCE_CONTACTS = False,  # Force contact calculation if already calculated (overridden by PROCESS_ALL)
        CONTACT_DISTANCE_CLUSTERING = 12,
        FACES_BY_COM = True,

        # SASA related (BROKEN)
        SASA = False, # Whether to run SASA calculations, currently broken
        FORCE_SASA=True, # DEPRECATED
        BALL_SIZE=1.6, # DEPRECATED


        # Compare GR clustering to EVA clustering
        COMPARE = True, # requiered
        FORCE_COMPARE= True,

        # Split by faces based on Eva
        SPLIT_FACES=False,
        SPLIT_FACES_ANYWAY = False,
        FORCE_SPLIT=True,

        # Clustering, from SM to plotting
        SKIP_CLUSTERING=False and "clustering" not in sys.argv, # Skip th entire block (overridden by PROCESS_ALL)
        FORCE_CLUSTERING=True,  # Force clustering if already calculated (overridden by PROCESS_ALL)
        ONLY_GR = False, # Whether to only cluster GR
        REMOVE_REDUNDANCY = True,
        CLUSTERING_METHOD = "MeanShift",
        QUANTILE= 0.1,
        N_SAMPLE_MULTIPLIER = None,
        BANDWIDTH = 0.03,

        N_CLUSTERS = 4,
        CLUSTER_BY_PCA = True,
        DIMENSIONS_PCA = [0,1,2],
        MINIMUM_SCORE = 0,

        HEATMAPS = True,
        GIFS = False,
        SNAPSHOTS = True,
        CHAINBOWS = False,
        GENERATE_CLUSTERS = False or "clusters" in sys.argv,
        DELETE_PREVIOUS = False or "delete" in sys.argv,
        REFRESH_PLOTS = False or "replot" in sys.argv,
        REPROCESS_CLUSTERS = False or "reprocess" in sys.argv,
        FACE_ALGORITHM = "weighted",
        USE_FACES = "generated", # "eva" or "generated"
        PLOT_DIHEDRALS = False or "dihedrals" in sys.argv,
        DIHEDRAL_ALGORITHM  = "HDBSCAN",


        )

    #quit()


