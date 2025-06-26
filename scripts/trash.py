

def cluster_by_face(reference, FORCE_ALL=False, DIMENSIONS=3, n_clusters = 4, minimum_score=0, pca = True,
                    pca_dimensions = [0,1,2], splitted=True, rem_red = True, method = "KMeans", quantile=0.1, n_sample_multiplier = 0.5,
                    bandwidth=0.1):



    if FORCE_ALL:
        FORCE_SM = True
        FORCE_CC = True
        FORCE_CLUSTER = True
        FORCE_PLOT = True
    else:
        FORCE_SM = False
        FORCE_CC = False
        FORCE_CLUSTER = False
        FORCE_PLOT = True

    if splitted:
        subfolder_name = "{}_" + reference.name
    else:
        subfolder_name = "{}"

    #print(root[subfolder_name.format("contacts")])

    for file in os.listdir(root[subfolder_name.format("contacts")]):

        if "filtered" in file:
            continue
        if not splitted:
            if "GR" not in file:
                continue
            if "contacts" in file:
                continue

        sprint(file)
        contacts_path = os.path.join(root[subfolder_name.format("contacts")], file)
        #print(contacts_path)
        if minimum_score > 0:
            contacts_df = pd.read_csv(contacts_path)
            original_len = len(contacts_df)
            classified_df = pd.read_csv(os.path.join(root.classified, "{}.csv".format(reference.name)))
            print(classified_df)
            cols = list(contacts_df.columns)
            for col in cols:
                if col in ["ResNum", "ResName"]:
                    continue
                #print(col)
                #print(classified_df[classified_df["ID"] == col])
                class_row = classified_df[classified_df["ID"] == col]
                #print(class_row)
                if len(class_row) < 1 :
                    cols.remove(col)
                    continue
                class_row = class_row.iloc[0]
                #print(class_row.Similarity)
                if class_row.Similarity < minimum_score:
                    cols.remove(col)
            filtered_len = len(cols)
            #print(contacts_df[cols])
            print1("Filtered {} dimers with a threshold of {}".format(original_len-filtered_len, minimum_score))
            contacts_path = os.path.join(root[subfolder_name.format("contacts")], "filtered_{}_".format(minimum_score) + file)
            contacts_df[cols].to_csv(contacts_path)



        if not pca:
            sms_path= generate_sm(reference, force=FORCE_SM, subfolder = subfolder_name, in_path = contacts_path)
            ccs_path = cc_analysis(reference, force=FORCE_CC, dimensions=DIMENSIONS, subfolder = subfolder_name, in_path = sms_path)
            if ccs_path is None:
                print("CC analysis failed")
                return
            clustered_path = clusterize_cc(reference, force=FORCE_CLUSTER, dimensions=DIMENSIONS, n_clusters=n_clusters,
                                           subfolder = subfolder_name, in_path = ccs_path)
            plot_path = plot_cc(reference, labels=False, labels_centres=True, # force = FORCE_PLOT
                                          dimensions=DIMENSIONS, subfolder = subfolder_name, in_path = clustered_path)
        else:
            from faces import get_pca_df, plot_pcas
            pca_path = get_pca_df(in_path=contacts_path, subfolder=subfolder_name, force=FORCE_CC, splitted=splitted)
            if rem_red:
                pca_path = remove_redundancy(pca_path)
            #plot_pcas(pcas, title="GR: {}  (N = {})".format(file.split(".")[0], len(pcas)))
            clustered_path = clusterize_pcas(method=method, subfolder=subfolder_name, in_path = pca_path, force=FORCE_CLUSTER,
                                             dimensions=pca_dimensions,splitted=splitted, quantile =quantile, n_sample_multiplier=n_sample_multiplier, bandwidth=bandwidth)
            plot_path = plot_clustered_pcas(reference, labels=False, labels_centres=True,  force = FORCE_PLOT,
                                dimensions=DIMENSIONS, subfolder=subfolder_name, in_path=clustered_path, pca = True,
                                            pca_dimensions=pca_dimensions, splitted=splitted)
    add_clusters_to_classified(reference, pca = pca, splitted=splitted)
    return

def cluster(reference, FORCE_ALL=False, DIMENSIONS = 3, score_id = "", thread = False):
    if FORCE_ALL:
        FORCE_SM = True
        FORCE_CC = True
        FORCE_CLUSTER = True
        FORCE_PLOT = True
    else:
        FORCE_SM = False
        FORCE_CC = False
        FORCE_CLUSTER = False
        FORCE_PLOT = False


    reference.sm_path = generate_sm(reference, force=FORCE_SM)
    reference.cc_path = cc_analysis(reference, force=FORCE_CC, dimensions=DIMENSIONS)
    reference.clustered_path = clusterize_cc(reference, force=FORCE_CLUSTER, dimensions=DIMENSIONS)
    reference.plot_path = plot_cc(reference, labels=False, labels_centres=True, force=FORCE_PLOT, dimensions=DIMENSIONS)
    return
    if reference.name == "GR":
        tprint("Comparing to Eva")
        df_cc = pd.read_csv(os.path.join(root.clustered, "GR_cc.csv"))
        scores = calculate_scores_GR(df_cc, score_id+str(DIMENSIONS))
        print1("Scores: cc: {}, eva: {}".format(scores[0], scores[1]))
        eprint("Compared successfully")




def cluster_snapshot_old(file, clusters, levels=None, color_clusters=False, chainbows = True):
    from imports import load_single_pdb, load_references
    from pyMol import pymol_start, pymol_load_path, pymol_colour,pymol_list_to_bfactors, pymol_align_chains, pymol_group, \
        pymol_open_saved_cluster, pymol_get_all_objects, pymol_save_temp_session, pymol_save_cluster, pymol_open_session_terminal, \
        colours,ncolours, pymol_reset, pymol_orient, pymol_save_snapshot, get_all_obj, pymol_disable, pymol_delete, \
        pymol_command_in_new_process, pymol_reinitialize
    sprint("Cluster snampshot")

    cluster_folders = ["angle_clusters2"]
    cluster_cols = ["angle_cluster2"]
    if levels is not None:
        l = levels -1
        cluster_folders = cluster_folders[l:]
        cluster_cols = cluster_cols[:l]


    options = []
    sele = []
    fname = file.split(".")[0]
    for n, (cluster_col, cluster_folder) in enumerate(zip(cluster_cols, cluster_folders)):
        dihedrals_path = os.path.join(root[cluster_folders[n]], fname + ".csv")
        df = pd.read_csv(file, index_col=0)
        #print(df)
        options.append([int(a) for a in sorted(set(df[cluster_cols[n]].values))])
        s = clusters[n]
        if s == "all":
            sele.append(options[n])
        else:
            sele.append([s])
        df.query(" | ".join(["{} == {}".format(cluster_col, n) for n in sele]), inplace=True)
        #print(df.to_string())
        fname = fname + "-{}".format(s)

    #pymol_start(show=False)
    print(file)
    filename = os.path.basename(fname)
    print(filename)
    ref = load_references(identifier=filename.split("-")[0])[0]
    pymol_load_path(ref.path, ref.name)
    pymol_colour("chainbow", ref.name)
    print("Sele:", sele)

    for c in sele[-1]:
        if c == -1 or c == "-1":
            continue
        if c != 1: #DeBUG
            continue
        print("##",get_all_obj())


        print("Cluster:", sele[:-1], c)
        subset = df[df[cluster_cols[-1]] == c]
        print(subset)
        chains_to_align = [[ref.name, ref.chain]]
        for row in subset.itertuples():
            dimer = load_single_pdb(identifier=row.id, pickle_folder=local.dimers, quiet=True)[0]
            name = pymol_load_path(dimer.replaced_path, row.id + str(row.is1to2))
            if row.is1to2:
                chains_to_align.append([name, row.mon1])
            else:
                chains_to_align.append([name, row.mon2])
            if chainbows:
                pymol_colour("chainbow", name)
            elif color_clusters:
                c = int(row.__getattribute__(cluster_cols[-1]))
                pymol_colour(colours[c % ncolours], name)
                pymol_colour("chainbow", name)
            else:
                resids = [res.id[1] for res in dimer.monomer1.replaced.get_residues()]
                sele1 = name + " and c. {}".format(dimer.monomer1.chain)
                sele2 = name + " and c. {}".format(dimer.monomer2.chain)
                list1 = [min(x) for x in dimer.contact_surface.d_s_matrix.T]
                list2 = [min(x) for x in dimer.contact_surface.d_s_matrix]
                pymol_list_to_bfactors(val_list=list1, obj_name=sele1, resids=resids)
                pymol_list_to_bfactors(val_list=list2, obj_name=sele2, resids=resids)
                pymol_colour("blue_yellow_red", name, spectrum="b")
        print(chains_to_align)
        print(get_all_obj())

        pymol_align_chains(chains_to_align)
        pymol_orient()
        local["snapshots"] = "snapshots"
        #session_path = pymol_save_temp_session()
        #pymol_open_session_terminal(session_path)
        pymol_save_snapshot(filename + "-{}".format(c), folder=local.snapshots)
        print("hi")
        for obj, _ in chains_to_align[1:]:
            pymol_delete(obj)
    pymol_delete()
    collect_garbage()
    pymol_reinitialize()
    collect_garbage()
    import time
    time.sleep(10)



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

