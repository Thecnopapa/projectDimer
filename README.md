# ProjectDimer
A tool for protein-protein interaction classification.
<br>
<br>

**Iain Visa** @ Universitat de Barcelona / IBMB-CSIC
<br>
<br>

## Installation: 
In a Linux `bash` terminal (Windows not fully tested): 
- Clone this repository wherever you like.
- Download and install [conda](https://www.anaconda.com).
- Follow instructions in `environment.txt` to setup the environment.
- Source `bashrc.project` by adding the following snippet at the end of your `.bashrc` file:
  ```
  
  FOLDER_PATH = "/path/to/.../projectDimer"
  projectDimer(){
      if [ -f $FOLDER_PATH/bashrc.project ]; then
          source $FOLDER_PATH/bashrc.project
      
      fi
  }
  ```
- Restart the shell (`bash` or `exec bash`).
- Type `projectDimer`.
- You should see printed:
  ```
  projectDimer path is: /path/to/.../projectDimer
  ```
- `run` and `show` commands should now be available.


## Configuration
To set up your desired dataset:
- Save a`comma`separated list of PDB entries (can be downloaded from a query in the PDB) to the `data/pdb_lsits` folder.
- Tweak all desired parameters in the config.txt file (WIP).
- Set up file system with `run setup`.
- To manually add PDB files for analysis add them to the `localdata/projectDimer/many_pdbs` folder if not changed.
- Execute main script with `run main`

## Visualising results
The main generated plots and images can be found in `localdata/projectDimer/`. The following file structure shows were to fain some of these:
```

.
├──projectDimer (repository)
│   ├── charts
│   └── dataframes
├──localdata
.   └── projectDimer
        ├── images
        │   ├── clusters_by_face (clusters PyMol snapshots, by their interacting regions)
        │   ├── dihedral_figs (angles of clusters)
        │   ├── heatmaps (interaction surface matrixes)
        │   └── res_cords (mpl: interaction surface 3D representation)
        └── snapshots
            ├── _faces_generated (PyMol snapshots of each cluster, default)
            └── ... (other requested cluster snapshots)

```
The `show` command allows visualisation of the main elements of interest.
Each element has a unique ID within its category, for example:
- `molecule`: 1M2Z
- `ref`: GR
- `monomer`: 1M2Z_A, 1M2Z_A_2_000 
- `dimer`: 1M2Z_Aa_2_000
- `cluster` GR-02-00

The general command format is:  
`show` `[category]` `ID` `option1` `option2` `...`  


Extra `options` to tweak the display are:
- `pymol` to display th entity with PyMol (with `cluster`, must add `session` to keep the interactive session open)
- `mpl` shows `matplotlib` representation (if available)
- `plot` might display some related plots
- `c=[colour]` to colour PyMol representations (if available), `colour` = any PyMol colour, or `matrix`, or `mutations` 


## Tweaking results
Results can be partially regenerated with new parameters.  
To do so modify the config.txt file with the desired parameters (WIP).  
And then execute the `run main` command with the any of the following keywords:
- `force` to redo EVERYTHING.
- `dimers` to reprocess dimers
- `clusters` to regenerate clusters
- `reprocess` to reanalyse clusters
- `replot` to regenerate all figures
- `delete` to delete all previously generated clusters, and regenerate them














