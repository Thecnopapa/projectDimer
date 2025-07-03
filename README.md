# ProjectDimer
A tool for protein-protein interaction classification.
<br>
<br>


**Iain Visa** @ Universitat de Barcelona / IBMB-CSIC
<br>
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

## To see data/results 
- Run `visualisation.py` with the following commands:
  * `dimer`/s + the dimer id, e.g. `1M2Z_AB'
  * `clusters-eva` to see clusters from Eva's classification
  * `clusters-cc` to see clusters from CC analysis + KMeans
  * `clusters-score`/s to see how (GR only) clusters compare between methods


## Current clustering for GR
![GR_cc_clustered.png](images/cc/GR_cc.png)

## Current Data distribution:
![GR_similarity.png](charts/Classified Similarities of GR.png)
![monomers_df.png](charts/monomers_df.png)
![failed_df.png](charts/failed_df.png)

