import os
import globals



globals.set_root("../")
if os.name == "nt":
    globals.set_local("C:/Users/iainv/localdata/_local/projectB")
from globals import root, local



from  molecules import *
molecule1 = PDB(os.path.join(root.references, "AR.pdb"))

references = []

for reference in os.listdir(root.references):
    if reference.endswith(".pdb"):
        references.append(Reference(os.path.join(root.references,reference)))
