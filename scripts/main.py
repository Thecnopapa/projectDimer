import os
import globals



globals.set_root("../")
if os.name == "nt":
    globals.set_local("C:/Users/iainv/localdata/_local/projectB")
from globals import root, local




from references import load_references
references = load_references(force_reload=True)

from  molecules import *












# Save and exit
for reference in references:
    reference.pickle()