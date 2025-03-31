
from utilities import *
import pandas as pd
import os
import setup
import time
import sys
# Imports that need globals initialised:
from Globals import root, local, vars

from imports import *

sprint("Loading References")
vars["references"] = load_references(force_reload=False)
print1("References loaded")


from clustering import split_by_faces
for reference in vars.references:
    if reference.name != "GR":
        continue
    split_by_faces(reference)