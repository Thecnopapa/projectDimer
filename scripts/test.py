
from utilities import *
import pandas as pd
import os
import setup
import time
import sys
# Imports that need globals initialised:
from Globals import root, local, vars


progress = ProgressBar(10)
for i in range(10):
    progress.add()
    time.sleep(0.5)

