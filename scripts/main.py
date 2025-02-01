import os
import globals





globals.set_root("../")
print(globals.root)


print(os.listdir(globals.root.path))

print(os.path.abspath(globals.root.path))

from globals import root

