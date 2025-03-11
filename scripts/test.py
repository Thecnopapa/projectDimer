

import gc

gc.enable()

print(len(gc.get_objects()))

gc.collect()
input()
print(len(gc.get_objects()))
