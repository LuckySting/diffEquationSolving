import ctypes
import json

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

library = ctypes.cdll.LoadLibrary('./solver.so')
solver = library.solver
solver.restype = ctypes.c_void_p
solved = solver()
solved_str = ctypes.string_at(solved).decode('utf-8')
solved_array = json.loads(solved_str)
solved_mat = np.array(solved_array)
fig = plt.figure()
matplotlib.use('TkAgg')
ax = fig.add_subplot()
pos = ax.imshow(solved_mat, cmap='hot', vmin=300, vmax=3000)
fig.colorbar(pos, ax=ax)
plt.show()
