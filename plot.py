from ctypes import c_void_p, cdll, string_at, c_double, c_int64
import json

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

library = cdll.LoadLibrary('./solver.so')
solver = library.solver
solver.restype = c_void_p
solved = solver(c_double(10), c_double(10), c_double(0.01), c_double(0.01), c_int64(300), c_int64(8))
solved_str = string_at(solved).decode('utf-8')
solved_array = json.loads(solved_str)
solved_mat = np.array(solved_array)
fig = plt.figure()
matplotlib.use('TkAgg')
ax = fig.add_subplot()
pos = ax.imshow(solved_mat, cmap='hot', vmin=300, vmax=3000)
fig.colorbar(pos, ax=ax)
plt.show()
