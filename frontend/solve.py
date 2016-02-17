#!/usr/bin/env python

import parse
import c_solver
import numpy as np
import matplotlib.pyplot as plt
import getopt
import sys
import math

methods = { "jacobi" : c_solver.SOLVER_JACOBI, "sor" : c_solver.SOLVER_SOR, "fft" : c_solver.SOLVER_FFT, "magic" : c_solver.SOLVER_MAGIC }

method = methods["jacobi"]

opts, args = getopt.getopt(sys.argv[1:], "m:")

if len(args) == 0:
    raise Exception("Provide an input file name to work on.")

for opt, arg in opts:
    if opt == '-m':
        try:
            method = methods[arg]
        except:
            raise Exception('Method ' + arg + ' unknown.')

grid = parse.create_grid(args[0])
if grid == 0:
    raise Exception('Input file is incomplete')

res = c_solver.solve(grid, method)
if math.isnan(res):
    sys.exit(1)
result = c_solver.reduce_to_array(grid)


fig = plt.figure(figsize=(6, 3.2))
ax = fig.add_subplot(111)
ax.set_title('colorMap')
plt.imshow(np.array(result))
ax.set_aspect('equal')

cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')


factor = complex(0, grid.len)

y, x = np.mgrid[0:100:factor, 0:100:factor]
# need to make things negative
negative_result = [[-result[x][y] for y in range(grid.len)] for x in range(grid.len)]
v, u = np.gradient(negative_result)
fig, ax = plt.subplots()
scale = 1
ax.quiver(x[::scale, ::scale], y[::scale, ::scale], u[::scale, ::scale], v[::scale, ::scale], scale=1)
plt.show()

