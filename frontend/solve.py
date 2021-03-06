#!/usr/bin/env python

import parse
import c_solver
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import getopt
import sys
import math
import pickle

methods = { "jacobi" : c_solver.SOLVER_JACOBI, "sor" : c_solver.SOLVER_SOR, "fft" : c_solver.SOLVER_FFT, "magic" : c_solver.SOLVER_MAGIC }

method = methods["jacobi"]
verifier = ""
opts, args = getopt.getopt(sys.argv[1:], "m:v:Vp")
vector = False
do_exit = False
def add_plot(obj, title):
    my_cmap = matplotlib.cm.get_cmap('magma')
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    ax.set_title(title)
    x = np.arange(0, len(obj[0]), 1)
    y = np.arange(0, len(obj[0]), 1)
    X, Y = np.meshgrid(x, y)
    levels = np.arange(-800.0,800,100)
    CS = plt.contour(X, Y, obj, levels)
    plt.imshow(np.array(obj), cmap=my_cmap)
    plt.colorbar(orientation='vertical')
    ax.set_aspect('equal')


if len(args) == 0:
    raise Exception("Provide an input file name to work on.")

for opt, arg in opts:
    if opt == '-m':
        try:
            method = methods[arg]
        except:
            raise Exception('Method ' + arg + ' unknown.')
    if opt == '-v':
        verifier = arg
    if opt == '-V':
        vector = True
    if opt == '-p':
        do_exit = True

grid = parse.create_grid(args[0])
if grid == 0:
    raise Exception('Input file is incomplete')

res = c_solver.solve(grid, method)
if math.isnan(res):
    sys.exit(1)
result = c_solver.reduce_to_array(grid)
with open('out.dat', 'wb') as f:
    pickle.dump(result, f)

if verifier != "":
    with open(verifier, 'rb') as f:
        print("Verifying using " + verifier)
        correct = pickle.load(f)
        # quick sanity check
        assert(len(correct) == len(result))
        assert(len(correct[0]) == len(result[0]))
        delta = [[correct[y][x] - result[y][x] for x in range(len(correct[0]))] for y in range(len(correct))]
        err = 0.0
        for y in range(len(correct)):
            for x in range(len(correct[0])):
                if delta[y][x] != 0.0:
                    err += pow(delta[y][x], 2.0)
        err = math.sqrt(err) / len(correct)
        print("err per cell = " + str(err))
        if not do_exit:
            add_plot(correct, "Exact result")
            add_plot(delta, "Difference")

if do_exit:
    sys.exit(0)

add_plot(result, "Computed result - " + args[0])



factor = complex(0, grid.len)
if vector:
    y, x = np.mgrid[0:100:factor, 0:100:factor]
    negative_result = [[-result[grid.len - (x+1)][y] / 5000 for y in range(grid.len)] for x in range(grid.len)]
    v, u = np.gradient(negative_result)
    fig, ax = plt.subplots()
    scale = 2
    ax.quiver(x[::scale, ::scale], y[::scale, ::scale], u[::scale, ::scale], v[::scale, ::scale], scale=2)

plt.show()

