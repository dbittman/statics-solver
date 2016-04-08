#!/usr/bin/env python

import parse
import c_solver
import numpy as np
import matplotlib.pyplot as plt
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
    fig = plt.figure(figsize=(6, 3.2))
    ax = fig.add_subplot(111)
    ax.set_title(title)
    plt.imshow(np.array(obj))
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
if do_exit:
    sys.exit(0)
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
                err += (pow(delta[y][x], 2.0)) / (pow(res, 2.0) * len(correct) ** 2)
        print(":: " + str(pow(res, 2.0)))
        err = err / (len(correct) ** 2)
        print("red. Chi-square = " + str(err))
        add_plot(correct, "Exact result")
        add_plot(delta, "Difference")


add_plot(result, "Computed result")



factor = complex(0, grid.len)
if vector:
    y, x = np.mgrid[0:100:factor, 0:100:factor]
    negative_result = [[-result[grid.len - (x+1)][y] / 10000 for y in range(grid.len)] for x in range(grid.len)]
    v, u = np.gradient(negative_result)
    fig, ax = plt.subplots()
    scale = 1
    ax.quiver(x[::scale, ::scale], y[::scale, ::scale], u[::scale, ::scale], v[::scale, ::scale], scale=1)

plt.show()

