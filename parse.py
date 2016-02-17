from ctypes import *
import c_solver

import math

def parse_coords(grid, token):
    coords = token.split(',')
    return (int(eval(coords[0])), int(eval(coords[1])))

def parse_cell(grid, tokens):
    x,y = parse_coords(grid, tokens[1])
    print(" * Set cell " + str(x) + "," + str(y) + " initial to " + tokens[3])
    grid.cells[x-1].contents[y-1].initial = float(tokens[3])

def parse_dirichlet(grid, tokens):
    # WARNING - off-by-one error central in here!
    sx,sy = parse_coords(grid, tokens[1])
    ex,ey = parse_coords(grid, tokens[2])
    dist = int(math.sqrt((ex - sx) ** 2 + (ey - sy) ** 2))
    ix = abs((ex - sx) / dist)
    iy = abs((ey - sy) / dist)
    expr = ' '.join(str(tok) for tok in tokens[4:])
    print(" * Setting Dirichlet " + str(sx) + "," + str(sy) + " through " + str(ex) + "," + str(ey) + " to " + expr)
    for x,y in set([(int(sx + inc*ix), int(sy + inc*iy)) for inc in range(0, dist+1)]):
        grid.cells[x-1].contents[y-1].dirichlet = eval(expr)
        grid.cells[x-1].contents[y-1].dirichlet_present = 1

parsers = {
        'cell':parse_cell,
        'dirichlet':parse_dirichlet
        }

def create_grid(input_path):
    grid = 0
    print("Building grid...")
    with open(input_path) as infile:
        for _, line in enumerate(infile):
            tokens = line.split('#', 1)[0].split()
            if len(tokens) == 0:
                continue
            # okay, let's parse shit.
            if tokens[0] == 'gridsize':
                size = int(tokens[1])
                print(" * Creating grid of size " + str(size) + "x" + str(size))
                grid = c_solver.SolverGrid_factory(size)
            else:
                if grid == 0:
                    raise Exception('must specify gridsize before any other command')
                parsers[tokens[0]](grid, tokens)
    return grid

