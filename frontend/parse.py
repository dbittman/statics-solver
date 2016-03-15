from ctypes import *
import c_solver

import math

def parse_coords(grid, token):
    coords = token.split(',')
    return (int(eval(coords[0])), int(eval(coords[1])))

def parse_cell(grid, tokens):
    x,y = parse_coords(grid, tokens[1])
    print(" * Set cell " + str(x) + "," + str(y) + " initial to " + tokens[3])
    grid.initials[x-1][y-1] = float(tokens[3])

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
        grid.dirichlets[x-1][y-1] = eval(expr)
        grid.dirichlet_presents[x-1][y-1] = 1

directions = {
        'top':0,
        'right':1,
        'bottom':2,
        'left':3
        }

def parse_neumann(grid, tokens):
    # WARNING - off-by-one error central in here!
    sx,sy = parse_coords(grid, tokens[1])
    ex,ey = parse_coords(grid, tokens[2])
    dist = int(math.sqrt((ex - sx) ** 2 + (ey - sy) ** 2))
    ix = abs((ex - sx) / dist)
    iy = abs((ey - sy) / dist)
    expr = ' '.join(str(tok) for tok in tokens[5:])
    direc = directions[tokens[3]]
    print(" * Setting Neumann " + str(sx) + "," + str(sy) + " through " + str(ex) + "," + str(ey) + " along " + tokens[3] + " to " + expr)
    for x,y in set([(int(sx + inc*ix), int(sy + inc*iy)) for inc in range(0, dist+1)]):
        grid.neumanns[direc].contents[x-1][y-1] = eval(expr)
        grid.neumann_presents[x-1][y-1] |= (1 << direc)
        print("NEUPR = " + str(grid.neumann_presents[x-1][y-1]))

        mirrorx = x-1
        mirrory = y-1
        if direc == 0:
            mirrory -= 1
        elif direc == 1:
            mirrorx += 1
        elif direc == 2:
            mirrory += 1
        else:
            mirrorx -= 1
        if mirrorx >= 0 and mirrorx < grid.len and mirrory >= 0 and mirrory < grid.len:
            grid.neumanns[(direc + 2) % 4].contents[mirrorx][mirrory] = eval(expr)
            grid.neumann_presents[mirrorx][mirrory] |= (1 << ((direc+2) % 4))


parsers = {
        'cell':parse_cell,
        'dirichlet':parse_dirichlet,
        'neumann':parse_neumann
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

