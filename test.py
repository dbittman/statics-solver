#!/usr/bin/env python
from ctypes import *

solver = CDLL("./solver.so")

class SolverCell(Structure):
    _fields_ = [
            ("value", c_float),
            ("neumann", c_float * 4),
            ("value_prev", c_float),
            ("dirichlet", c_float),
            ("initial", c_float),
            ("error", c_float),
            ("char", c_char * 4),
            ("char", c_char)]


def SolverGrid_factory(length):
    class SolverGrid(Structure):
        _fields_ = [
                ("len", c_int),
                ("cells", POINTER(SolverCell * length) * length)]
    g = SolverGrid()
    g.len = length
    return g


Line = SolverCell * 10
grid = SolverGrid_factory(10)


for i in range(10):
    l = Line()
    for j in range(10):
        l[j].value = j * i
    grid.cells[i] = pointer(l)

solver.python_call(byref(grid))

print(grid.cells[1].contents[1].value)

