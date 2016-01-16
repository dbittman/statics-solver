from ctypes import *

SOLVER_JACOBI = c_int(1)

class SolverCell(Structure):
    _fields_ = [
            ("value", c_float),
            ("value_prev", c_float),
            ("dirichlet", c_float),
            ("initial", c_float),
            ("error", c_float),
            ("dirichlet_present", c_char),
            ("neumann_present", c_char * 4),
            ("neumann", c_float * 4)]


def SolverGrid_factory(length):
    class SolverGrid(Structure):
        _fields_ = [
                ("len", c_int),
                ("cells", POINTER(SolverCell * length) * length)]
    g = SolverGrid()
    g.len = length
    Line = SolverCell * length
    for i in range(length):
        g.cells[i] = pointer(Line())
    return g

_solver = CDLL("./solver.so")
def solve(grid, method):
    _solver.solve(byref(grid), method)

def reduce_to_array(grid):
    arr = [[0.0 for x in range(grid.len)] for x in range(grid.len)]
    for x in range(0, grid.len):
        for y in range(0, grid.len):
            arr[x][y] = grid.cells[x].contents[y].value
    return arr
