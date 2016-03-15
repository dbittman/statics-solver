from ctypes import *

SOLVER_JACOBI = c_int(0)
SOLVER_SOR    = c_int(1)
SOLVER_FFT    = c_int(2)
SOLVER_MAGIC  = c_int(3)

_solver = CDLL("./engine/solver.so")
def SolverGrid_factory(length):
    class SolverGrid(Structure):
        _fields_ = [
                ("len", c_int),
                ("iters", c_int),
                ("values", POINTER(POINTER(c_float))),
                ("value_prevs", POINTER(POINTER(c_float))),
                ("initials", POINTER(POINTER(c_float))),
                ("dirichlet_presents", POINTER(POINTER(c_byte))),
                ("dirichlets", POINTER(POINTER(c_float))),
                ("neumann_presents", POINTER(POINTER(c_byte))),
                ("neumanns", (POINTER(POINTER(c_float)) * 4))
                ]
    g = SolverGrid()
    g.len = length
    _solver.init_grid(byref(g))
    return g

solve_call = _solver.solve
solve_call.restype = c_double
def solve(grid, method):
    return float(solve_call(byref(grid), method))

def reduce_to_array(grid):
    arr = [[0.0 for x in range(grid.len)] for x in range(grid.len)]
    for x in range(0, grid.len):
        for y in range(0, grid.len):
            arr[x][y] = grid.values[x][y]
    return arr
