from ctypes import *

SOLVER_JACOBI = c_int(0)
SOLVER_SOR    = c_int(1)
SOLVER_FFT    = c_int(2)
SOLVER_MAGIC  = c_int(3)

class SolverCell(Structure):
    _fields_ = [
            ("value", c_float),
            ("value_prev", c_float),
            ("dirichlet", c_float),
            ("initial", c_float),
            ("dirichlet_present", c_char),
            ("neumann_present", c_byte * 4),
            ("neumann", c_float * 4)]


def SolverGrid_factory(length):
    class SolverGrid(Structure):
        _fields_ = [
                ("len", c_int),
                ("iters", c_int),
                ("cells", POINTER(SolverCell * length) * length)]
    g = SolverGrid()
    g.len = length
    Line = SolverCell * length
    for i in range(length):
        g.cells[i] = pointer(Line())
    return g

_solver = CDLL("./engine/solver.so")
solve_call = _solver.solve
solve_call.restype = c_double
def solve(grid, method):
    return float(solve_call(byref(grid), method))

def reduce_to_array(grid):
    arr = [[0.0 for x in range(grid.len)] for x in range(grid.len)]
    for x in range(0, grid.len):
        for y in range(0, grid.len):
            arr[x][y] = grid.cells[x].contents[y].value
    return arr
