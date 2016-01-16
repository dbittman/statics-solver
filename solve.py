import parse
import c_solver
import numpy as np
import matplotlib.pyplot as plt

grid = parse.create_grid("example_sin.txt")
if grid == 0:
    raise Exception('Input file is incomplete')

c_solver.solve(grid, c_solver.SOLVER_JACOBI)

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

