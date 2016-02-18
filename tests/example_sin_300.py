import math
import numpy as np
import matplotlib.pyplot as plt
import pickle

def u(x,y,l):
    return (math.sinh(((l-y)/l) * math.pi) / math.sinh(math.pi)) * math.sin(math.pi * (x/l))

l = 300
grid = [[u(x,y,l) for x in range(l)] for y in range(l)]

result = grid


with open('correct_sin_300.dat', 'wb') as f:
    pickle.dump(result, f)

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


factor = complex(0, l)

y, x = np.mgrid[0:100:factor, 0:100:factor]
# need to make things negative
negative_result = [[-result[x][y] for y in range(l)] for x in range(l)]
v, u = np.gradient(negative_result)
fig, ax = plt.subplots()
scale = 1
ax.quiver(x[::scale, ::scale], y[::scale, ::scale], u[::scale, ::scale], v[::scale, ::scale], scale=1)
plt.show()

