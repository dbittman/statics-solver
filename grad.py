#!/usr/bin/env python2
import numpy as np
import matplotlib.pyplot as plt

w = 60
h = 60

arr = [[0 for x in range(w)] for x in range(h)]

for line in open("out"):
    a = line.split()
    if(len(a) > 0):
        arr[int(a[0]) // 10][int(a[1]) // 10] = -float(a[2])

y, x = np.mgrid[0:100:60j, 0:100:60j]
z = x * np.exp(-x**2 - y**2)
v, u = np.gradient(arr)
fig, ax = plt.subplots()
ax.quiver(x, y, u, v, scale=100)
plt.show()

