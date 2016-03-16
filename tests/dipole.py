import numpy as np
import matplotlib.pyplot as plt
import getopt
import sys
import math
import pickle

DIPOLE_DIST = 10
DIPOLE_X = 105
DIPOLE_Y = 100

with open('out.dat', 'rb') as f:
    output = pickle.load(f)
    for i in range(0, 20):
        print(str((i + .1)) + " " + str(output[DIPOLE_X+i][DIPOLE_Y]) + " .03")
