#!/usr/bin/env python3

import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


vertices = np.array([[0, 0, 0],
                     [1, 0, 0],
                     [1, 1, 0],
                     [0, 1, 0],
                     [0, 0, 1],
                     [1, 0, 1],
                     [1, 1, 1],
                     [0, 1, 1]])

edges = [[0, 1], [1, 2], [2, 3], [3, 0],  # Bottom face
         [4, 5], [5, 6], [6, 7], [7, 4],  # Top face
         [0, 4], [1, 5], [2, 6], [3, 7]]  # Connecting edges




fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

c1 = np.copy(vertices)
# c2 = np.copy(vertices)
# c2[:, 0] += 1
# c3 = np.copy(vertices)
# c3[:, 1] += 1
# c4 = np.copy(vertices)
# c4[:, 0] += 1
# c4[:, 1] += 1
# c5 = np.copy(c1)
# c5[:, 2] += 1
# c6 = np.copy(c2)
# c6[:, 2] += 1
# c7 = np.copy(c3)
# c7[:, 2] += 1
# c8 = np.copy(c4)
# c8[:, 2] += 1

for e in edges:
    ax.plot(*zip(*c1[e]), color='k', alpha=0.25)
    # ax.plot(*zip(*c2[e]), color='k', alpha=0.25)
    # ax.plot(*zip(*c3[e]), color='k', alpha=0.25)
    # ax.plot(*zip(*c4[e]), color='k', alpha=0.25)
    # ax.plot(*zip(*c5[e]), color='k', alpha=0.25)
    # ax.plot(*zip(*c6[e]), color='k', alpha=0.25)
    # ax.plot(*zip(*c7[e]), color='k', alpha=0.25)
    # ax.plot(*zip(*c8[e]), color='k', alpha=0.25)


nx, nz = 2, 2

dx = dy = dz = np.linspace(0, nx - 1, nx)
hdx = hdy = hdz = np.linspace(0.5, nx - 1.5, nx - 1)

px, py, pz = np.meshgrid(dx, dy, dz)

Ex_x, Ex_y, Ex_z = np.meshgrid(hdx, dy, dz)
Ey_x, Ey_y, Ey_z = np.meshgrid(dx, hdy, dz)
Ez_x, Ez_y, Ez_z = np.meshgrid(dx, dy, hdz)
Bx_x, Bx_y, Bx_z = np.meshgrid(dx, hdy, hdz)
By_x, By_y, By_z = np.meshgrid(hdx, dy, hdz)
Bz_x, Bz_y, Bz_z = np.meshgrid(hdx, hdy, dz)

ax.scatter(Ex_x, Ex_y, Ex_z, c='g', marker='$E_x$', s=180)
ax.scatter(Ey_x, Ey_y, Ey_z, c='b', marker='$E_y$', s=180)
ax.scatter(Ez_x, Ez_y, Ez_z, c='c', marker='$E_z$', s=180)

ax.scatter(Bx_x, Bx_y, Bx_z, c='r', marker='$B_x$', s=180)
ax.scatter(By_x, By_y, By_z, c='y', marker='$B_y$', s=180)
ax.scatter(Bz_x, Bz_y, Bz_z, c='m', marker='$B_z$', s=180)

# wsx, wsy, wsz = np.meshgrid(dx, dy, dz)
# ax.scatter(wsx, wsy, wsz, c='tab:purple', marker='h', s=180)

# particle = np.array([1.25, 0.5, 0.5])
# ax.scatter(particle[0], particle[1], particle[2], c='tab:orange', marker='X', s=180)

# xs = [[0.5, 0.5], [1.5, 1.5]]
# ys = [[0, 2]]
# zs = [[0, 0], [1, 1], [2, 2]]

# for lx in xs:
#     for ly in ys:
#         for lz in xs:
#             ax.plot(lx, ly, lz, c='k', ls='--', alpha=0.25)
#             ax.plot(ly, lx, lz, c='k', ls='--', alpha=0.25)
#             ax.plot(lx, lz, ly, c='k', ls='--', alpha=0.25)


plt.show()





