#!/usr/bin/env python3

import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# control colour
colors = np.empty((3, 3, 3, 4))
alpha = 0.9

colors[0] = [1, 0, 0, alpha] # red
colors[1] = [0, 1, 0, alpha] # green
colors[2] = [0, 0, 1, alpha] # blue
# colors[3] = [1, 1, 0, alpha] # yellow
# colors[4] = [1, 1, 1, alpha] # grey

points = np.ones((3, 3, 3))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.voxels(points, facecolors=colors, edgecolors='grey')


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

# nx, nz = 4, 4
#
# dx = np.linspace(0, nx - 1, nx)
# hdx = np.linspace(0.5, nx - 1.5, nx - 1)
#
# bydx = np.linspace(0, nx - 1, nx)
# byhdx = np.linspace(0.5, nx - 1.5, nx - 1)
#
# Ex_x, Ex_z = np.meshgrid(hdx, dx)
# Ez_x, Ez_z = np.meshgrid(dx, hdx)
# By_x, By_z = np.meshgrid(byhdx, byhdx)
#
# wx, wz = np.meshgrid(dx[:-1], dx[:-1])
# ex_x_stencil, ex_z_stencil = np.meshgrid(hdx[:-1], dx[:-1])
# # by_x_stencil, by_z_stencil = np.meshgrid(byhdx[:-1], byhdx[:-1])
#
# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_axisbelow(True)
# ax.set_aspect('equal')
#
# ax.scatter(Ex_x, Ex_z, c='g', marker='>', s=180, label='Ex')
# ax.scatter(Ez_x, Ez_z, c='b', marker='^', s=180, label='Ez')
# # ax.scatter(By_x, By_z, c='r', marker='s', s=120, label='By')
#
# ax.scatter(wx, wz, c='orange', marker='h', s=180)
# ax.scatter([1.0], [1.0], c='m', marker='X', s=180)
# ax.scatter(ex_x_stencil, ex_z_stencil, c='tab:olive', marker='>', s=180)
# # ax.scatter(by_x_stencil, by_z_stencil, c='orange', marker='s', s=180)
#
#
# ax.set_xlim([-1, nx])
# ax.set_ylim([-1, nz])
#
# ax.set_xticks(dx)
# ax.set_xticks(hdx, minor=True)
# ax.set_yticks(dx)
# ax.set_yticks(hdx, minor=True)
#
# ax.grid(which='major', linestyle='-', linewidth=1.5, color='k')
# ax.grid(which='minor', linestyle='--', linewidth=1.0, color='gray')
#
# ax.legend()
#
# plt.show()





