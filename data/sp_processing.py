#!/usr/bin/env python3

import numpy as np
from scipy import constants
import matplotlib.pyplot as plt

# t2 = np.load('/home/cepheid/TriForce/game_engine/data/test_0000000001.npy')
# print(t2.shape, t2.dtype)
# print(t2)
# gamma = 1 / np.sqrt(1 - ((v * v).sum(axis=1)) / constants.c**2)
#
# ke = (weight * (gamma - 1) * constants.m_e * constants.c**2).sum()
#
# print(f'v = {v}')
# print(f'Gamma = {gamma}')
# print(f'KE = {ke} J')

# xmin, xmax = -0.5, 0.5
# zmin, zmax = -0.5, 0.5
#
# xs = np.linspace(xmin, xmax, 100)
# zs = np.linspace(zmin, zmax, 100)
#
# X, Z = np.meshgrid(xs, zs)
#
# Er = -1.109 / (X**2 + Z**2)
#
# theta = np.arctan2(Z, X)
#
# Ex = Er * np.cos(theta)
# Ez = Er * np.sin(theta)
#
# fig, ax = plt.subplots()
# im = ax.contourf(X, Z, np.sqrt(Ex**2 + Ez**2))
# fig.colorbar(im, ax=ax)
# plt.show()

# nx, nz = 3, 3
#
# dx = np.linspace(0, 3, nx + 1)
# hdx = np.linspace(0.5, 2.5, nx)
#
# bydx = np.linspace(0, 3, nx + 1)
# byhdx = np.linspace(0.5, 1.5, nx)
#
# Ex_x, Ex_z = np.meshgrid(dx, hdx)
# Ez_x, Ez_z = np.meshgrid(hdx, dx)
# By_x, By_z = np.meshgrid(bydx, bydx)
#
# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_axisbelow(True)
# ax.set_aspect('equal')
#
# ax.scatter(Ex_x, Ex_z, c='g', marker='>', s=180, label='Ex')
# ax.scatter(Ez_x, Ez_z, c='b', marker='^', s=180, label='Ez')
# ax.scatter(By_x, By_z, c='r', marker='s', s=120, label='By')
#
# ax.set_xlim([-0.5, 3.5])
# ax.set_ylim([-0.5, 3.5])
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





