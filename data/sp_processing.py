#!/usr/bin/env python3

import numpy as np
# from scipy import constants
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from time import time

Nx = 5
Ny = 10
Nz = 15

# ncells = np.zeros((Nx - 1, Ny - 1, Nz - 1))
# full = np.zeros((Nx, Ny, Nz))
ex = np.zeros((Nx - 1, Ny, Nz))
ey = np.zeros((Nx, Ny - 1, Nz))
ez = np.zeros((Nx, Ny, Nz - 1))
hx = np.zeros((Nx, Ny - 1, Nz - 1))
hy = np.zeros((Nx - 1, Ny, Nz - 1))
hz = np.zeros((Nx - 1, Ny - 1, Nz))
ex_update = ex[:, 1:-1, 1:-1]
ey_update = ey[1:-1, :, 1:-1]
ez_update = ez[1:-1, 1:-1, :]

eydz = ey[:, :, 1:]
ezdy = ez[:, 1:, :]

ezdx = ez[1:, :, :]
exdz = ex[:, :, 1:]

eydx = ey[1:, :, :]
exdy = ex[:, 1:, :]

hydx = hy[1:, 1:-1, :]
hxdy = hx[1:-1, 1:, :]

hzdx = hz[1:, :, 1:-1]
hxdz = hx[1:-1, :, 1:]

hzdy = hz[:, 1:, 1:-1]
hydz = hy[:, 1:-1, 1:]


strides_map = {
    1: '1zu',
    (Nz - 1): 'Nz - 1',
    Nz: 'Nz',
    (Ny - 1) * (Nz - 1): '(Ny - 1) * (Nz - 1)',
    (Ny - 1) *  Nz     : '(Ny - 1) * Nz',
     Ny      * (Nz - 1): 'Ny * (Nz - 1)',
     Ny      *  Nz     : 'Ny * Nz',
}

shape_map = {
    Nx: 'Nx',
    Nx - 1: 'Nx - 1',
    Nx - 2: 'Nx - 2',
    Ny: 'Ny',
    Ny - 1: 'Ny - 1',
    Ny - 2: 'Ny - 2',
    Nz: 'Nz',
    Nz - 1: 'Nz - 1',
    Nz - 2: 'Nz - 2'
}

def map_strides(stride):
    stride = [i // 8 for i in stride]
    return f'{strides_map[stride[0]]}, {strides_map[stride[1]]}, {strides_map[stride[2]]}'

def map_shape(shape):
    return f'{shape_map[shape[0]]}, {shape_map[shape[1]]}, {shape_map[shape[2]]}'


print(f'Ez: shape: {map_shape(ez.shape)}, strides: {map_strides(ez.strides)}')
print(f'Hx: shape: {map_shape(hx.shape)}, strides: {map_strides(hx.strides)}')
print(f'Hy: shape: {map_shape(hy.shape)}, strides: {map_strides(hy.strides)}')
print(f'Ez_update: shape: {map_shape(ez_update.shape)}, strides: {map_strides(ez_update.strides)}')

print(f'Eydz: shape: {map_shape(eydz.shape)}, strides: {map_strides(eydz.strides)}')
print(f'Ezdy: shape: {map_shape(ezdy.shape)}, strides: {map_strides(ezdy.strides)}')
print(f'Ezdx: shape: {map_shape(ezdx.shape)}, strides: {map_strides(ezdx.strides)}')
print(f'Exdz: shape: {map_shape(exdz.shape)}, strides: {map_strides(exdz.strides)}')

print(f'Hxdy: shape: {map_shape(hxdy.shape)}, strides: {map_strides(hxdy.strides)}')
print(f'Hydx: shape: {map_shape(hydx.shape)}, strides: {map_strides(hydx.strides)}')

# eydx_shape = eydx.shape
# ezdx_shape = ezdx.shape
# eydx_strides = eydx.strides
# ezdx_strides = ezdx.strides
#
# exdy_shape = exdy.shape
# ezdy_shape = ezdy.shape
# exdy_strides = exdy.strides
# ezdy_strides = ezdy.strides
#
# exdz_shape = exdz.shape
# eydz_shape = eydz.shape
# exdz_strides = exdz.strides
# eydz_strides = eydz.strides
#
# hydx_shape = hydx.shape
# hzdx_shape = hzdx.shape
# hydx_strides = hydx.strides
# hzdx_strides = hzdx.strides
#
# hxdy_shape = hxdy.shape
# hzdy_shape = hzdy.shape
# hxdy_strides = hxdy.strides
# hzdy_strides = hzdy.strides
#
# hxdz_shape = hxdz.shape
# hydz_shape = hydz.shape
# hxdz_strides = hxdz.strides
# hydz_strides = hydz.strides
#
# ncells_shape = ncells.shape
# full_shape = full.shape
# ex_shape = ex.shape
# ey_shape = ey.shape
# ez_shape = ez.shape
# hx_shape = hx.shape
# hy_shape = hy.shape
# hz_shape = hz.shape
#
# ex_update_shape = ex_update.shape
# ey_update_shape = ey_update.shape
# ez_update_shape = ez_update.shape
#
# ncells_strides = ncells.strides
# full_strides = full.strides
# ex_strides = ex.strides
# ey_strides = ey.strides
# ez_strides = ez.strides
# hx_strides = hx.strides
# hy_strides = hy.strides
# hz_strides = hz.strides
#
# ex_update_strides = ex_update.strides
# ey_update_strides = ey_update.strides
# ez_update_strides = ez_update.strides
#
# repeats = [
#     ('full', full_shape, full_strides),
#     ('ncells', ncells_shape, ncells_strides),
#     ('ex', ex_shape, ex_strides),
#     ('ey', ey_shape, ey_strides),
#     ('ez', ez_shape, ez_strides),
#     ('hx', hx_shape, hx_strides),
#     ('hy', hy_shape, hy_strides),
#     ('hz', hz_shape, hz_strides),
#     ('ex_update', ex_update_shape, ex_update_strides),
#     ('ey_update', ey_update_shape, ey_update_strides),
#     ('ez_update', ez_update_shape, ez_update_strides),
#     ('eydz', eydz_shape, eydz_strides),
#     ('ezdy', ezdy_shape, ezdy_strides),
#     ('exdz', exdz_shape, exdz_strides),
#     ('ezdx', ezdx_shape, ezdx_strides),
#     ('exdy', exdy_shape, exdy_strides),
#     ('eydx', eydx_shape, eydx_strides),
#     ('hydz', hydz_shape, hydz_strides),
#     ('hzdy', hzdy_shape, hzdy_strides),
#     ('hxdy', hxdy_shape, hxdy_strides),
#     ('hydx', hydx_shape, hydx_strides),
#     ('hxdz', hxdz_shape, hxdz_strides),
#     ('hzdx', hzdx_shape, hzdx_strides),
# ]
#
# shapes = [(n, d) for n, d, _ in repeats]
# strides = [(n, s) for n, _, s in repeats]
#
# shps = [d for _, d, _ in repeats]
# stds = [s for _, _, s in repeats]
#
# print(f'Number of Shapes: {len(set(shps))} / {len(shps)}')
# print(f'Number of Strides: {len(set(stds))} / {len(stds)}')
#
# def create_extent(name, shape):
#     return f'using {name}_shape_pml = std::extents<std::size_t, {map_shape(shape)}>;'
#
# def create_stride(name, stride):
#     return f'constexpr auto {name}_stride_pml = std::array{{{map_strides(stride)}}};'
#
# def create_str(name, shape, stride):
#     return f'constexpr auto {name} = mdspan_t{{std::extents{{{map_shape(shape)}}}, std::array{{{map_strides(stride)}}}}};'
#
#
# print('\n'.join(create_extent(n, d) for n, d in shapes))
# print('\n'.join(create_stride(n, s) for n, s in strides))
# print('\n'.join(create_str(n, d, s) for n, d, s in repeats))

# def test(arr, ans, name):
#     start = time()
#     for i in range(1000):
#         result = arr[:, :, 1:] - arr[:, :, :-1]
#         assert np.all(result == ans)
#     print(f'constexpr auto {name} Total time: {time() - start}')

# # Test dx
# test(jki, 40000, 'jki')
# test(kji, 40000, 'kji')

# # Test dy
# test(ikj, 200, 'ikj')
# test(kij, 200, 'kij')

# # Test dz
# test(ijk, 1, 'ijk')
# test(jik, 1, 'jik')


# vertices = np.array([[0, 0, 0],
#                      [1, 0, 0],
#                      [1, 1, 0],
#                      [0, 1, 0],
#                      [0, 0, 1],
#                      [1, 0, 1],
#                      [1, 1, 1],
#                      [0, 1, 1]])
#
# edges = [[0, 1], [1, 2], [2, 3], [3, 0],  # Bottom face
#          [4, 5], [5, 6], [6, 7], [7, 4],  # Top face
#          [0, 4], [1, 5], [2, 6], [3, 7]]  # Connecting edges
#
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(111, projection='3d')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
#
# c1 = np.copy(vertices)
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
#
# for e in edges:
#     ax.plot(*zip(*c1[e]), color='k', alpha=0.25)
#     ax.plot(*zip(*c2[e]), color='k', alpha=0.25)
#     ax.plot(*zip(*c3[e]), color='k', alpha=0.25)
#     ax.plot(*zip(*c4[e]), color='k', alpha=0.25)
#     ax.plot(*zip(*c5[e]), color='k', alpha=0.25)
#     ax.plot(*zip(*c6[e]), color='k', alpha=0.25)
#     ax.plot(*zip(*c7[e]), color='k', alpha=0.25)
#     ax.plot(*zip(*c8[e]), color='k', alpha=0.25)
#
#
# nx, nz = 3, 3
#
# dx = dy = dz = np.linspace(0, nx - 1, nx)
# hdx = hdy = hdz = np.linspace(0.5, nx - 1.5, nx - 1)
#
# px, py, pz = np.meshgrid(dx, dy, dz)
#
# Ex_x, Ex_y, Ex_z = np.meshgrid(hdx, dy, dz)
# Ey_x, Ey_y, Ey_z = np.meshgrid(dx, hdy, dz)
# Ez_x, Ez_y, Ez_z = np.meshgrid(dx, dy, hdz)
# Bx_x, Bx_y, Bx_z = np.meshgrid(dx, hdy, hdz)
# By_x, By_y, By_z = np.meshgrid(hdx, dy, hdz)
# Bz_x, Bz_y, Bz_z = np.meshgrid(hdx, hdy, dz)
#
# ax.scatter(Ex_x, Ex_y, Ex_z, c='g', marker='$E_x$', s=180)
# ax.scatter(Ey_x, Ey_y, Ey_z, c='b', marker='$E_y$', s=180)
# ax.scatter(Ez_x, Ez_y, Ez_z, c='c', marker='$E_z$', s=180)
#
# # ax.scatter(Bx_x, Bx_y, Bx_z, c='r', marker='$B_x$', s=180)
# # ax.scatter(By_x, By_y, By_z, c='y', marker='$B_y$', s=180)
# # ax.scatter(Bz_x, Bz_y, Bz_z, c='m', marker='$B_z$', s=180)
#
# # wsx, wsy, wsz = np.meshgrid(dx, dy, dz)
# # ax.scatter(wsx, wsy, wsz, c='tab:purple', marker='h', s=180)
#
# particle = np.array([0.50223, 0.4999, 0.50223])
# ax.scatter(particle[0], particle[1], particle[2], c='tab:orange', marker='X', s=180)
#
# xs = [[0.5, 0.5], [1.5, 1.5]]
# ys = [[0, 2]]
# zs = [[0, 0], [1, 1], [2, 2]]
#
# for lx in xs:
#     for ly in ys:
#         for lz in xs:
#             ax.plot(lx, ly, lz, c='k', ls='--', alpha=0.25)
#             ax.plot(ly, lx, lz, c='k', ls='--', alpha=0.25)
#             ax.plot(lx, lz, ly, c='k', ls='--', alpha=0.25)
#
#
# plt.show()





