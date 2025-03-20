#!/usr/bin/env python3

from adios2 import FileReader
# import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable


data_dir = '/home/cepheid/TriForce/game_engine/data'

EPS0 = 8.8541878188E-12 # F/m
MU0 = 1.25663706127E-6 # N/A^2

dt = 5.000040741175102e-12

def load_field(n, name):
    filename = f'/fields_{n:010d}.bp'
    with FileReader(data_dir + filename) as f:
        return f.read(name)

def plot3d(n, name, step):
    print(f'Processing file {n}')
    vmin, vmax = -0.1, 0.1
    # norm = Normalize(vmin=vmin, vmax=vmax)

    frame = load_field(n, name)

    fig, ax = plt.subplots(figsize=(10, 10))

    nx = frame.shape[0] // 2
    ny = frame.shape[1] // 2
    nz = frame.shape[2] // 2

    im = ax.contourf(frame[:, :, nz].T, levels=100, vmin=vmin, vmax=vmax)
    #im = ax.pcolormesh(frame[:, :, nz].T, vmin=vmin, vmax=vmax)

    fig.colorbar(ScalarMappable(norm=im.norm), ax=ax)
    # fig.colorbar(im)

    ax.set_xlabel('x (arb)')
    ax.set_ylabel('y (arb)')
    ax.set_title(f'{name}[:, ny/2, :]')

    plt.savefig(data_dir + f'/pngs/{name}_{n // step:08}.png')
    plt.clf()
    plt.close(fig)

# def total_field_energy(n):
#     Ex = load_frame(n, 'Ex')[:-1, :, :]
#     Ey = load_frame(n, 'Ey')[:, :-1, :]
#     Ez = load_frame(n, 'Ez')[:, :, :-1]
#     Bx = load_frame(n, 'Bx')[:, :-1, :-1]
#     By = load_frame(n, 'By')[:-1, :, :-1]
#     Bz = load_frame(n, 'Bz')[:-1, :-1, :]
#
#     E = 0.5 * EPS0 * (Ex**2 + Ey**2 + Ez**2)
#     B = (0.5 / MU0) * (Bx**2 + By**2 + Bz**2)
#
#     return E.sum() * (5e-3**3), B.sum() * (5e-3**3), (Ex**2).sum() * 0.5 * EPS0 * (5e-3**3), (Ey**2).sum() * 0.5 * EPS0 * (5e-3**3), (Ez**2).sum() * 0.5 * EPS0 * (5e-3**3)

# def plot_total_energy(start, stop, step):
#     def plot_regular(arr):
#         times = np.linspace(0, stop * dt, len(arr))
#         fig, ax = plt.subplots(figsize=(8,8))
#         # ax.plot(times, arr[:, 1], 'b', label='B-Field')
#         # ax.plot(times, arr[:, 0], 'r', label='E-Field')
#         ax.plot(times, arr[:, 2], 'r', label='Ex-Field')
#         ax.plot(times, arr[:, 3], 'g', label='Ey-Field')
#         ax.plot(times, arr[:, 4], 'b', label='Ez-Field')
#         ax.grid(True, which='major', axis='both', linestyle='--')
#         # ax.set_xlim([0, 0.85e-8])
#         ax.set_xlabel('Time (s)')
#         ax.set_ylabel('Energy (J)')
#         ax.set_title(f'Total Field Energy')
#         ax.legend()
#         plt.savefig(data_dir + f'/pngs/total_energy.png')
#         plt.close(fig)
#
#     def plot_log(arr):
#         fig, ax = plt.subplots()
#         ax.set_yscale('log')
#         ax.plot(arr)
#         plt.savefig(data_dir + f'/pngs/total_energy_log.png')
#         plt.close(fig)
#
#     targs = [n for n in range(start, stop + step, step)]
#
#     with mp.Pool(16) as p:
#         result = p.map(total_field_energy, targs)
#
#     plot_regular(np.asarray(result))

# def load_particles(n, name, metric):
#     filename = f'/{name}_{metric}_{n:010}.h5'
#     with h5py.File(data_dir + filename, 'r') as f:
#         return f[metric][:, :, :]
#
#
# def plot_metric(n, name, metric, step):
#     print(f'Processing file {n}')
#
#     frame = load_particles(n, name, metric)
#     fig, ax = plt.subplots(figsize=(10, 10))
#
#     im = ax.pcolormesh(frame[:, 0, :])
#     fig.colorbar(im)
#
#     ax.set_xlabel('x (arb)')
#     ax.set_ylabel('y (arb)')
#     ax.set_title(f'{name}[:, ny/2, :]')
#
#     plt.savefig(data_dir + f'/pngs/{name}_{n // step:010}.png')
#     plt.clf()
#     plt.close(fig)


def main():
    start = 0
    stop = 4000
    step = 40

    # targs = [(n, 'electrons', 'density', step) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(plot_metric, targs)

    field = 'Ez'
    targs = [(n, field, step) for n in range(start, stop + step, step)]
    with mp.Pool(16) as p:
        p.starmap(plot3d, targs)


if __name__ == '__main__':
    main()