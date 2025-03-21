#!/usr/bin/env python3

from adios2 import FileReader
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable


data_dir = '/home/cepheid/TriForce/game_engine/data'

EPS0 = 8.8541878188E-12 # F/m
MU0 = 1.25663706127E-6 # N/A^2

dt = 1.829541541469147e-11

def load_field(n, name, file_dir):
    filename = f'/fields_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
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

def total_field_energy(n, file_dir):
    Ex = load_field(n, 'Ex', file_dir)[:, :-1, :-1]
    Ey = load_field(n, 'Ey', file_dir)[:-1, :, :-1]
    Ez = load_field(n, 'Ez', file_dir)[:-1, :-1, :]
    Hx = load_field(n, 'Hx', file_dir)[:-1, :, :]
    Hy = load_field(n, 'Hy', file_dir)[:, :-1, :]
    Hz = load_field(n, 'Hz', file_dir)[:, :, :-1]

    E = 0.5 * EPS0 * (Ex**2 + Ey**2 + Ez**2)
    H = (0.5 * MU0) * (Hx**2 + Hy**2 + Hz**2)

    total = E + H
    return E.sum() * (0.01**3), H.sum() * (0.01**3), total.sum() * (0.01**3)

def plot_total_energy(start, stop, step):
    def plot_regular(arr1, arr2, log=False):
        times = np.linspace(0, stop * dt, len(arr1))
        fig, ax = plt.subplots(figsize=(8,8))

        # ax.plot(times, arr1[:, 0], label='E-Field float')
        # ax.plot(times, arr1[:, 1], label='B-Field float')
        # ax.plot(times, arr2[:, 0], label='E-Field double')
        # ax.plot(times, arr2[:, 1], label='B-Field double')
        ax.plot(times, arr1[:, 2], label='Total single')
        ax.plot(times, arr2[:, 2], label='Total double')

        if log:
            ax.set_yscale('log')

        ax.grid(True, which='major', axis='both', linestyle='--')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Energy (J)')
        ax.set_title(f'Total Field Energy')
        ax.legend()
        plt.savefig(data_dir + f'/pngs/total_energy.png')
        plt.close(fig)

    targs1 = [(n, '/single_test') for n in range(start, stop + step, step)]
    targs2 = [(n, '/double_test') for n in range(start, stop + step, step)]

    with mp.Pool(16) as p:
        result1 = p.starmap(total_field_energy, targs1)
        result2 = p.starmap(total_field_energy, targs2)

    plot_regular(np.asarray(result1), np.asarray(result2), True)


def compare_fields(n, name, step):
    print(f'Processing file {n}')
    # vmin, vmax = -0.1, 0.1

    single = load_field(n, name, '/single_test')
    double = load_field(n, name, '/double_test')

    diff = double - single

    fig, ax = plt.subplots(figsize=(10, 10))

    nx = diff.shape[0] // 2
    ny = diff.shape[1] // 2
    nz = diff.shape[2] // 2

    im = ax.contourf(diff[:, :, nz].T, levels=100)#, vmin=vmin, vmax=vmax)
    #im = ax.pcolormesh(frame[:, :, nz].T, vmin=vmin, vmax=vmax)

    # fig.colorbar(ScalarMappable(norm=im.norm), ax=ax)
    fig.colorbar(im)

    ax.set_xlabel('x (arb)')
    ax.set_ylabel('y (arb)')
    ax.set_title(f'{name}[:, ny/2, :]')

    plt.savefig(data_dir + f'/pngs/diff_{name}_{n // step:08}.png')
    plt.clf()
    plt.close(fig)



def main():
    start = 0
    stop = 4000
    step = 40

    # targs = [(n, 'electrons', 'density', step) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(plot_metric, targs)

    # field = 'Ez'
    # targs = [(n, field, step) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(plot3d, targs)

    # field = 'Ez'
    # targs = [(n, field, step) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(compare_fields, targs)

    plot_total_energy(start, stop, step)


if __name__ == '__main__':
    main()