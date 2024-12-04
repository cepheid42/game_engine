#!/usr/bin/env python3

import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt

data_path = '/home/cepheid/TriForce/game_engine/data'

EPS0 = 8.8541878188E-12 # F/m
MU0 = 1.25663706127E-6 # N/A^2

def plot1d(n):
    print(f'Plotting file {n:06d}')
    file = data_path + f'/Ez_{n:06d}.csv'
    nx = ny = nz = 120

    fig, ax = plt.subplots()

    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',')

    # ax.plot(data, label='data')
    # ax.set_ylim([-1, 1])

    # nx, ny = data.shape
    # # ax.plot(data[:, ny // 2], label='data')
    # ax.plot(data[nx // 2, :], label='data')
    # ax.set_ylim([-0.01, 0.01])

    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx - 1, ny, nz)) # Ex
    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny - 1, nz)) # Ey
    data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny, nz - 1)) # Ez

    # ax.plot(data[:, ny // 2, nz // 2], label='data')
    # ax.plot(data[nx // 2, :, nz // 2], label='data')
    ax.plot(data[nx // 2, ny // 2, :], label='data')
    ax.set_ylim([-0.01, 0.01])

    plt.savefig(data_path + f'/pngs/ez_{n:06d}.png')
    plt.clf()
    plt.close(fig)


def plot2d(n):
    print(f'Plotting file {n:06d}')
    file = data_path + f'/Ez_{n:06d}.csv'

    data = np.genfromtxt(file, dtype=np.float64, delimiter=',')
    fig, ax = plt.subplots()
    ax.contourf(data, levels=100)

    plt.savefig(data_path + f'/pngs/ez_{n:06d}.png')
    plt.clf()
    plt.close(fig)


def plot3d(n):
    print(f'Plotting file {n:06d}')
    file = data_path + f'/Ez_{n:06d}.csv'
    nx = ny = nz = 120
    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx - 1, ny, nz)) # Ex
    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny - 1, nz)) # Ey
    data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny, nz - 1)) # Ez
    fig, ax = plt.subplots()

    # ax.contourf(data[nx // 2, :, :], levels=100)
    # ax.contourf(data[:, ny // 2, :], levels=100)
    ax.contourf(data[:, :, nz // 2], levels=100)

    plt.savefig(data_path + f'/pngs/ez_{n:06d}.png')
    plt.clf()
    plt.close(fig)


def calculate_total_energy(n):
    print(f'Plotting file {n:06d}')

    nx = ny = nz = 120
    exdata = np.genfromtxt(data_path + f'/Ex_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny, nz))
    eydata = np.genfromtxt(data_path + f'/Ey_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny - 1, nz))
    ezdata = np.genfromtxt(data_path + f'/Ez_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny, nz - 1))
    hxdata = np.genfromtxt(data_path + f'/Hx_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny - 1, nz - 1))
    hydata = np.genfromtxt(data_path + f'/Hy_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny, nz - 1))
    hzdata = np.genfromtxt(data_path + f'/Hz_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny - 1, nz))

    E = exdata[:, :-1, :-1]**2 + eydata[:-1, :, :-1]**2 + ezdata[:-1, :-1, :]**2
    H = hxdata[:-1, :, :]**2 + hydata[:, :-1, :]**2 + hzdata[:, :, :-1]**2

    # return (0.5 * MU0 * H).sum()
    return (0.5 * (EPS0 * E + MU0 * H)).sum()


def plot_total_field_energy(start, nsteps, step):
    fig, ax = plt.subplots()

    targs = [n  for n in range(start, nsteps // step)]

    with mp.Pool(32) as p:
        result = p.map(calculate_total_energy, targs)

    result = np.asarray(result)

    ax.plot(result)
    ax.set_yscale('log')
    plt.show()


def main():
    start = 0
    nsteps = 1000
    step = 10

    plot_total_field_energy(start, nsteps, step)

    # targs = [n for n in range(start, nsteps // step)]
    #
    # with mp.Pool(16) as p:
    #     # p.map(plot1d, targs)
    #     # p.map(plot2d, targs)
    #     p.map(plot3d, targs)



if __name__ == '__main__':
    main()
