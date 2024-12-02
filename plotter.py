#!/usr/bin/env python3

import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt



# nt = 400
# step = 4

data_path = '/home/cepheid/TriForce/game_engine/data'

def plot1d(n):
    print(f'Plotting file {n:04d}')
    file = data_path + f'/Ez_{n:04d}.csv'
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

    plt.savefig(data_path + f'/pngs/ez_{n:04d}.png')
    plt.clf()
    plt.close(fig)


def plot2d(n):
    print(f'Plotting file {n:04d}')
    file = data_path + f'/Ez_{n:04d}.csv'

    data = np.genfromtxt(file, dtype=np.float64, delimiter=',')
    fig, ax = plt.subplots()
    ax.contourf(data, levels=100)

    plt.savefig(data_path + f'/pngs/ez_{n:04d}.png')
    plt.clf()
    plt.close(fig)


def plot3d(n):
    print(f'Plotting file {n:04d}')
    file = data_path + f'/Ez_{n:04d}.csv'
    nx = ny = nz = 120
    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx - 1, ny, nz)) # Ex
    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny - 1, nz)) # Ey
    data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny, nz - 1)) # Ez
    fig, ax = plt.subplots()

    ax.contourf(data[nx // 2, :, :], levels=100)
    # ax.contourf(data[:, ny // 2, :], levels=100)
    # ax.contourf(data[:, :, nz // 2], levels=100)

    plt.savefig(data_path + f'/pngs/ez_{n:04d}.png')
    plt.clf()
    plt.close(fig)


def main():
    start = 0
    nsteps = 400
    step = 4

    targs = [n for n in range(start, nsteps // step)]

    with mp.Pool(16) as p:
        # p.map(plot1d, targs)
        # p.map(plot2d, targs)
        p.map(plot3d, targs)



if __name__ == '__main__':
    main()
