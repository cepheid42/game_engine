#!/usr/bin/env python3

import numpy as np

import matplotlib.pyplot as plt



nt = 400
step = 4

data_path = '/home/cepheid/TriForce/game_engine/data'

def plot1d(file, n):
    data = np.genfromtxt(file, dtype=np.float64, delimiter=',')
    fig, ax = plt.subplots()

    ax.plot(data, label='data')
    ax.set_ylim([-1, 1])

    plt.savefig(data_path + f'/pngs/ez_{n:04d}.png')
    plt.clf()
    plt.close(fig)


def plot2d(file, n):
    data = np.genfromtxt(file, dtype=np.float64, delimiter=',')
    fig, ax = plt.subplots()
    ax.contourf(data)

    plt.savefig(data_path + f'/pngs/ez_{n:04d}.png')
    plt.clf()
    plt.close(fig)


def plot3d(file, n):
    nx = ny = nz = 104
    data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny, nz - 1)) # Ez
    fig, ax = plt.subplots()

    # ax.contourf(data[nx // 2, :, :])
    # ax.contourf(data[:, ny // 2, :])
    ax.contourf(data[:, :, nz // 2])

    plt.savefig(data_path + f'/pngs/ez_{n:04d}.png')
    plt.clf()
    plt.close(fig)


def main():
    for n in range(0, nt // step):
        print(f'Step {n}')
        file = data_path + f'/Ez_{n:04d}.csv'

        plot1d(file, n)
        # plot2d(file, n)
        # plot3d(file, n)


if __name__ == '__main__':
    main()
