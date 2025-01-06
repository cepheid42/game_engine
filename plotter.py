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

    data = np.genfromtxt(file, dtype=np.float64, delimiter=',')

    ax.plot(data, label='data')
    ax.set_ylim([-1.1, 1.1])

    # nx, ny = data.shape
    # ax.plot(data[:, ny // 2], label='data')
    # # ax.plot(data[nx // 2, :], label='data')
    # ax.set_ylim([-1.1, 1.1])

    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx - 1, ny, nz)) # Ex
    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny - 1, nz)) # Ey
    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny, nz - 1)) # Ez

    # ax.plot(data[:, ny // 2, nz // 2], label='data')
    # ax.plot(data[nx // 2, :, nz // 2], label='data')
    # ax.plot(data[nx // 2, ny // 2, :], label='data')
    # ax.set_ylim([-0.01, 0.01])

    plt.savefig(data_path + f'/pngs/Ez_{n:06d}.png')
    plt.close(fig)


def plot2d(n):
    print(f'Plotting file {n:06d}')
    file = data_path + f'/Ez_{n:06d}.csv'
    # file1 = '/home/cepheid/TriForce/game_engine/periodic_y_data' + f'/Ez_{n:06d}.csv'
    # file2 = '/home/cepheid/TriForce/game_engine/py_nohy_data' + f'/Ez_{n:06d}.csv'

    # prime = np.genfromtxt(file1, dtype=np.float64, delimiter=',')
    # test = np.genfromtxt(file2, dtype=np.float64, delimiter=',')

    # data = prime - test

    data = np.genfromtxt(file, dtype=np.float64, delimiter=',')
    fig, ax = plt.subplots()
    im = ax.pcolormesh(data, shading='gouraud')#, vmin=-1, vmax=1)
    plt.colorbar(im)

    plt.savefig(data_path + f'/pngs/Ez_{n:06d}.png')
    plt.close(fig)


def plot3d(n):
    print(f'Plotting file {n:06d}')
    file = data_path + f'/Ez_{n:06d}.csv'
    nx = 124
    ny = 124
    nz = 124

    # f1 = data_path + f'/Ez_psi_{n:06d}.csv'
    # f2 = f'/home/cepheid/TriForce/game_engine/data_bak/Ez_psi_{n:06d}.csv'
    # d1 = np.genfromtxt(f1, dtype=np.float64, delimiter=',').reshape((nx, ny, nz - 1))
    # d2 = np.genfromtxt(f2, dtype=np.float64, delimiter=',').reshape((nx, ny, nz - 1))
    # data = np.abs(d1 - d2)
    # fig, ax = plt.subplots()
    # # im = ax.contourf(data[nx // 2, :, :], levels=100)
    # # im = ax.contourf(data[:, ny // 2, :], levels=100)
    # im = ax.contourf(data[:, :, nz // 2], levels=100)
    # fig.colorbar(im)
    #
    # plt.savefig(data_path + f'/pngs/Ez_{n:06d}.png')
    # plt.close(fig)

    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx - 1, ny, nz)) # Ex
    # data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny - 1, nz)) # Ey
    data = np.genfromtxt(file, dtype=np.float64, delimiter=',').reshape((nx, ny, nz - 1)) # Ez
    fig, ax = plt.subplots()

    # ax.contourf(data[nx // 2, :, :], levels=100)
    # ax.contourf(data[:, ny // 2, :], levels=100)
    im = ax.contourf(data[:, :, nz // 2], levels=100)
    fig.colorbar(im)

    plt.savefig(data_path + f'/pngs/Ez_{n:06d}.png')
    plt.close(fig)


def calculate_total_energy(n):
    print(f'Plotting file {n:06d}')

    nx = ny = nz = 124
    exdata = np.genfromtxt(data_path + f'/Ex_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny, nz))
    eydata = np.genfromtxt(data_path + f'/Ey_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny - 1, nz))
    ezdata = np.genfromtxt(data_path + f'/Ez_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny, nz - 1))
    hxdata = np.genfromtxt(data_path + f'/Hx_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny - 1, nz - 1))
    hydata = np.genfromtxt(data_path + f'/Hy_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny, nz - 1))
    hzdata = np.genfromtxt(data_path + f'/Hz_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny - 1, nz))
    E = exdata[:, :-1, :-1]**2 + eydata[:-1, :, :-1]**2 + ezdata[:-1, :-1, :]**2
    H = hxdata[:-1, :, :]**2 + hydata[:, :-1, :]**2 + hzdata[:, :, :-1]**2

    # ezdata = np.genfromtxt(data_path + f'/Ez_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny))
    # hxdata = np.genfromtxt(data_path + f'/Hx_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny - 1))
    # hydata = np.genfromtxt(data_path + f'/Hy_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny))
    # E = ezdata[:-1, :-1]**2
    # H = hxdata[:-1, :]**2 + hydata[:, :-1]**2

    # exdata = np.genfromtxt(data_path + f'/Ex_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny))
    # eydata = np.genfromtxt(data_path + f'/Ey_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny - 1))
    # hzdata = np.genfromtxt(data_path + f'/Hz_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny - 1))
    # E = exdata[:, :-1]**2 + eydata[:-1, :]**2
    # H = hzdata[:, :]**2

    # # exdata = np.genfromtxt(data_path + f'/Ex_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny))
    # # eydata = np.genfromtxt(data_path + f'/Ey_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny - 1))
    # ezdata = np.genfromtxt(data_path + f'/Ez_{n:06d}.csv', dtype=np.float64, delimiter=',')
    # # hxdata = np.genfromtxt(data_path + f'/Hx_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx, ny - 1))
    # hydata = np.genfromtxt(data_path + f'/Hy_{n:06d}.csv', dtype=np.float64, delimiter=',')
    # # hzdata = np.genfromtxt(data_path + f'/Hz_{n:06d}.csv', dtype=np.float64, delimiter=',').reshape((nx - 1, ny - 1))
    #
    # E = ezdata[:-1]**2
    # H = hydata**2

    return (0.5 * (EPS0 * E + MU0 * H)).sum()


def plot_total_field_energy(start, nsteps, step):
    def plot_regular(arr):
        fig, ax = plt.subplots()
        ax.plot(arr)
        plt.savefig('/home/cepheid/TriForce/game_engine/total_energy.png')
        plt.close(fig)

    def plot_log(arr):
        fig, ax = plt.subplots()
        ax.set_yscale('log')
        ax.plot(arr)
        plt.savefig('/home/cepheid/TriForce/game_engine/total_energy_log.png')
        plt.close(fig)

    targs = [n  for n in range(start, nsteps // step)]

    with mp.Pool(32) as p:
        result = p.map(calculate_total_energy, targs)

    result = np.asarray(result)

    np.savetxt('original_te.csv', result, delimiter=',')

    # plot_regular(result)
    # plot_log(result)

    # fig, ax = plt.subplots()
    # ax.plot(result)
    # ax.set_yscale('log')
    # plt.show()


def compare_total_energy():
    path = '/home/cepheid/TriForce/game_engine/'
    d1 = np.genfromtxt(path + '/original_te.csv', dtype=np.float64, delimiter=',')
    d2 = np.genfromtxt(path + '/special_te.csv', dtype=np.float64, delimiter=',')
    d3 = np.genfromtxt(path + '/extra_te.csv', dtype=np.float64, delimiter=',')

    fig, ax = plt.subplots()

    ax.plot(d1, 'b-', label='Original')
    ax.plot(d2, 'r:', label='Works, but shouldnt')
    ax.plot(d3, 'k', label='No worky')

    # ax.vlines([13, 34], ymin=1.0E-12, ymax=1.0E-3, colors='k')

    ax.set_xlabel('step #')
    ax.set_ylabel(r'$\log{(U)} \quad (\text{kg} / \text{m}^2 \text{s}^2)$')

    ax.set_yscale('log')
    ax.legend()
    plt.show()


def main():
    start = 0
    nsteps = 400
    step = 4

    # compare_total_energy()
    # plot_total_field_energy(start, nsteps, step)

    targs = [n for n in range(start, nsteps // step)]

    with mp.Pool(16) as p:
        # p.map(plot1d, targs)
        p.map(plot2d, targs)
        # p.map(plot3d, targs)


if __name__ == '__main__':
    main()
