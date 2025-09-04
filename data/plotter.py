#!/usr/bin/env python3

from adios2 import FileReader
import numpy as np
import multiprocessing as mp

import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
from matplotlib import colors
from matplotlib.collections import LineCollection

from matplotlib.patches import Circle
from scipy import constants, special


data_dir = '/home/cepheid/TriForce/game_engine/data'

EPS0 = 8.8541878188E-12 # F/m
MU0 = 1.25663706127E-6 # N/A^2

# dt = 4e-17
# dx, dy, dz = 2e-8, 0.01, 2e-8
# cell_volume = dx * dy * dz
# xmin, xmax = -15e-6, 15e-6
# ymin, ymax = 0.0, 0.01
# zmin, zmax = -15e-6, 15e-6
# nx, ny, nz = 1500, 1, 1500

# dt = 9.655312943196189e-14
# dx, dy, dz = 5.78918e-05, 5.78918e-05, 5.78918e-05
# cell_volume = dx * dy * dz
# xmin, xmax = 0.0, 0.0110573338
# ymin, ymax = 0.0, 5.78918e-05
# zmin, zmax = 0.0, 0.0110573338
# nx, ny, nz = 192, 2, 192

dt = 1.7332498813918236e-12
dx, dy, dz = 0.0001, 0.0001, 0.0001
cell_volume = dx * dy * dz
xmin, xmax = 0.0, 0.005
ymin, ymax = 0.0, 0.005
zmin, zmax = 0.0, 0.005
nx, ny, nz = 51, 51, 51

# s_to_ns = 1.0e9
# T_to_G = 1.0e4
# Vm_to_kVcm = 1.0e-5
# Am_to_Acm = 1.0e-4
# J_to_kJ = 1.0e-3

def colored_line(x, y, c, ax, **lc_kwargs):
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)

    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    lc = LineCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment

    return ax.add_collection(lc)

def particle_positions(start, stop, step, group_name, file_dir):
    frames = []
    for n in range(start, stop, step):
        file = f'/{group_name}_dump_{n:010d}.bp'
        # print(data_dir + file_dir + file)
        with FileReader(data_dir + file_dir + file) as f:
            frames.append(f.read('Position'))

    file = f'/fields_{0:010d}.bp'
    with FileReader(data_dir + file_dir + file) as f:
        Ex = f.read('Ex')#[:, :, :-1]
        Ez = f.read('Ez')#[:-1, :, :]
        # field = f.read('Hy')

    Ex = (Ex[:, :, 1:] + Ex[:, :, :-1]) / dx
    Ez = (Ez[1:, :, :] + Ez[:-1, :, :]) / dz
    field = np.sqrt(Ex**2 + Ez**2)[30:-30, :, 30:-30]
    positions = np.asarray(frames).squeeze()

    fig, ax = plt.subplots(figsize=(12, 12))
    ax.set_aspect('equal')
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title('Electron in Radial E, Uniform B')

    # xs = np.linspace(xmin + dx/2, xmax - dx/2, nx - 1, endpoint=True)
    # zs = np.linspace(zmin + dz/2, zmax - dz/2, nz - 1, endpoint=True)

    xs = np.linspace(xmin, xmax, nx - 1, endpoint=True)[30:-30]
    zs = np.linspace(zmin, zmax, nz - 1, endpoint=True)[30:-30]


    im = ax.contourf(zs, xs, field[:, 0, :].T, levels=200, cmap='viridis')
    # im = ax.pcolormesh(zs, xs, field[:, 0, :], cmap='coolwarm')
    fig.colorbar(im, ax=ax, shrink=0.8)

    nt = positions.shape[0]
    color = np.linspace(0, 1, nt)
    line = colored_line(positions[:, 0], positions[:, 2], color, ax, linewidth=2, cmap='autumn')
    ax.scatter(positions[:, 0], positions[:, 2], c=color, s=40, cmap='autumn')

    cs = ax.contour(zs, xs, field[:, 0, :], levels=3)
    ax.clabel(cs)

    # ax.set_xlim([xmin, xmax])
    # ax.set_ylim([zmin, zmax])
    ax.grid()
    plt.show()
    # plt.savefig(data_dir + f'/positions.png')
    # plt.close(fig)


def plot_distributions(start, stop, step, group_name, file_dir):
    # n = 0
    # file1 = f'/{group_name}_dump_{n:010d}.bp'
    # with FileReader(data_dir + file_dir + file1) as f:
    #     # weight1 = f.read('Weight')
    #     gamma1 = f.read('Gamma')

    n = 0
    n2 = stop
    file1 = f'/{group_name}_dump_{n:010d}.bp'
    file2 = f'/{group_name}_dump_{n2:010d}.bp'
    with FileReader(data_dir + '/lsi_test' + file1) as f:
        weight1 = f.read('Weight')
        gamma1 = f.read('Gamma')

    with FileReader(data_dir + '/lsi_test' + file2) as f:
        weight2 = f.read('Weight')
        gamma2 = f.read('Gamma')

    energy1 = (gamma1 - 1) * constants.m_e * constants.c**2 / constants.elementary_charge
    energy2 = (gamma2 - 1) * constants.m_e * constants.c**2 / constants.elementary_charge
    fig, ax = plt.subplots(figsize=(12,10))

    n, bins, p = ax.hist(energy2, weights=weight2, color='r', bins=1000, density=True, histtype='step', log=False)
    ax.hist(energy1, weights=weight1, color='b', bins=bins, density=True, histtype='step', log=False)

    # plt.grid(ls="--")
    # plt.ylabel("fraction")
    # plt.xlabel(f"energy (eV)")

    filename= data_dir + f'/{group_name}_energy_distribution.png'
    plt.savefig(filename)
    plt.close()


def load_field(n, name, file_dir):
    filename = f'/fields_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        return f.read(name)

def load_particle_data(n, name, group_name, file_dir):
    filename = f'/{group_name}_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        return f.read(name)


def plot_metric(n, step, metric, group_name, file_dir):
    print(f'Processing file {n}')

    filename = f'/{group_name}_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        time = f.read_attribute('Time')

    data = load_particle_data(n, metric, group_name, file_dir)
    nnx, nny, nnz = data.shape

    xs = np.linspace(xmin, xmax, nx - 1)
    zs = np.linspace(zmin, zmax, nz - 1)

    fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')
    # ax.set_aspect('equal')

    # norm = colors.LogNorm(vmin=1e26, vmax=1e28)
    # im = ax.contourf(zs, xs, data[:, 0, :], levels=np.logspace(26, 28, 50), norm=norm, cmap='jet')
    # fig.colorbar(ScalarMappable(norm=norm, cmap='jet'), ax=ax, shrink=0.82)

    # im = ax.contourf(data[:, :, nnz // 2], levels=100)
    # fig.colorbar(im, ax=ax)

    ax.plot(data[:, nny // 2, nnz // 2])

    ax.set_title(f'{group_name.capitalize()} {metric} @ {time:.4e} ns')
    ax.set_ylabel(r'x ($\mu$m)')
    ax.set_xlabel(r'z ($\mu$m)')

    # ax.set_xlim([xmin, xmax])
    # ax.set_ylim([zmin, zmax])

    plt.savefig(data_dir + f'/pngs/{group_name}_{metric}_{n // step:010}.png')
    plt.clf()
    plt.close(fig)


def plot_single_field(n, step, name, file_dir):
    # print(f'Processing file {n}...')

    filename = f'/fields_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        time = f.read_attribute('Time')

    field = load_field(n, name, file_dir)

    nnx, nny, nnz = field.shape

    # xs = np.linspace(xmin, xmax, field.shape[0])
    # xs[xs == 0] = np.nextafter(0, 1)
    # zs = np.linspace(zmin, zmax, field.shape[2])

    fig, ax = plt.subplots(figsize=(10, 10), layout='constrained')
    # fig.supxlabel(r'z ($\mu$m)')
    # fig.supylabel(r'x ($\mu$m)')
    fig.suptitle(f'{name} @ {time:.4e} ns')

    # ax.plot(xs, field[:, 0, 50], label=f'{name}')
    # ax.hlines([-1.0, 1.0], xmin=xs[10], xmax=xs[-10], linestyles='--')

    # vmin, vmax = -6e12, 6e12
    # norm = colors.Normalize(vmin=vmin, vmax=vmax)
    im = ax.contourf(field[:, nny // 2, :], cmap='plasma')
    fig.colorbar(im, ax=ax, format='{x:3.1e}', pad=0.01, shrink=0.8)
    # fig.colorbar(ScalarMappable(norm=norm, cmap='jet'), ax=ax, format='{x:3.1e}', pad=0.01, shrink=1.0)
    ax.set_aspect('equal')

    # plt.show()
    plt.savefig(data_dir + f'/pngs/{name}_{n // step:010}.png')
    plt.clf()
    plt.close(fig)


def plot_fields(n, step, file_dir):
    print(f'Processing file {n}...')

    def plot(name, ax, figure, vmin=None, vmax=None):
        field = load_field(n, name, file_dir)
        nnx, nny, nnz = field.shape
        field = field[:, :, nnz // 2]

        # xs = np.linspace(xmin, xmax, field.shape[0])
        # zs = np.linspace(zmin, zmax, field.shape[1])
        # norm = colors.Normalize(vmin=vmin, vmax=vmax)
        im = ax.pcolormesh(field, cmap='plasma')
        # im = ax.contourf(field, levels=100, cmap='plasma', vmin=vmin, vmax=vmax)
        figure.colorbar(im, ax=ax, format='{x:3.1e}', pad=0.01, shrink=1.0)
        # figure.colorbar(ScalarMappable(norm=norm, cmap='plasma'), ax=ax, format='{x:3.1e}', pad=0.01, shrink=1.0)
        ax.set_aspect('equal')
        ax.set_title(f'{name}')

    filename = f'/fields_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        time = f.read_attribute('Time')

    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(14, 10), layout='constrained', sharex=True, sharey=False)
    fig.supxlabel(r'z ($\mu$m)')
    fig.supylabel(r'x ($\mu$m)')
    fig.suptitle(f'Fields @ {time:.4e} ns')
    plot('Ex', axes[0, 0], fig)
    plot('Ey', axes[0, 1], fig)
    plot('Ez', axes[0, 2], fig)
    plot('Hx', axes[1, 0], fig)
    plot('Hy', axes[1, 1], fig)
    plot('Hz', axes[1, 2], fig)
    plot('Jx', axes[2, 0], fig)
    plot('Jy', axes[2, 1], fig)
    plot('Jz', axes[2, 2], fig)

    plt.savefig(data_dir + f'/pngs/fields_{n // step:010}.png')
    plt.clf()
    plt.close(fig)

def total_field_energy(n, file_dir):
    print(f'Calculating field energy for step {n}')
    Ex = load_field(n, 'Ex', file_dir)[:, :-1, :-1]
    Ey = load_field(n, 'Ey', file_dir)[:-1, :, :-1]
    Ez = load_field(n, 'Ez', file_dir)[:-1, :-1, :]
    Hx = load_field(n, 'Hx', file_dir)[:-1, :, :]
    Hy = load_field(n, 'Hy', file_dir)[:, :-1, :]
    Hz = load_field(n, 'Hz', file_dir)[:, :, :-1]

    E = 0.5 * EPS0 * (Ex**2 + Ey**2 + Ez**2)
    H = 0.5 * MU0 * (Hx**2 + Hy**2 + Hz**2)
    total = E + H
    return total.sum() * cell_volume
    # return E.sum() * cell_volume, H.sum() * cell_volume, total.sum() * cell_volume

def plot_field_energy(start, stop, step, file_dir):
    print('Plotting field energy...')
    print('Processing fields...', end=' ')
    targs1 = [(n, file_dir) for n in range(start, stop, step)]
    with mp.Pool(16) as p:
        result = p.starmap(total_field_energy, targs1)
    result = np.asarray(result)
    print('Done.')

    time = np.linspace(0, stop * dt * 10**15, stop // step)

    fig, ax = plt.subplots(figsize=(8, 8))

    # ax.plot(time, result[:, 0], label='E-Field')
    # ax.plot(time, result[:, 1], label='H-Field')
    # ax.plot(time, result[:, 2], '-.', label='Total')
    ax.plot(result)

    # ax.set_yticks([0.0, 0.0005, 0.001, 0.0015, 0.002])
    # ax.grid(True, which='major', axis='both', linestyle='--')
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('Energy (J)')
    ax.set_title(f'Total Field Energy')
    # ax.legend()
    plt.savefig(data_dir + f'/total_field_energy.png')
    plt.close(fig)
    print('Finished plotting field energy.')


def calculate_KE(n, group_name, file_dir):
    weight = load_particle_data(n, 'Weight', group_name + '_dump', file_dir)
    gamma = load_particle_data(n, 'Gamma', group_name + '_dump', file_dir)
    mass = constants.m_e if group_name == 'electrons' else constants.m_e
    return (weight * (gamma - 1.0) * constants.c**2 * mass).sum()

def plot_KE(start, stop, step, file_dir):
    print('Plotting particle KE...')
    print(f'Processing electrons...', end=' ')
    targs = [(n, 'electrons', file_dir) for n in range(start, stop, step)]
    with mp.Pool(16) as p:
        electron_energy = p.starmap(calculate_KE, targs)
    print('Done.')

    print('Processing ions...', end=' ')
    targs = [(n, 'ions', file_dir) for n in range(start, stop, step)]
    with mp.Pool(16) as p:
        ion_energy = p.starmap(calculate_KE, targs)
    print('Done')

    electron_energy = np.asarray(electron_energy)
    ion_energy = np.asarray(ion_energy)

    # time = np.linspace(0, stop * dt, stop // step)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(electron_energy, label='Electrons')
    ax.plot(ion_energy, label='Ions')

    # ax.set_xlabel('Time (fs)')
    ax.set_ylabel('KE (J)')
    ax.set_title(f'Total Particle KE')
    ax.legend()
    # plt.show()
    plt.savefig(data_dir + f'/total_particle_energy.png')
    plt.clf()
    plt.close(fig)
    print('Finished plotting particle KE.')


def divE(start, stop, step, file_dir):
    error = []
    for n in range(start, stop + step, step):
        Ex = load_field(n, 'Ex', file_dir)[:, :-1, :-1]
        # Ey = load_field(n, 'Ey', file_dir)[:-1, :, :-1]
        Ez = load_field(n, 'Ez', file_dir)[:-1, :-1, :]

        rho = load_particle_data(n, 'Density', 'electrons', file_dir).squeeze()[:-1, :-1]
        rho0 = constants.e / (dx * dy * dz)

        divergence = np.diff(Ex, axis=0).squeeze()[:, :-1] + np.diff(Ez, axis=2).squeeze()[:-1, :]

        R = (divergence - rho + rho0) / rho0
        error.append(np.sqrt((R**2).mean()))

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(error)
    plt.show()


def main():
    # step = 150
    # start = 0
    # stop = 30000

    step = 160
    start = 0
    stop = 16000

    file_dir = '/positron_test'

    # divE(start, stop, step, file_dir)

    # plot_distributions(start, stop, step, 'electrons', file_dir)
    # plot_distributions(start, stop, step, 'ions', file_dir)

    # targs = [(n, step, 'Density', 'electrons', file_dir) for n in range(start, stop, step)]
    # with mp.Pool(16) as p:
    #    p.starmap(plot_metric, targs)

    # targs = [(n, step, 'Temperature', 'electrons', file_dir) for n in range(start, stop, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(plot_metric, targs)

    # targs = [(n, step, 'Density', 'ions', file_dir) for n in range(start, stop, step)]
    # with mp.Pool(16) as p:
    #    p.starmap(plot_metric, targs)

    # targs = [(n, step, 'Temperature', 'ions', file_dir) for n in range(start, stop, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(plot_metric, targs)

    # targs = [(n, step, file_dir) for n in range(start, stop, step)]
    # with mp.Pool(16) as p:
    #    p.starmap(plot_fields, targs)

    # targs = [(n, step, 'Jz', file_dir) for n in range(start, stop, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(plot_single_field, targs)

    # particle_positions(start, stop, step, 'electrons', file_dir)

    plot_KE(start, stop, step, file_dir)
    plot_field_energy(start, stop, step, file_dir)

    # plot_single_field(0, 1, 'Ex', file_dir)
    # plot_single_field(0, 1, 'Ez', file_dir)


if __name__ == '__main__':
    main()
