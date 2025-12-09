#!/usr/bin/env python3

from adios2 import FileReader
import numpy as np
import multiprocessing as mp

import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
from matplotlib import colors, ticker
from matplotlib.collections import LineCollection

from scipy import constants, special

data_dir = '/home/cepheid/TriForce/game_engine/data'

EPS0 = 8.8541878188E-12 # F/m
MU0 = 1.25663706127E-6 # N/A^2

s_to_ns = 1.0e9
s_to_ps = 1.0e12
T_to_G = 1.0e4
Vm_to_kVcm = 1.0e-5
Am_to_Acm = 1.0e-4
J_to_kJ = 1.0e-3
H_to_B = MU0

# def colored_line(x, y, c, ax, **lc_kwargs):
#     default_kwargs = {"capstyle": "butt"}
#     default_kwargs.update(lc_kwargs)
#
#     x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
#     y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))
#
#     coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
#     coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
#     coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
#     segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)
#
#     lc = LineCollection(segments, **default_kwargs)
#     lc.set_array(c)  # set the colors of each segment
#
#     return ax.add_collection(lc)

# def particle_positions(start, stop, step, group_name, file_dir):
#     frames = []
#     for n in range(start, stop, step):
#         file = f'/{group_name}_dump_{n:010d}.bp'
#         with FileReader(data_dir + file_dir + file) as f:
#             frames.append(f.read('Position'))
#
#     file = f'/fields_{0:010d}.bp'
#     with FileReader(data_dir + file_dir + file) as f:
#         Ex = f.read('Ex')#[:, :, :-1]
#         Ez = f.read('Ez')#[:-1, :, :]
#         # field = f.read('Hy')
#
#     Ex = (Ex[:, :, 1:] + Ex[:, :, :-1]) / dx
#     Ez = (Ez[1:, :, :] + Ez[:-1, :, :]) / dz
#     field = np.sqrt(Ex**2 + Ez**2)[30:-30, :, 30:-30]
#     positions = np.asarray(frames).squeeze()
#
#     fig, ax = plt.subplots(figsize=(12, 12))
#     ax.set_aspect('equal')
#     ax.set_xlabel("x (m)")
#     ax.set_ylabel("z (m)")
#     ax.set_title('Electron in Radial E, Uniform B')
#
#     # xs = np.linspace(xmin + dx/2, xmax - dx/2, nx - 1, endpoint=True)
#     # zs = np.linspace(zmin + dz/2, zmax - dz/2, nz - 1, endpoint=True)
#
#     xs = np.linspace(xmin, xmax, nx - 1, endpoint=True)[30:-30]
#     zs = np.linspace(zmin, zmax, nz - 1, endpoint=True)[30:-30]
#
#
#     im = ax.contourf(zs, xs, field[:, 0, :].T, levels=200, cmap='viridis')
#     # im = ax.pcolormesh(zs, xs, field[:, 0, :], cmap='coolwarm')
#     fig.colorbar(im, ax=ax, shrink=0.8)
#
#     nt = positions.shape[0]
#     color = np.linspace(0, 1, nt)
#     line = colored_line(positions[:, 0], positions[:, 2], color, ax, linewidth=2, cmap='autumn')
#     ax.scatter(positions[:, 0], positions[:, 2], c=color, s=40, cmap='autumn')
#
#     cs = ax.contour(zs, xs, field[:, 0, :], levels=3)
#     ax.clabel(cs)
#
#     # ax.set_xlim([xmin, xmax])
#     # ax.set_ylim([zmin, zmax])
#     ax.grid()
#     plt.show()
#     # plt.savefig(data_dir + f'/positions.png')
#     # plt.close(fig)


def plot_distributions(start, stop, step, group_name, file_dir):
    n = 0
    n2 = stop
    file1 = f'/{group_name}_dump_{n:010d}.bp'
    file2 = f'/{group_name}_dump_{n2:010d}.bp'
    with FileReader(data_dir + file_dir + file1) as f:
        weight1 = f.read('Weight')
        gamma1 = f.read('Gamma')

    with FileReader(data_dir + file_dir + file2) as f:
        weight2 = f.read('Weight')
        gamma2 = f.read('Gamma')

    energy1 = (gamma1 - 1) * constants.m_e * constants.c**2 / constants.elementary_charge
    energy2 = (gamma2 - 1) * constants.m_e * constants.c**2 / constants.elementary_charge
    fig, ax = plt.subplots(figsize=(12,10))

    n, bins, p = ax.hist(energy2, color='r', bins=1000, histtype='step', log=True)
    ax.hist(energy1, color='b', bins=bins, histtype='step', log=True)

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


def plot_metric(n, step, metric, group_name, file_dir):
    filename = f'/{group_name}_dump_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        data = f.read('Weight')
        # sim_step = f.read_attribute('step')
        # dt = f.read_attribute('dt')
        # dims = f.read_attribute('dims')
        # x_range = f.read_attribute('x_range')
        # y_range = f.read_attribute('y_range')
        # z_range = f.read_attribute('z_range')


    # time = sim_step * dt * s_to_ns
    # nx, ny, nz = np.asarray(dims, dtype=int)
    # xmin, xmax = x_range
    # ymin, ymax = y_range
    # zmin, zmax = z_range

    # xs = np.linspace(xmin, xmax, nx)
    # ys = np.linspace(ymin, ymax, ny)
    # zs = np.linspace(zmin, zmax, nz)

    fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')
    # ax.set_aspect('equal')

    norm = colors.LogNorm(vmin=1e24, vmax=1e28)
    im = ax.contourf(data[:, 0, :], levels=np.logspace(24, 28, 50), norm=norm, cmap='jet')
    fig.colorbar(ScalarMappable(norm=norm, cmap='jet'), ax=ax, shrink=0.82)

    # ax.set_title(f'{group_name.capitalize()} {metric} @ {time:.4e} ns')
    ax.set_ylabel(r'x ($\mu$m)')
    ax.set_xlabel(r'z ($\mu$m)')

    # ax.set_xlim([xmin, xmax])
    # ax.set_ylim([zmin, zmax])

    plt.savefig(data_dir + f'/pngs/{group_name}_{metric}_{n // step:010}.png')
    plt.clf()
    plt.close(fig)


def plot_single_field(n, step, name, file_dir):
    # filename = f'/fields_{n:010d}.bp'
    # with FileReader(data_dir + file_dir + filename) as f:
    #     time = f.read_attribute('Time')

    field = load_field(n, name, file_dir)

    nnx, nny, nnz = field.shape

    fig, ax = plt.subplots(figsize=(10, 10), layout='constrained')
    # fig.supxlabel(r'z ($\mu$m)')
    # fig.supylabel(r'x ($\mu$m)')
    # fig.suptitle(f'{name} @ {time * s_to_ns:.4e} ns')

    im = ax.pcolormesh(field[:, 0, :])
    fig.colorbar(im, ax=ax, format='{x:3.1e}', pad=0.01, shrink=0.8)
    ax.set_aspect('equal')

    # plt.show()
    plt.savefig(data_dir + f'/pngs/{name}_{n // step:010}.png')
    plt.clf()
    plt.close(fig)


def plot_fields(n, step, file_dir):
    def plot(name, ax, figure):
        field = load_field(n, name, file_dir)
        nnx, nny, nnz = field.shape
        field = field[:, :, nnz // 2]
        # field = field[:, 0, :]

        # if name[0] == 'H':
        #     field *= H_to_B

        # im = ax.pcolormesh(field)
        im = ax.contourf(field, levels=100)
        figure.colorbar(im, ax=ax, format='{x:3.1e}', pad=0.01, shrink=0.5)
        ax.set_aspect('equal')
        ax.set_title(f'{name}')

    # filename = f'/fields_{n:010d}.bp'
    # with FileReader(data_dir + file_dir + filename) as f:
    #     step = f.read_attribute('step')

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 8), layout='constrained', sharex=True, sharey=False)
    fig.supxlabel(r'z ($\mu$m)')
    fig.supylabel(r'x ($\mu$m)')
    # fig.suptitle(f'Fields @ {time:.4e} ns')
    plot('Ex', axes[0, 0], fig)
    plot('Ey', axes[0, 1], fig)
    plot('Ez', axes[0, 2], fig)
    plot('Hx', axes[1, 0], fig)
    plot('Hy', axes[1, 1], fig)
    plot('Hz', axes[1, 2], fig)
    # plot('Jx', axes[2, 0], fig)
    # plot('Jy', axes[2, 1], fig)
    # plot('Jz', axes[2, 2], fig)

    plt.savefig(data_dir + f'/pngs/fields_{n // step:010}.png')
    plt.clf()
    plt.close(fig)

def total_field_energy(n, file_dir):
    filename = f'/fields_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        ex = f.read('Ex')[:, :-1, :-1]
        ey = f.read('Ey')[:-1, :, :-1]
        ez = f.read('Ez')[:-1, :-1, :]
        hx = f.read('Hx')[:-1, :, :]
        hy = f.read('Hy')[:, :-1, :]
        hz = f.read('Hz')[:, :, :-1]
        # cell_volume = f.read_attribute('cell_volume')

    cell_volume = 0.01**3
    e_field = 0.5 * EPS0 * (ex**2 + ey**2 + ez**2)
    h_field = 0.5 * MU0 * (hx**2 + hy**2 + hz**2)
    return (e_field + h_field).sum() * cell_volume

def plot_field_energy(start, stop, step, file_dir):
    print('Calculating total field energy...', end=' ')
    filename = f'/fields_{0:010d}.bp'
    # with FileReader(data_dir + file_dir + filename) as f:
    #     dt = f.read_attribute('dt')

    targs = [(n, file_dir) for n in range(start, stop + step, step)]
    with mp.Pool(16) as p:
        result = p.starmap(total_field_energy, targs)
    result = np.asarray(result)

    # time = np.linspace(0, stop * dt * s_to_ns, (stop + step) // step)
    fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')
    ax.plot(result)
    # ax.plot(time, result)
    # ax.set_yticks([0, 0.0005, 0.0010, 0.0015, 0.0020])
    # ax.set_xlim([time[0], time[-1]])
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Energy (J)')
    ax.set_title(f'Total Field Energy')
    # ax.legend()
    plt.savefig(data_dir + f'/total_field_energy.png')
    plt.close(fig)
    print('Done.')

def calculate_Temp(n, group_name, file_dir):
    filename = f'/{group_name}_dump_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        mass = f.read_attribute('Mass')
        weight = f.read('Weight')
        velocity = f.read('Velocity')

    ttl_weight = weight.sum(axis=0)[0]
    avg_velocity = (weight * velocity).sum(axis=0) / ttl_weight
    dv = velocity - avg_velocity[None, :]
    ttl_sum_dv2 = (weight * (dv**2).sum(axis=1)[:, None]).sum(axis=0)[0]
    avg_temp = ttl_sum_dv2 * mass / (3.0 * constants.e * ttl_weight)
    return avg_temp

def calculate_KE(n, group_name, file_dir):
    filename = f'/{group_name}_dump_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        mass = f.read_attribute('Mass')
        weight = f.read('Weight')
        gamma = f.read('Gamma')
    return (weight * (gamma - 1.0) * constants.c**2 * mass).sum()

def plot_KE(groups, start, stop, step, file_dir):
    print('Calculating total particle energies...', end=' ')
    filename = f'/{groups[0]}_dump_{0:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        dt = f.read_attribute('dt')

    group_data = {}
    for g in groups:
        targs = [(n, g, file_dir) for n in range(start, stop + step, step)]
        with mp.Pool(16) as p:
            ke = p.starmap(calculate_KE, targs)
        group_data[g] = np.asarray(ke)

    time = np.linspace(0, stop * dt * s_to_ns, (stop + step) // step)
    fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')

    for name, data in group_data.items():
        ax.plot(time, data, label=name)

    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.set_xlim([time[0], time[-1]])
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('KE (J)')
    ax.set_title(f'Total Particle KE')
    ax.legend()
    plt.savefig(data_dir + f'/total_particle_energy.png')
    plt.clf()
    plt.close(fig)
    print('Done.')


def plot_Temp(groups, start, stop, step, file_dir):
    print('Calculating total particle temperatures...', end=' ')
    filename = f'/{groups[0]}_dump_{stop:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        dt = f.read('dt')
        ts = f.read('Time')

    group_data = {}
    for g in groups:
        targs = [(n, g, file_dir) for n in range(start, stop + step, step)]
        with mp.Pool(16) as p:
            temp = p.starmap(calculate_Temp, targs)
        group_data[g] = np.asarray(temp)

    time = np.linspace(0, ts, (stop + step) // step)
    fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')

    for name, data in group_data.items():
        ax.plot(time * s_to_ps, data, label=name)

    eV_to_erg = 1.60218e-12
    T1 = 1250
    T2 = 250
    m = 1.9945e-23 # g
    n = 1e20 # cm^-3
    Z = 6
    lnL = 10.0

    coef = (8 / 3) * np.sqrt(np.pi)
    Ts = 0.5 * (T1 + T2) * eV_to_erg
    mu = coef * (Z * 4.8e-10)**4 * lnL * n / (m**0.5 * Ts**1.5)

    Tc1 = 0.5 * (T1 + T2) + 0.5 * (T1 - T2) * np.exp(-mu * time)
    Tc2 = 0.5 * (T1 + T2) + 0.5 * (T2 - T1) * np.exp(-mu * time)

    # Tc1 = Tc1 * constants.e / constants.k
    # Tc2 = Tc2 * constants.e / constants.k

    ax.plot(time * s_to_ps, Tc1, label='Carbon1 Theory')
    ax.plot(time * s_to_ps, Tc2, label='Carbon2 Theory')

    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('T (eV)')
    ax.set_title(f'Total Particle Temperature')
    ax.legend()
    plt.savefig(data_dir + f'/total_particle_temp.png')
    plt.clf()
    plt.close(fig)
    print('Done.')


def main():
    step = 4
    start = 0
    stop = 400

    file_dir = '/em_test'

    # plot_Temp(['carbon1', 'carbon2'], start, stop, step, file_dir)

    # plot_distributions(start, stop, step, 'carbon1', file_dir)
    # plot_distributions(start, stop, step, 'carbon2', file_dir)

    # targs = [(n, step, 'Density', 'electrons', file_dir) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #    p.starmap(plot_metric, targs)

    # targs = [(n, step, 'Density', 'ions', file_dir) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #    p.starmap(plot_metric, targs)

    targs = [(n, step, file_dir) for n in range(start, stop + step, step)]
    with mp.Pool(8) as p:
       p.starmap(plot_fields, targs)

    # targs = [(n, step, 'Ez', file_dir) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(plot_single_field, targs)

    # particle_positions(start, stop, step, 'electrons', file_dir)

    # plot_Temp(['electrons', 'ions'], start, stop, step, file_dir)
    # plot_KE(['electrons', 'ions'], start, stop, step, file_dir)
    # plot_field_energy(start, stop, step, file_dir)
    # for i in range(start, stop + step, step):
    #     plot_single_field(i, step, 'Hx', file_dir)
    # plot_single_field(0, 1, 'Ez', file_dir)


if __name__ == '__main__':
    main()
