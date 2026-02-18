#!/usr/bin/env python3

from adios2 import FileReader
import numpy as np
import multiprocessing as mp

import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
from matplotlib import colors, ticker
# from matplotlib.collections import LineCollection

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
    with FileReader(file_dir + filename) as f:
        return f.read(name)


def plot_metric(n, step, metric, group_name, file_dir):

    filename = f'/{group_name}_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        data = f.read('Weight')
        sim_step = f.read_attribute('step')
        dt = f.read_attribute('dt')
        dims = f.read_attribute('dims')
        x_range = f.read_attribute('x_range')
        y_range = f.read_attribute('y_range')
        z_range = f.read_attribute('z_range')


    time = sim_step * dt * s_to_ns
    nx, ny, nz = np.asarray(dims, dtype=int)
    xmin, xmax = x_range
    ymin, ymax = y_range
    zmin, zmax = z_range

    xs = np.linspace(xmin, xmax, nx)
    ys = np.linspace(ymin, ymax, ny)
    zs = np.linspace(zmin, zmax, nz)

    fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')
    # ax.set_aspect('equal')

    norm = colors.LogNorm(vmin=1e24, vmax=1e28)
    im = ax.contourf(zs, xs, data[:, 0, :], levels=np.logspace(24, 28, 50), norm=norm, cmap='jet')
    fig.colorbar(ScalarMappable(norm=norm, cmap='jet'), ax=ax, shrink=0.82)

    ax.set_title(f'{group_name.capitalize()} {metric} @ {time:.4e} ns')
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

    im = ax.pcolormesh(field[:, nny // 2, :], cmap='coolwarm')
    fig.colorbar(im, ax=ax, format='{x:3.1e}', pad=0.01, shrink=0.8)
    ax.set_aspect('equal')

    # plt.show()
    plt.savefig(data_dir + f'/pngs/{name}_{n // step:010}.png')
    plt.clf()
    plt.close(fig)
    # return np.max(field[:, nny // 2, :])


def plot_fields(n, step, file_dir):
    def plot(name, ax, figure):
        field = load_field(n, name, file_dir + '/lsi_test')
        nnx, nny, nnz = field.shape
        # field = field[:, :, nnz // 2]
        field = field[:, 0, :]

        # if name[0] == 'H':
        #     field *= H_to_B

        # xs = np.linspace(xmin, xmax, field.shape[0])
        # zs = np.linspace(zmin, zmax, field.shape[1])
        # norm = colors.Normalize(vmin=-10**15, vmax=10**15)
        im = ax.pcolormesh(field)
        figure.colorbar(im, ax=ax, format='{x:3.1e}', pad=0.01, shrink=0.5)
        # figure.colorbar(ScalarMappable(norm=norm, cmap='coolwarm'), ax=ax, format='{x:3.1e}', pad=0.01, shrink=1.0)
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
    plot('Bx', axes[1, 0], fig)
    plot('By', axes[1, 1], fig)
    plot('Bz', axes[1, 2], fig)
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
        hx = f.read('Bx')[:-1, :, :]
        hy = f.read('By')[:, :-1, :]
        hz = f.read('Bz')[:, :, :-1]
        cell_volume = f.read_attribute('cell_volume')

    e_field = 0.5 * EPS0 * (ex**2 + ey**2 + ez**2)
    h_field = 0.5 * MU0 * (hx**2 + hy**2 + hz**2)
    return (e_field + h_field).sum() * cell_volume

def plot_field_energy(start, stop, step, file_dir):
    print('Calculating total field energy...', end=' ')
    filename = f'/fields_{0:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        dt = f.read_attribute('dt')

    targs = [(n, file_dir) for n in range(start, stop + step, step)]
    with mp.Pool(16) as p:
        result = p.starmap(total_field_energy, targs)
    result = np.asarray(result)

    time = np.linspace(0, stop * dt * s_to_ns, (stop + step) // step)
    fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')

    ax.plot(time, result)
    ax.set_yticks([0, 0.0005, 0.0010, 0.0015, 0.0020])
    ax.set_xlim([time[0], time[-1]])
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
    T1 = 1250 / np.sqrt(3)
    T2 = 250 / np.sqrt(3)
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
    plt.show()



def ionization(start, stop, step, file_dir):
    ionized = []
    ion_num = []
    for n in range(start, stop + step, step):
        ionized_al = f'/Al+_dump_{n:010d}.bp'
        with FileReader(data_dir + file_dir + ionized_al) as f:
            data = f.read('Weight')
            ionized.append(data.sum())
            ion_num.append(data.shape[0])

    ionized = np.asarray(ionized)

    v_beam = 1.32523e7
    sigma = 1.428e-20
    e_den = 1.1e27
    time_thry = np.linspace(0, 6, 100) * 0.53e-15
    charge_thry = 1.0 - np.exp(-v_beam * e_den * sigma * time_thry)

    dt = 5e-18
    time = np.arange(start, stop + step, step) * dt
    mean_ion_charge = ionized / 6.6e10

    fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')
    ax.plot(time_thry, charge_thry)
    ax.plot(time, mean_ion_charge)

    plt.show()
    # plt.savefig(data_dir + f'/total_particle_temp.png')
    # plt.close(fig)

def fusion(stop, groups):
    def thermal_reactivity_BH1992(T, bg, mrc2, c):
        T2 = T**2
        T3 = T**3
        theta = T / (1.0 - ((T * c[1] + T2 * c[3] + T3 * c[5]) / (1.0 + T * c[2] + T2 * c[4] + T3 * c[6])))
        xsi = (bg**2 / (4.0 * theta))**(1.0/3.0)
        return c[0] * theta * np.sqrt(xsi / (mrc2 * T3)) * np.exp(-3.0 * xsi)

    # theory solution for second moment of neutron spectrum (expressed as FWHM) is from
    # Ballabio, Kallne, Gorini 1998, Relativistic calculation of fusion product spectra for thermonuclear plasmas
    def calc_FWHM_Ballabio1998(T, alpha, w0):
        delta_w = alpha[0] * T ** (2 / 3) / (1.0 + alpha[1] * T ** alpha[2]) + (alpha[3] * T)
        whalf = w0 * (1.0 + delta_w) * np.sqrt(T)

        return whalf

    # def find_nearest(array, value):
    #     return (np.abs(array - value)).argmin()

    def calc_sigma_fwhm_dyde(weight, gamma, t_final, delta_ij):
        ttl_weight = weight.sum()
        dYdt = ttl_weight / t_final

        sim_volume = 1.0e-6 ** 3
        density = 1e26
        sigma = 1e6 * dYdt * (1.0 + delta_ij) / (density**2 * sim_volume) # m^3 -> cm^3

        E_ttl = (weight * (gamma - 1.0) * constants.m_n * constants.c ** 2 * 1e-3 / constants.e).sum() / ttl_weight
        energy = ((gamma - 1.0) * constants.m_n * constants.c**2 * 1e-3 / constants.e)
        dE2_ttl = ((E_ttl - energy) ** 2 * weight).sum()
        fwhm = np.sqrt(8.0 * np.log(2) * dE2_ttl / ttl_weight)

        return sigma, fwhm

    T_keV_theory = np.linspace(0.1, 20, 100)

    plot_DD = 'DD' in groups
    plot_DT = 'DT' in groups

    if plot_DD:
        # DD Collisions
        dd_reactivity_theory = thermal_reactivity_BH1992(
            T_keV_theory,
            31.3970,  # Bg, sqrt(keV)
            937814,    # mrc2, keV
            [5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0]
        )

        dd_fwhm_theory = calc_FWHM_Ballabio1998(
            T_keV_theory,
            [1.7013e-3, 0.16888, 0.49, 7.9460e-4],
            82.542 # omega0, sqrt(keV)
        )

    if plot_DT:
        # DT Collisions
        dt_reactivity_theory = thermal_reactivity_BH1992(
            T_keV_theory,
            34.3827,  # Bg, sqrt(keV)
            1124656,    # mrc2, keV
            [1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-6]
        )

        dt_fwhm_theory = calc_FWHM_Ballabio1998(
            T_keV_theory,
            [5.1068e-4, 7.6223e-3, 1.78, 8.7691e-5],
            177.259 # omega0, sqrt(keV)
        )

    # start plot
    plt.style.use('triforce.mplstyle')
    good_colors = ['#db6d00', '#006ddb', '#920000', '#52a736', '#9B30FF']

    fig, ax = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')

    ax[0].set_xlabel(r'$T_{\mathrm{i}}$ (keV)')
    ax[0].set_ylabel(r'$\langle \sigma v \rangle$ (cm$^3$/s)')
    ax[0].set_xlim([0, 20])
    ax[0].set_yscale('log')
    ax[0].set_ylim([1.0e-22, 1e-15])
    ax[0].grid()

    ax[1].set_xlabel(r'$T_{\mathrm{i}}$ (keV)')
    ax[1].set_ylabel(r'FWHM (keV)')
    ax[1].set_xlim([0, 20])
    ax[1].set_yscale('log')
    ax[1].set_ylim([99.9, 1.5e3])  # set min at 99.9 instead of 100 so intermediate ticks not printed
    ax[1].yaxis.set_label_coords(-0.1, 0.5)
    ax[1].grid(which='both')

    if plot_DD:
        ax[0].plot(T_keV_theory, dd_reactivity_theory, '--k', label='DD Theory')
        ax[1].plot(T_keV_theory, dd_fwhm_theory, '--k', label='DD Theory')

    if plot_DT:
        ax[0].plot(T_keV_theory, dt_reactivity_theory, '-k', label='DT Theory')
        ax[1].plot(T_keV_theory, dt_fwhm_theory, '-k', label='DT Theory')

    temps = [5, 10, 19]
    DD_reactivity_sim = np.zeros(3)
    DD_fwhm = np.zeros(3)
    DT_reactivity_sim = np.zeros(3)
    DT_fwhm = np.zeros(3)

    for i, t in enumerate(temps):
        if plot_DD:
            # plot DD
            dd_file = f'/DD_fusion_{t}keV/neutrons_dump_{stop:010d}.bp'
            with FileReader(data_dir + dd_file) as f:
                final_time = f.read("Time")[0]
                weights = f.read('Weight')
                gammas = f.read('Gamma')

            sigmaDD, fwhmDD = calc_sigma_fwhm_dyde(weights, gammas, final_time, 1.0)
            DD_reactivity_sim[i] = sigmaDD
            DD_fwhm[i] = fwhmDD

        if plot_DT:
            # plot DT
            dt_file = f'/DT_fusion_{t}keV/neutrons_dump_{stop:010d}.bp'
            with FileReader(data_dir + dt_file) as f:
                final_time = f.read("Time")[0]
                weights = f.read('Weight')
                gammas = f.read('Gamma')

            sigmaDT, fwhmDT = calc_sigma_fwhm_dyde(weights, gammas, final_time, 0.0)
            DT_reactivity_sim[i] = sigmaDT
            DT_fwhm[i] = fwhmDT

    if plot_DD:
        ax[0].plot(temps, DD_reactivity_sim, marker='o', c=good_colors[3], linestyle='none', ms=10, markeredgewidth=1, markeredgecolor='k', label='D-D')
        ax[1].plot(temps, DD_fwhm, marker='o', c=good_colors[3], linestyle='none', ms=10, markeredgewidth=1, markeredgecolor='k', label='D-D')

    if plot_DT:
        ax[0].plot(temps, DT_reactivity_sim, marker='o', c=good_colors[2], linestyle='none', ms=10, markeredgewidth=1, markeredgecolor='k', label='D-T')
        ax[1].plot(temps, DT_fwhm, marker='o', c=good_colors[2], linestyle='none', ms=10, markeredgewidth=1, markeredgecolor='k', label='D-T')

    ax[0].legend(loc=4)
    ax[0].text(0.04, 0.9, '(a)', fontsize=18, weight='bold', transform=ax[0].transAxes)
    ax[1].text(0.04, 0.9, '(b)', fontsize=18, weight='bold', transform=ax[1].transAxes)
    plt.show()


def main():
    # plot_distributions(start, stop, step, 'carbon1', file_dir)
    # plot_distributions(start, stop, step, 'carbon2', file_dir)

    # targs = [(n, step, 'Density', 'electrons', file_dir) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #    p.starmap(plot_metric, targs)

    # targs = [(n, step, 'Density', 'ions', file_dir) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #    p.starmap(plot_metric, targs)

    file_dir = '/home/cepheid/TriForce/game_engine/data'
    start = 0
    stop = 2000
    step = 75
    targs = [(n, step, file_dir) for n in range(start, stop + step, step)]
    with mp.Pool(8) as p:
       p.starmap(plot_fields, targs)

    # targs = [(n, step, 'Jz', file_dir) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(plot_single_field, targs)

    # particle_positions(start, stop, step, 'electrons', file_dir)

    # plot_Temp(['electrons', 'ions'], start, stop, step, file_dir)
    # plot_KE(['electrons', 'ions'], start, stop, step, file_dir)
    # plot_field_energy(start, stop, step, file_dir)

    # plot_single_field(0, 1, 'Ex', file_dir)
    # plot_single_field(0, 1, 'Ez', file_dir)

    # fname = 'SB_G4_Z29_kdsdk_MeV_barns.csv'
    # oname = 'SB_G4_Z29_kdsdk_eV_m2.txt'
    # data = np.genfromtxt('/home/cepheid/TriForce/game_engine/tests/mikes_files/' + fname, dtype=np.float64)
    # data[:, 0] *= 1e6
    # data[:, 1:] *= 1e-28
    # np.savetxt('/home/cepheid/TriForce/game_engine/tests/cross_section_data/' + oname, data)


if __name__ == '__main__':
    main()
