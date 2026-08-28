import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from adios2 import FileReader
from matplotlib.cm import ScalarMappable
from matplotlib import colors, ticker
from pathlib import Path

to_kilo = 1.0e-3
to_mega = 1.0e-6
s_to_fs = 1.0e15
s_to_ns = 1.0e9
Vm_to_kVcm = 1.0e-5
T_to_gauss = 1.0e4
Am2_to_Acm2 = 1.0e-4
m_to_cm = 1.0e-2
m3_to_cm3 = 1.0e-6

def plot_density(groups, step, data_path, block=True, save=False):
    from matplotlib import gridspec

    if isinstance(groups, str):
        groups = [groups]

    cols = 2 if len(groups) != 1 else 1
    rows = len(groups) // cols
    if len(groups) % cols != 0:
        rows += 1

    fig = plt.figure(figsize=(6 * cols, 6 * rows), layout='constrained')
    gs = gridspec.GridSpec(rows, cols, figure=fig)

    for i, g in enumerate(groups):
        filename = f'/{g}_{step:010d}.bp'
        if not Path(data_path + filename).exists():
            continue

        with FileReader(data_path + filename) as f:
            density = f.read("Density")[::-1, 0, :] * m3_to_cm3
            times = f.read('Time') * s_to_ns
        density = np.ma.masked_where(density <= 0, density)

        ax = fig.add_subplot(gs[i])
        ax.set_title(f'{g.capitalize()} Density')

        # cmap = 'plasma' if g == 'protons' else 'viridis'
        # den_norm = colors.LogNorm(vmin=1e18, vmax=1e22)
        # im = ax.contourf(zs[:-1], xs[:-1], density, levels=np.logspace(18, 22, 50), cmap=cmap, norm=den_norm)
        # cbar = plt.colorbar(ScalarMappable(norm=den_norm, cmap=cmap), ax=ax)

        im = ax.contourf(density, levels=100, cmap='jet')

        cbar = plt.colorbar(im, ax=ax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(r'$\text{cm}^{-1}$')

    # fig.suptitle(f'nt = {times:.7f} ns')
    if save:
        plt.savefig(data_path + '/density.png')
    else:
        plt.show(block=block)


def plot_temperature(groups, step, data_path, xs, zs, block=True, save=False):
    from matplotlib import gridspec

    if isinstance(groups, str):
        groups = [groups]

    cols = 2 if len(groups) != 1 else 1
    rows = len(groups) // cols
    if len(groups) % cols != 0:
        rows += 1

    fig = plt.figure(figsize=(6 * cols, 6 * rows), layout='constrained')
    gs = gridspec.GridSpec(rows, cols, figure=fig)

    for i, g in enumerate(groups):
        filename = f'/{g}_{step:010d}.bp'
        with FileReader(data_path + filename) as f:
            temperature = f.read("Temperature")[:, 0, :].T
            times = f.read('Time') * s_to_ns
        temperature = np.ma.masked_where(temperature <= 0, temperature)

        ax = fig.add_subplot(gs[i])
        ax.set_title(f'{g.capitalize()} Temperature')

        # if g == 'deuterium':
        #     temp_norm = colors.LogNorm(vmin=10, vmax=1e5)
        #     im = ax.contourf(xs[:-1], zs[:-1], temperature, levels=np.logspace(1, 5, 50), cmap='jet', norm=temp_norm)
        #     cbar = plt.colorbar(ScalarMappable(norm=temp_norm, cmap='jet'), ax=ax)
        # elif g == 'electrons':
        #     temp_norm = colors.LogNorm(vmin=10, vmax=1e5)
        #     im = ax.contourf(xs[:-1], zs[:-1], temperature, levels=np.logspace(1, 5, 50), cmap='jet', norm=temp_norm)
        #     cbar = plt.colorbar(ScalarMappable(norm=temp_norm, cmap='jet'), ax=ax)
        # else:
        im = ax.contourf(xs[:-1], zs[:-1], temperature, levels=100, cmap='jet')
        cbar = plt.colorbar(im, ax=ax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('eV')

    fig.suptitle(f'nt = {times:.7f} ns')
    if save:
        plt.savefig(data_path + '/temperature.png')
    else:
        plt.show(block=block)

    # temp_norm = colors.LogNorm(vmin=1, vmax=1e4)
    # im = ax.contourf(xs[:-1], zs[:-1], temp, levels=np.logspace(1, 4.5, 50), cmap='jet', norm=temp_norm)
    # cbar = plt.colorbar(ScalarMappable(norm=temp_norm, cmap='jet'), ax=ax)

    # temp_norm = colors.Normalize(vmin=1.0, vmax=10000)
    # im = ax.contourf(xs[:-1], zs[:-1], temp, levels=50, cmap='jet', norm=temp_norm)
    # cbar = plt.colorbar(ScalarMappable(norm=temp_norm, cmap='jet'), ax=ax)

    # im = ax.contourf(xs[:-1], zs[:-1], temp, locator=ticker.LogLocator(), cmap='jet')
    # cbar = plt.colorbar(im, ax=ax)

    # im = ax.contourf(xs[:-1], zs[:-1], temp, levels=100, cmap='jet')
    # cbar = plt.colorbar(im, ax=ax)


def plot_field_energy(data_path, smith_data, block=True, save=False):
    with FileReader(data_path + '/fields_energy.bp') as f:
        variables = f.available_variables()
        steps = int(variables['Time']['AvailableStepsCount'])
        time = f.read('Time', step_selection=[0, steps])
        ex = f.read('Ex Energy', step_selection=[0, steps])
        ey = f.read('Ey Energy', step_selection=[0, steps])
        ez = f.read('Ez Energy', step_selection=[0, steps])
        bx = f.read('Bx Energy', step_selection=[0, steps])
        by = f.read('By Energy', step_selection=[0, steps])
        bz = f.read('Bz Energy', step_selection=[0, steps])

    time *= s_to_fs
    field_energy = (ex + ey + ez + bx + by + bz) * to_kilo / 0.01 # kJ/dy # Just joules
    fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(r'Joules')
    ax.set_title('Field Energy')
    ax.plot(time, field_energy, label='TriForce')

    smith_field_data = np.genfromtxt(smith_data, delimiter=',')
    ax.plot(smith_field_data[:, 0], smith_field_data[:, 1], 'r--', label='Smith')

    if save:
        plt.savefig(data_path + '/field_energy.png')
    else:
        plt.show(block=block)


def plot_particle_energy(data_path, smith_data, block=True, save=False):
    particles = dict()
    with FileReader(data_path + '/particles_energy.bp') as f:
        variables = f.available_variables()
        steps = int(variables['Time']['AvailableStepsCount'])
        time = f.read('Time', step_selection=[0, steps])
        for v in variables.keys():
            if v in ['Step', 'Time', 'dt']:
                continue
            particles[v] = f.read(v, step_selection=[0, steps])

    time *= s_to_fs

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(r'Joules')
    ax.set_title('Particle Energy')

    for k, v in particles.items():
        ax.plot(time, v * to_kilo / 0.01, label=f'{k.capitalize()} Energy')

    electron_data = np.genfromtxt(smith_data + f'/smith_lsi_electron_energy.csv', delimiter=',')
    proton_data = np.genfromtxt(smith_data + f'/smith_lsi_proton_energy.csv', delimiter=',')

    ax.plot(electron_data[:, 0], electron_data[:, 1], 'r--', label='Smith Electron')
    ax.plot(proton_data[:, 0], proton_data[:, 1], 'g--', label='Smith Ion')

    ax.legend()
    if save:
        plt.savefig(data_path + '/particle_energy.png')
    else:
        plt.show(block=block)



def plot_total_particle_yield(data_path, groups, bounds, block=True):
    start, stop, step = bounds

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Y')
    ax.set_title('Yield')

    plot_params = {
        'electrons': ('--', 'k', 'x', 8, 'full', 25),
        'photons':   ('--', 'k', '>', 8, 'full', 20),
        'protons':   ('--', 'r', 'P', 8, 'full', 20), # proton
        'neutrons':  ('--', 'b', 'D', 8, 'none', 25), # neutron
        'helium3':   ('-.', 'g', 'v', 8, 'full', 20), # He3
        'tritium':   ('-.', 'm', '*', 8, 'none', 25), # tritium
    }

    # particles = dict()
    for g in groups:
        density = []
        times = []
        for n in range(start, stop + step, step):
            filename = f'/{g}_{n:010d}.bp'
            with FileReader(data_path + filename) as f:
                cell_vol = f.read_attribute('Cell Volume')
                density.append(f.read('Density').sum() * cell_vol)
                times.append(f.read('Time'))
        density = np.array(density)
        times = np.array(times)

        _, c, m, ms, fs, mark_every = plot_params[g]
        ax.plot(times * s_to_ns, density, c=c, marker=m, ms=ms, markevery=mark_every, fillstyle=fs, label=g.capitalize())

    ax.legend()
    plt.show()


def plot_fields(data_path, step):
    filename = f'/fields_{step:010d}.bp'
    with FileReader(data_path + filename) as f:
        Ex = f.read('Ey')[:, 0, :].T

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
    im = ax.contourf(Ex)
    plt.colorbar(im, ax=ax)
    plt.show()


def plot_energies(data_paths, smith_path):
    lsp_energies = np.genfromtxt(smith_path + '/Energy_Diagnostics/LSP_E/history.txt', comments='#', usecols=(1, 3, 5, 6))
    lsp_energies[:, 0] *= 1.0e6 # time ns -> fs
    
    def load_energies(ipath):
        # Read EM Field Energies
        with FileReader(ipath + '/fields_energy.bp') as f:
            varbles = f.available_variables()
            nsteps = int(varbles['Time']['AvailableStepsCount'])
            ex = f.read('Ex Energy', step_selection=[0, nsteps])
            ey = f.read('Ey Energy', step_selection=[0, nsteps])
            ez = f.read('Ez Energy', step_selection=[0, nsteps])
            bx = f.read('Bx Energy', step_selection=[0, nsteps])
            by = f.read('By Energy', step_selection=[0, nsteps])
            bz = f.read('Bz Energy', step_selection=[0, nsteps])
        # Read Particle Energies
        with FileReader(ipath + '/particles_energy.bp') as f:
            variables = f.available_variables()
            steps = int(variables['Time']['AvailableStepsCount'])
            t = f.read('Time', step_selection=[0, steps])
            e = f.read('electrons', step_selection=[0, steps])
            ii = f.read('protons', step_selection=[0, steps])
            # ii = f.read('deuterium', step_selection=[0, steps])

        t *= s_to_fs
        return (ex + ey + ez + bx + by + bz), e, ii, t

    data = dict()
    for path in data_paths:
        name = ' '.join(p.capitalize() for p in Path(path).name.split('_')[1:])
        data[name] = load_energies(path)

    fig, ax = plt.subplots(3, 1, figsize=(6, 10), layout='constrained')
    ax[0].plot(lsp_energies[:, 0], lsp_energies[:, 2], 'r--', label='LSP(E)')
    ax[1].plot(lsp_energies[:, 0], lsp_energies[:, 1], 'r--', label='LSP(E)')
    ax[2].plot(lsp_energies[:, 0], lsp_energies[:, 3], 'r--', label='LSP(E)')
    for label, (fields, e_energy, i_energy, time) in data.items():
        ax[0].plot(time, fields, '-.', label=label)
        ax[1].plot(time, e_energy, '-.', label=label)
        ax[2].plot(time, i_energy, '-.', label=label)

    y_labels = [r'Field (J cm$^{-1}$)', r'Electron (J cm$^{-1}$)',  r'Proton (J cm$^{-1}$)']
    for i, a in enumerate(ax):
        a.grid()
        a.set_xlabel('Time (fs)')
        a.set_ylabel(y_labels[i])
        a.legend()

    plt.show()


def plot_spectra(data_paths, smith_path):
    import os
    lsp_spectra = np.genfromtxt(smith_path + '/Particle_Spectra/LSPE_SPEC/30_Spectrum.dat')
    lsp_spectra[:, 1:] *= 100.0

    def get_spectra(pathy, group):
        step = 15000 if os.path.exists(pathy + f'/{group}_dump_{15000:010d}.bp') else 7500
        with FileReader(pathy + f'/{group}_dump_{step:010d}.bp') as f:
            mass = f.read_attribute('Mass')
            charge = f.read_attribute('Charge')
            weight = f.read('Weight')
            gamma = f.read('Gamma')
            velocity = f.read('Velocity')

        v_mask = np.where(velocity[:, 0] > 0.0) # only particles moving in +x
        weight = weight[v_mask]
        gamma = gamma[v_mask]

        energies = (gamma - 1.0) * mass * constants.c**2 / constants.e
        energies *= to_mega

        limit = 10.0 if group == 'electrons' else 22.5
        energy_bins = np.arange(0, limit, 0.1) # bins in MeV
        count, bins = np.histogram(energies, bins=energy_bins, weights=weight)

        charge = np.abs(charge) * count / 0.01
        return charge, energy_bins

    data = dict()
    for path in data_paths:
        name = ' '.join(p.capitalize() for p in Path(path).name.split('_')[1:])
        i_spectra, i_bins = get_spectra(path, 'protons')
        # i_spectra, i_bins = get_spectra(path, 'deuterium')
        e_spectra, e_bins = get_spectra(path, 'electrons')
        data[name] = (i_spectra, i_bins, e_spectra, e_bins)

    fig, ax = plt.subplots(2, 1, figsize=(6, 10), layout='constrained')
    ax[0].plot(lsp_spectra[:100, 0], lsp_spectra[:100, 2], label='LSP(E)')
    ax[1].plot(lsp_spectra[:209, 0], lsp_spectra[:209, 1], label='LSP(E)')

    for label, (i, ib, e, eb) in data.items():
        ax[0].plot(eb[:-1], e, label=label)
        ax[1].plot(ib[:-1], i, label=label)

    for a in ax:
        a.set_yscale('log')
        a.grid()
        a.legend()

    plt.show()
