import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from adios2 import FileReader
from dataclasses import dataclass, field
from matplotlib.cm import ScalarMappable
from matplotlib import colors, ticker

J_to_kJ = 1.0e-3
s_to_fs = 1.0e15
s_to_ns = 1.0e9
Vm_to_kVcm = 1.0e-5
T_to_gauss = 1.0e4
Am2_to_Acm2 = 1.0e-4
m_to_cm = 1.0e-2
inv_cubic_m_to_cubic_cm = 1.0e-6

@dataclass
class ParticlePlotData:
    velocities : np.ndarray = field(default_factory=lambda: np.zeros(shape=(1, 3), dtype=np.float64))
    positions : np.ndarray = field(default_factory=lambda: np.zeros(shape=(1, 3), dtype=np.float64))
    weights : np.ndarray = field(default_factory=lambda: np.zeros(shape=(1,), dtype=np.float64))
    gammas : np.ndarray = field(default_factory=lambda: np.zeros(shape=(1,), dtype=np.float64))
    times : np.ndarray = field(default_factory=lambda: np.zeros(shape=(1,), dtype=np.float64))


def plot_density(groups, step, data_path, xs, zs, block=True):
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
            density = f.read("Density")[:, 0, :].T * inv_cubic_m_to_cubic_cm
            times = f.read('Time') * s_to_ns
        density = np.ma.masked_where(density <= 0, density)

        ax = fig.add_subplot(gs[i])
        ax.set_title(f'{g.capitalize()} Density')

        im = ax.contourf(xs[:-1], zs[:-1], density, levels=100, cmap='jet')

        cbar = plt.colorbar(im, ax=ax)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(r'$\text{cm}^{-1}$')

    fig.suptitle(f'nt = {times:.7f} ns')
    plt.show(block=block)

    # den_norm = colors.LogNorm(vmin=1e21, vmax=2e23)
    # im = ax.contourf(xs[:-1], zs[:-1], density, levels=np.logspace(21, 24, 50), cmap='jet', norm=den_norm)
    # cbar = plt.colorbar(ScalarMappable(norm=den_norm, cmap='jet'), ax=ax)

    # den_norm = colors.Normalize(vmin=4.848e21, vmax=1.805e23)
    # im = ax.contourf(xs[:-1], zs[:-1], density, levels=100, cmap='jet', norm=den_norm)
    # cbar = plt.colorbar(ScalarMappable(norm=den_norm, cmap='jet'), ax=ax)

    # im = ax.contourf(xs[:-1], zs[:-1], density, levels=100, cmap='jet')
    # cbar = plt.colorbar(im, ax=ax)


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


def plot_field_energy(data_path, block=True):
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

    time *= s_to_ns
    field_energy = (ex + ey + ez + bx + by + bz) # Just joules
    fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(r'Joules')
    ax.set_title('Field Energy')
    ax.plot(time, field_energy)
    plt.show(block=block)


def plot_particle_energy(data_path, block=True):
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
        ax.plot(time, v, label=f'{k.capitalize()} Energy')

    ax.legend()
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







