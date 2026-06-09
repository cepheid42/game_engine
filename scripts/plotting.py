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


def plot_density(group, step, data_path, xs, zs, block=True):
    filename = f'/{group}_{step:010d}.bp'
    with FileReader(data_path + filename) as f:
        density = f.read("Density")[:, 0, :].T * inv_cubic_m_to_cubic_cm
        time = f.read('Time') * s_to_ns

    density = np.ma.masked_where(density <= 0, density)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title(f'{group.capitalize()} Density\nt = {time:.7f} ns')

    # den_norm = colors.LogNorm(vmin=1e21, vmax=2e23)
    # im = ax.contourf(xs[:-1], zs[:-1], density, levels=np.logspace(21, 24, 50), cmap='jet', norm=den_norm)
    # cbar = plt.colorbar(ScalarMappable(norm=den_norm, cmap='jet'), ax=ax)

    den_norm = colors.Normalize(vmin=4.848e21, vmax=1.805e23)
    im = ax.contourf(xs[:-1], zs[:-1], density, levels=100, cmap='jet', norm=den_norm)
    cbar = plt.colorbar(ScalarMappable(norm=den_norm, cmap='jet'), ax=ax)

    # im = ax.contourf(xs[:-1], zs[:-1], density, levels=100, cmap='jet')
    # cbar = plt.colorbar(im, ax=ax)

    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel(r'$\text{cm}^{-1}$')

    plt.show(block=block)


def plot_temperature(group, step, data_path, xs, zs, block=True, save=False):
    filename = f'/{group}_{step:010d}.bp'
    with FileReader(data_path + filename) as f:
        temp = f.read("Temperature")[:, 0, :].T
        time = f.read('Time') * s_to_ns

    temp = np.ma.masked_where(temp <= 0, temp)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title(f'{group.capitalize()} Temperature\nt = {time:.7f} ns')

    temp_norm = colors.LogNorm(vmin=1, vmax=1e4)
    im = ax.contourf(xs[:-1], zs[:-1], temp, levels=np.logspace(1, 4.5, 50), cmap='jet', norm=temp_norm)
    cbar = plt.colorbar(ScalarMappable(norm=temp_norm, cmap='jet'), ax=ax)

    # temp_norm = colors.Normalize(vmin=1.0, vmax=10000)
    # im = ax.contourf(xs[:-1], zs[:-1], temp, levels=50, cmap='jet', norm=temp_norm)
    # cbar = plt.colorbar(ScalarMappable(norm=temp_norm, cmap='jet'), ax=ax)

    # im = ax.contourf(xs[:-1], zs[:-1], temp, locator=ticker.LogLocator(), cmap='jet')
    # cbar = plt.colorbar(im, ax=ax)

    # im = ax.contourf(xs[:-1], zs[:-1], temp, levels=100, cmap='jet')
    # cbar = plt.colorbar(im, ax=ax)

    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('eV')

    if save:
        plt.savefig(data_path + f'/{group}_temperature_{step}.png')
        plt.close(fig)
    else:
        plt.show(block=block)

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

def plot_total_neutron_yield(data_path, bounds, block=True):
    start, stop, step = bounds

    times = []
    density = []
    for n in range(start, stop + step, step):
        filename = f'/neutrons_{n:010d}.bp'
        with FileReader(data_path + filename) as f:
            density.append(f.read('Density').sum())
            times.append(f.read('Time') * s_to_ns)

    times = np.array(times)
    density = np.array(density)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Y')
    ax.set_title('Neutron Yield')

    ax.plot(times, density)
    plt.show()



