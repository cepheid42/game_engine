#!/usr/bin/env python3

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from adios2 import FileReader, Stream
from scipy import constants

from scripts.particle_generation import create_particles
from scripts.domain_params import *

sim_name = 'seinfeld3D'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + '/data'

shape = (57, 57, 113)

xmin, xmax = -0.14, 0.14
ymin, ymax = -0.14, 0.14
zmin, zmax = -0.28, 0.28

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 5.0e-12
t_end = 80e-9
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = 20

# ===== Particles =====
px_range = (-4e-2, 4e-2)
py_range = (-4e-2, 4e-2)
pz_range = (-15e-2, 15e-2)
e_temp = 4.0 / np.sqrt(3)
i_temp = 1.0 / np.sqrt(3)

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=(e_temp, e_temp, e_temp), # eV
    density=1.0e17, # m^-3,
    ppc=(3, 3, 3),
    # distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

ions = Particles(
    name='ions',
    mass=constants.m_p,
    charge=1,
    atomic_number=1,
    temp=(i_temp, i_temp, i_temp), # eV
    density=1.0e17, # m^-3,
    ppc=(2, 2, 2),
    # distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ===== Collisions and Particle Params =====
particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs='outflow',
    interp_order=1,
    particle_data=(electrons, ions)
)

# ===== Electromagnetic Params =====
Bz_applied = np.full((shape[0] - 1, shape[1] - 1, shape[2]), 0.1)
with Stream(data_path + f'/{sim_name}_applied_fields.bp', 'w') as f:
    f.write('Bz', Bz_applied, Bz_applied.shape, (0, 0, 0), Bz_applied.shape)

em_params = EMParams(
    save_interval=save_interval,
    pml_depth=6,
    em_bcs=(1, 1, 1, 1, 1, 1),
    applied_fields=data_path + f'/{sim_name}_applied_fields.bp'
)

# ===== Metrics Params =====
metric_params = Metrics(
    data_path,
    (MetricType.ParticleEnergy, MetricType.FieldEnergy)
)

sim_params = Simulation(
    name=sim_name,
    shape=shape,
    nthreads=48,
    dt=dt,
    t_end=t_end,
    nt=nt,
    cfl=cfl,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    em_params=em_params,
    particle_params=particle_params,
    metric_params=metric_params,
    coll_enabled=False
)

print(f'Setting up "{sim_name}"')
create_particles(sim_params, electrons, data_path)
create_particles(sim_params, ions, data_path)
update_header(sim_params, project_path=project_path)

subprocess.run(
    ['meson', 'compile', '-C', build_path, '-j4'],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL
).check_returncode()

subprocess.run(build_path + '/game_engine').check_returncode()

with FileReader(data_path + f'/{sim_name}/fields_energy.bp') as f:
    variables = f.available_variables()
    steps = int(variables['Time']['AvailableStepsCount'])
    ex = f.read('Ex Energy', step_selection=[0, steps])
    ey = f.read('Ey Energy', step_selection=[0, steps])
    ez = f.read('Ez Energy', step_selection=[0, steps])
    bx = f.read('Bx Energy', step_selection=[0, steps])
    by = f.read('By Energy', step_selection=[0, steps])
    bz = f.read('Bz Energy', step_selection=[0, steps])

field_energy = ex + ey + ez + bx + by + bz

with FileReader(data_path + f'/{sim_name}/particles_energy.bp') as f:
    variables = f.available_variables()
    steps = int(variables['Time']['AvailableStepsCount'])
    time = f.read('Time', step_selection=[0, steps]) * 1e9
    e_energy = f.read('electrons', step_selection=[0, steps])
    i_energy = f.read('ions', step_selection=[0, steps])

particle_energy = e_energy + i_energy

lsp_particle_data= np.genfromtxt('./data/seinfeld-3d-lsp-particles-energy.txt', skip_header=1)
lsp_particles = np.interp(time, lsp_particle_data[:, 0], lsp_particle_data[:, 1])

lsp_field_data = np.genfromtxt('./data/seinfeld-3d-lsp-fields-energy.txt', skip_header=1)
lsp_fields = np.interp(time, lsp_field_data[:, 0], lsp_field_data[:, 1])

fig, ax1 = plt.subplots(figsize=(10, 6), layout='constrained')
ax1.set_xlim([0.0, 80.0])
ax1.grid()
ax2 = ax1.twinx()

ax1.plot(time, lsp_particles, 'r--', label='LSP Particles')
ax1.plot(time, particle_energy, color='r', label='Particle KE')

ax2.plot(time, lsp_fields, 'b--', label='LSP Fields')
ax2.plot(time, field_energy, color='b', label='Field Energy')

line1, label1 = ax1.get_legend_handles_labels()
line2, label2 = ax2.get_legend_handles_labels()

ax1.set_ylim([2.27e-4, 2.31e-4])
ax1.set_yticks(np.arange(0.000226, 0.000231, 0.0000005))
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('Particle KE (J)')
ax1.set_title(f'Total Energy Seinfeld 3D')

ax2.set_ylim([0.0, 3.25e-6])
ax2.set_yticks(np.arange(0.0, 3.2e-6, 1.0e-6))
ax2.set_ylabel('Field Energy (J)')
ax2.legend(line1 + line2, label1 + label2)

plt.show()
