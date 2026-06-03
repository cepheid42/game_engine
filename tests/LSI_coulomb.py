#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from adios2 import FileReader

from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================

sim_name = f'lsi_coulomb_10eV_1e29m2_coll10steps'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

shape = (1551, 2, 1551)

xmin, xmax = -15.5e-6, 15.5e-6
ymin, ymax = 0.0, 0.01
zmin, zmax = -15.5e-6, 15.5e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 4.0e-17
t_end = 3.0e-13
nt = int(t_end / dt) + 1

save_interval = 10

# =====================
# ===== Particles =====
# =====================
px_range = (-5e-7, 5e-7) # meters
py_range = (ymin, ymax)
pz_range = (-1e-5, 1e-5)

ppc = (5, 1, 5)
density = 1.0e29 #8.5e27
temp = tuple(3 * [10 / np.sqrt(3)])

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=temp,
    density=density,
    ppc=ppc,
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

ions = Particles(
    name='ions',
    mass=constants.m_p,
    charge=1,
    atomic_number=1,
    temp=temp,
    density=density,
    ppc=ppc,
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
coulombParam = CoulombParams(10, 1) # set LnLambda = 0 to calculate it on the fly

particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs=ParticleBCType.Outflow,
    push_type=ParticlePushType.Boris,
    interp_order=1,
    particle_data=(electrons, ions),
    collisions=(
        Collision(groups=(electrons, ions), channels=('coulomb',), step_interval=10),
        Collision(groups=(ions, ions), channels=('coulomb',), self_scatter=True, step_interval=10),
        Collision(groups=(electrons, electrons), channels=('coulomb',), self_scatter=True, step_interval=10),
    )
)

# ==================================
# ===== Electromagnetic Params =====
# ==================================
em_params = EMParams(
    save_interval=save_interval,
    pml_depth=15,
    em_bcs=(1, 1, 2, 2, 1, 1),
    laser_enabled=True,
)

# ==========================
# ===== Metrics Params =====
# ==========================
metric_params = Metrics(
    data_path,
    (
        # MetricType.ParticleEnergy,
        # MetricType.FieldEnergy,
        # MetricType.FieldDump,
        # MetricType.ParticleDump,
        MetricType.ParticleDiagnostics,
    )
)

# ============================
# ===== Simulation Class =====
# ============================
sim_params = Simulation(
    name=sim_name,
    shape=shape,
    nthreads=48,
    dt=dt,
    t_end=t_end,
    nt=nt,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    em_params=em_params,
    particle_params=particle_params,
    metric_params=metric_params,
    # collisions_enabled=False,
    # push_enabled=False,
    # jdep_enabled=False
)

# ===========================
# ===== Compile and Run =====
# ===========================
print(f'Setting up "{sim_name}"')
create_data_dir(data_path)
create_particles(sim_params, electrons, data_path)
create_particles(sim_params, ions, data_path)
update_header(sim_params, project_path=project_path)

compile_project(build_path, output=True)
run_project(build_path + '/game_engine', output=True)

# # ===========================
# # ===== Post Processing =====
# # ===========================
with FileReader(data_path + '/fields_energy.bp') as f:
    variables = f.available_variables()
    steps = int(variables['Time']['AvailableStepsCount'])
    ex = f.read('Ex Energy', step_selection=[0, steps])
    ey = f.read('Ey Energy', step_selection=[0, steps])
    ez = f.read('Ez Energy', step_selection=[0, steps])
    bx = f.read('Bx Energy', step_selection=[0, steps])
    by = f.read('By Energy', step_selection=[0, steps])
    bz = f.read('Bz Energy', step_selection=[0, steps])

with FileReader(data_path + '/particles_energy.bp') as f:
    variables = f.available_variables()
    steps = int(variables['Time']['AvailableStepsCount'])
    time = f.read('Time', step_selection=[0, steps])
    e_energy = f.read('electrons', step_selection=[0, steps])
    i_energy = f.read('ions', step_selection=[0, steps])

time *= s_to_fs
field_energy = (ex + ey + ez + bx + by + bz) * J_to_kJ / dy
e_energy = e_energy * J_to_kJ / dy
i_energy = i_energy * J_to_kJ / dy

smith_field_data = np.genfromtxt('./data/smith_lsi_field_energy.csv', delimiter=',')
smith_electron_data = np.genfromtxt('./data/smith_lsi_electron_energy.csv', delimiter=',')
smith_proton_data = np.genfromtxt('./data/smith_lsi_proton_energy.csv', delimiter=',')

y_labels = [
    (r'Field (kJ m$^{-1}$)', smith_field_data, [0, 100]),
    (r'Electron (kJ m$^{-1}$)', smith_electron_data, [0, 15]),
    (r'Proton (kJ m$^{-1}$)', smith_proton_data, [0, 15])
]

fig, ax = plt.subplots(3, 1, figsize=(6, 10), layout='constrained')

for i, a in enumerate(ax):
    label, _, _ = y_labels[i]
    a.grid()
    a.set_xlabel('Time (fs)')
    a.set_ylabel(label)

ax[0].plot(time, field_energy, 'b-', label='TriForce')
ax[1].plot(time, e_energy, 'b-')
ax[2].plot(time, i_energy, 'b-')

ax[0].plot(smith_field_data[:, 0], smith_field_data[:, 1], 'r--', label='Smith')
ax[1].plot(smith_electron_data[:, 0], smith_electron_data[:, 1], 'r--')
ax[2].plot(smith_proton_data[:, 0], smith_proton_data[:, 1], 'r--')

# ax[0].set_ylim([0, 100])
# ax[1].set_ylim([0, 15])
# ax[2].set_ylim([0, 15])

ax[0].legend()

# plt.savefig(data_path + f'/lsi_comp_normal.png')
# plt.close(fig)

plt.show()
