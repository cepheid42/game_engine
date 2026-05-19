#!/usr/bin/env python3

import subprocess
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from adios2 import FileReader

from scripts.particle_generation import create_particles
from scripts.domain_params import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'lsi_test'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

shape = (551, 2, 351)

xmin, xmax = -15.5e-6, -4.5e-6
ymin, ymax = 0.0, 0.01
zmin, zmax = -15.5e-6, 15.5e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 4.0e-17
t_end = 3.0e-14
nt = int(t_end / dt) + 1

save_interval = 5

# =====================
# ===== Particles =====
# =====================
px_range = (-13.0e-6, -7.0e-6) # meters
py_range = (ymin, ymax)
pz_range = (-10.0e-6, 10.0e-6)

ppc = (10, 1, 10)
density = 8.5e25
temp = tuple(3 * [10000 / np.sqrt(3)])

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=temp,
    density=density,
    ppc=ppc,
    distribution='relativistic',
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
    distribution='relativistic',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs='outflow',
    interp_order=1,
    particle_data=(electrons, ions)
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
        MetricType.ParticleEnergy,
        MetricType.FieldEnergy,
        MetricType.FieldDump
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
    collisions_enabled=False,
    # push_enabled=False,
    # jdep_enabled=False
)

# ===========================
# ===== Compile and Run =====
# ===========================
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

# ===========================
# ===== Post Processing =====
# ===========================
J_to_kJ = 1.0e-3
s_to_fs = 1.0e15
s_to_ns = 1.0e9
Vm_to_kVcm = 1.0e-5
T_to_gauss = 1.0e4
Am2_to_Acm2 = 1.0e-4

Jxs = []
Jys = []
Jzs = []
time = 0
with FileReader(data_path + f'/fields_{750:010d}.bp') as f:
    Jxs.append(f.read("Jx")[:, 0, :].T)
    Jys.append(f.read("Jy")[:, 0, :].T)
    Jzs.append(f.read("Jz")[:, 0, :].T)
    time = f.read("Time")

fig, ax = plt.subplots(3, 1, figsize=(6, 12), layout='constrained')

fig.suptitle(f't = {time * s_to_ns} ns')

a0 = ax[0].contourf(Jxs[-1], levels=70, cmap='jet')
ax[0].set_xlabel('x')
ax[0].set_ylabel('z')
ax[0].set_title('Jx')

a1 = ax[1].contourf(Jys[-1], levels=70, cmap='jet')
ax[1].set_xlabel('x')
ax[1].set_ylabel('z')
ax[1].set_title('Jy')

a2 = ax[2].contourf(Jzs[-1], levels=70, cmap='jet')
ax[2].set_xlabel('x')
ax[2].set_ylabel('z')
ax[2].set_title('Jz')

plt.colorbar(a0, ax=ax[0])
plt.colorbar(a1, ax=ax[1])
plt.colorbar(a2, ax=ax[2])

plt.show()

# with FileReader(data_path + '/fields_energy.bp') as f:
#     variables = f.available_variables()
#     steps = int(variables['Time']['AvailableStepsCount'])
#     ex = f.read('Ex Energy', step_selection=[0, steps])
#     ey = f.read('Ey Energy', step_selection=[0, steps])
#     ez = f.read('Ez Energy', step_selection=[0, steps])
#     bx = f.read('Bx Energy', step_selection=[0, steps])
#     by = f.read('By Energy', step_selection=[0, steps])
#     bz = f.read('Bz Energy', step_selection=[0, steps])
#
# with FileReader(data_path + '/particles_energy.bp') as f:
#     variables = f.available_variables()
#     steps = int(variables['Time']['AvailableStepsCount'])
#     time = f.read('Time', step_selection=[0, steps])
#     e_energy = f.read('electrons', step_selection=[0, steps])
#     i_energy = f.read('ions', step_selection=[0, steps])
#
# time *= s_to_fs
# field_energy = (ex + ey + ez + bx + by + bz) * J_to_kJ / dy
# e_energy = e_energy * J_to_kJ / dy
# i_energy = i_energy * J_to_kJ / dy
#
# y_labels = [
#     (r'Field (kJ m$^{-1}$)', [0, 100]),
#     (r'Electron (kJ m$^{-1}$)', [0, 15]),
#     (r'Proton (kJ m$^{-1}$)', [0, 15])
# ]
#
# fig, ax = plt.subplots(3, 1, figsize=(6, 10), layout='constrained')
#
# for i, a in enumerate(ax):
#     label, _ = y_labels[i]
#     a.grid()
#     a.set_xlabel('Time (fs)')
#     a.set_ylabel(label)
#
# a0 = ax[0].plot(time, field_energy, 'b-', label='TriForce')
# a1 = ax[1].plot(time, e_energy, 'b-')
# a2 = ax[2].plot(time, i_energy, 'b-')
#
# # plt.savefig(data_path + f'/lsi_comp_normal.png')
# # plt.close(fig)
#
# plt.show()
