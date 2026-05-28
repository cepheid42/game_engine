#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from adios2 import FileReader, Stream

from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'opposites'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

shape = (501, 2, 301)

xmin, xmax = -15.0e-6, -5.0e-6
ymin, ymax = 0.0, 0.01
zmin, zmax = -15.0e-6, 15.0e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 4.0e-17
t_end = 1.0e-14
nt = int(t_end / dt) + 1

save_interval = 5

# =====================
# ===== Particles =====
# =====================
py_range = (ymin, ymax)

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
    px_range=(-11.0e-6, -9.0e-6),
    py_range=py_range,
    pz_range=(-10.0e-6, 10.0e-6)
)

ions = Particles(
    name='positrons',
    mass=constants.m_e,
    charge=1,
    atomic_number=0,
    temp=temp,
    density=density,
    ppc=ppc,
    distribution='relativistic',
    px_range=(-11.0e-6, -9.0e-6),
    py_range=py_range,
    pz_range=(-10.0e-6, 10.0e-6)
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs=ParticleBCType.Outflow,
    push_type=ParticlePushType.Boris,
    interp_order=1,
    particle_data=(electrons, ions)
)

# ==================================
# ===== Electromagnetic Params =====
# ==================================
fields_path = data_path + f'/{sim_name}_applied_fields.bp'
em_params = EMParams(
    save_interval=save_interval,
    em_bcs=(2, 2, 2, 2, 2, 2),
    applied_fields=fields_path
)

Ex = np.full((shape[0] - 1, shape[1], shape[2]), 1.0e13)
Bz = np.full((shape[0] - 1, shape[1] - 1, shape[2]), 1.0e4)

with Stream(fields_path, 'w') as f:
    f.write('Ex', Ex, Ex.shape, (0, 0, 0), Ex.shape)
    f.write('Bz', Bz, Bz.shape, (0, 0, 0), Bz.shape)

# ==========================
# ===== Metrics Params =====
# ==========================
metric_params = Metrics(
    data_path,
    (
        MetricType.ParticleEnergy,
        MetricType.FieldEnergy,
        MetricType.FieldDump,
        MetricType.ParticleDump,
        MetricType.ParticleDiagnostics,
    )
)

# ============================
# ===== Simulation Class =====
# ============================
sim_params = Simulation(
    name=sim_name,
    shape=shape,
    nthreads=32,
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

# # ===========================
# # ===== Compile and Run =====
# # ===========================
# print(f'Setting up "{sim_name}"')
# create_data_dir(data_path)
# create_particles(sim_params, electrons, data_path)
# create_particles(sim_params, ions, data_path)
# update_header(sim_params, project_path=project_path)
#
# compile_project(build_path, output=True)
# run_project(build_path + '/game_engine', output=True)

# ===========================
# ===== Post Processing =====
# ===========================
xs = np.linspace(xmin, xmax, shape[0])
zs = np.linspace(zmin, zmax, shape[2])

# with FileReader(data_path + f'/electrons_{250:010d}.bp') as f:
#     density = f.read("Density")[:, 0, :].T
#     # temp1 = f.read('Temperature')[:, 0, :].T
#
#
# fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
# im = ax.contourf(xs[:-1], zs[:-1], density, levels=70, cmap='jet')
# plt.colorbar(im, ax=ax)
# plt.show()

Jxs = []
Jys = []
Jzs = []
time = 0
with FileReader(data_path + f'/fields_{250:010d}.bp') as f:
    # Jxs.append(f.read("Jx")[:, 0, :].T)
    # Jys.append(f.read("Jy")[:, 0, :].T)
    # Jzs.append(f.read("Jz")[:, 0, :].T)

    Jxs.append(f.read("Ex")[:, 0, :].T + Ex[:, 0, :].T)
    Jys.append(f.read("Ey")[:, 0, :].T)
    Jzs.append(f.read("Ez")[:, 0, :].T)

    # Jxs.append(f.read("Hx")[:, 0, :].T * constants.mu_0)
    # Jys.append(f.read("Hy")[:, 0, :].T * constants.mu_0)
    # Jzs.append(f.read("Hz")[:, 0, :].T * constants.mu_0 + Bz[:, 0, :].T)

    time = f.read("Time")

fig, ax = plt.subplots(3, 1, figsize=(6, 12), layout='constrained')

fig.suptitle(f't = {time * s_to_ns} ns')

a0 = ax[0].contourf(xs[:-1], zs, Jxs[-1], levels=70, cmap='jet')
# a0 = ax[0].contourf(xs, zs[:-1], Jxs[-1], levels=70, cmap='jet')
ax[0].set_xlabel('x')
ax[0].set_ylabel('z')
ax[0].set_title('Jx')

a1 = ax[1].contourf(xs, zs, Jys[-1], levels=70, cmap='jet')
# a1 = ax[1].contourf(xs[:-1], zs[:-1], Jys[-1], levels=70, cmap='jet')
ax[1].set_xlabel('x')
ax[1].set_ylabel('z')
ax[1].set_title('Jy')

a2 = ax[2].contourf(xs, zs[:-1], Jzs[-1], levels=70, cmap='jet')
# a2 = ax[2].contourf(xs[:-1], zs, Jzs[-1], levels=70, cmap='jet')

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

# Exs = []
# Eys = []
# Ezs = []
# Hxs = []
# Hys = []
# Hzs = []
# for n in range(0, nt + save_interval, save_interval):
#     with FileReader(data_path + f'/fields_{n:010d}.bp') as f:
#         Exs.append(np.sum(f.read('Ex')[:, 0, :]**2))
#         Eys.append(np.sum(f.read('Ey')[:, 0, :]**2))
#         Ezs.append(np.sum(f.read('Ez')[:, 0, :]**2))
#         Hxs.append(np.sum(f.read('Hx')[:, 0, :]**2))
#         Hys.append(np.sum(f.read('Hy')[:, 0, :]**2))
#         Hzs.append(np.sum(f.read('Hz')[:, 0, :]**2))
#
# Exs = np.array(Exs)
# Eys = np.array(Eys)
# Ezs = np.array(Ezs)
# Hxs = np.array(Hxs)
# Hys = np.array(Hys)
# Hzs = np.array(Hzs)
#
# ex = 0.5 * constants.epsilon_0 * Exs * dx * dy * dz
# ey = 0.5 * constants.epsilon_0 * Eys * dx * dy * dz
# ez = 0.5 * constants.epsilon_0 * Ezs * dx * dy * dz
# bx = 0.5 * constants.mu_0 * Hxs * dx * dy * dz
# by = 0.5 * constants.mu_0 * Hys * dx * dy * dz
# bz = 0.5 * constants.mu_0 * Hzs * dx * dy * dz
#
# with FileReader(data_path + '/particles_energy.bp') as f:
#     variables = f.available_variables()
#     steps = int(variables['Time']['AvailableStepsCount'])
#     time = f.read('Time', step_selection=[0, steps])
#     e_energy = f.read('electrons', step_selection=[0, steps])
#     i_energy = f.read('positrons', step_selection=[0, steps])
#
# time *= s_to_fs
# field_energy = (ex + ey + ez + bx + by + bz) * m_to_cm / dy
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
