#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from adios2 import FileReader

from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'adams_lsi_coulomb'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

shape = (551, 2, 151)

xmin, xmax = -15.5e-6, -4.5e-6
ymin, ymax = 0.0, 0.01
zmin, zmax = -15.5e-6, 15.5e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 2.0e-17
t_end = 1.5e-13
nt = int(t_end / dt) + 1

save_interval = 10

# =====================
# ===== Particles =====
# =====================
px_range = (-11.0e-6, -9.0e-6) # meters
py_range = (ymin, ymax)
pz_range = (-10.0e-6, 10.0e-6)

ppc = (5, 1, 5)
density = 1.0e29 #m^-3
Tcold = 1
Thot = 10

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=tuple(3 * [Thot]),
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
    temp=tuple(3 * [Tcold]),
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
coll_interval = 1
coulombParam = CoulombParams(0.0, 1) # set LnLambda = 0 to calculate it on the fly

particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs=ParticleBCType.Outflow,
    push_type=ParticlePushType.Boris,
    interp_order=1,
    particle_data=(electrons, ions),
    collisions=(
        Collision(groups=(electrons, ions), channels=('coulomb',), step_interval=coll_interval),
        Collision(groups=(ions, ions), channels=('coulomb',), self_scatter=True, step_interval=coll_interval),
        Collision(groups=(electrons, electrons), channels=('coulomb',), self_scatter=True, step_interval=coll_interval),
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
    collisions_enabled=True,
    push_enabled=True,
    jdep_enabled=True,
    em_enabled=True,
    velocity_backstep_enabled=True,
)

# # ===========================
# # ===== Compile and Run =====
# # ===========================
# print(f'Setting up "{sim_name}"')
# create_data_dir(data_path)
# create_particles(sim_params, [electrons, ions], data_path)
# update_header(sim_params, project_path=project_path, data_path=data_path)
#
# compile_project(build_path, output=True)
# run_project(build_path + '/game_engine', output=True)

# # ===========================
# # ===== Post Processing =====
# # ===========================
xs = np.linspace(xmin, xmax, shape[0])
zs = np.linspace(zmin, zmax, shape[2])
# # Jxs = []
# Jys = []
# # Jzs = []
# time = 0
# with FileReader(data_path + f'/fields_{1500:010d}.bp') as f:
#     # Jxs.append(f.read("Ex")[:, 0, :].T)
#     Jys.append(f.read("Ey")[:, 0, :].T)
#     # Jzs.append(f.read("Ez")[:, 0, :].T)
#     time = f.read("Time")
#
# fig, ax = plt.subplots(3, 1, figsize=(6, 12), layout='constrained')
#
# fig.suptitle(f't = {time * s_to_ns} ns')
#
# # a0 = ax[0].contourf(Jxs[-1], levels=70, cmap='jet', vmin=-3e15, vmax=2e15)
# # ax[0].set_xlabel('x')
# # ax[0].set_ylabel('z')
# # ax[0].set_title('Ex')
#
# a1 = ax[1].contourf(xs, zs, Jys[-1], levels=70, cmap='jet')
# ax[1].set_xlabel('x')
# ax[1].set_ylabel('z')
# ax[1].set_title('Ey')
#
# # a2 = ax[2].contourf(Jzs[-1], levels=70, cmap='jet', vmin=-1e15, vmax=1e15)
# # ax[2].set_xlabel('x')
# # ax[2].set_ylabel('z')
# # ax[2].set_title('Ez')
#
# # plt.colorbar(a0, ax=ax[0])
# plt.colorbar(a1, ax=ax[1])
# # plt.colorbar(a2, ax=ax[2])
#
# plt.show()

from matplotlib.cm import ScalarMappable
from matplotlib import colors

filename = f'/electrons_{7500:010d}.bp'
with FileReader(data_path + filename) as f:
    density = f.read("Density")[:, 0, :].T
    temp = f.read('Temperature')[:, 0, :].T
    time = f.read('Time')

den_norm = colors.LogNorm(vmin=1e27, vmax=2e29)
temp_norm = colors.LogNorm(vmin=Tcold, vmax=4.2e4)

fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')

# im = ax.contourf(xs[:-1], zs[:-1], temp, levels=np.logspace(1, 5, 50), cmap='jet', norm=temp_norm)
# plt.colorbar(ScalarMappable(norm=temp_norm, cmap='jet'), ax=ax)
ax.set_xlabel('x')
ax.set_ylabel('z')
ax.set_title(f't = {time * s_to_ns}')

# im = ax.contourf(xs[:-1], zs[:-1], density, levels=50, cmap='jet')

im = ax.contourf(density, levels=np.logspace(27, 29, 50), cmap='jet', norm=den_norm)

# plt.colorbar(ScalarMappable(norm=den_norm, cmap='jet'), ax=ax)
plt.colorbar(im, ax=ax)

plt.show()

# def calculate_Temp(nn, m, group_name):
#     filename = f'/{group_name}_dump_{nn:010d}.bp'
#     with FileReader(data_path + filename) as f:
#         weight = f.read('Weight')
#         velocity = f.read('Velocity')
#
#     ttl_weight = weight.sum()
#     avg_velocity = (weight * velocity).sum() / ttl_weight # is multiply by weight, summing, and then dividing doing anything?
#     dv = velocity - avg_velocity
#     ttl_sum_dv2 = (weight * (dv**2).sum(axis=1)[:, None]).sum(axis=0)[0]
#     avg_temp = ttl_sum_dv2 * m / (3.0 * constants.e * ttl_weight)
#     return avg_temp
#
# step = save_interval
# stop = nt + step
#
# e_data = []
# for n in range(0, stop, step):
#     e_data.append(calculate_Temp(n, constants.m_e, 'electrons'))
# e_data = np.array(e_data)
#
# i_data = []
# for n in range(0, stop, step):
#     i_data.append(calculate_Temp(n, constants.m_p, 'ions'))
# i_data = np.array(i_data)
#
# with FileReader(data_path + '/particles_energy.bp') as f:
#     variables = f.available_variables()
#     steps = int(variables['Time']['AvailableStepsCount'])
#     time = f.read('Time', step_selection=[0, steps])
#     e_energy = f.read('electrons', step_selection=[0, steps])
#     i_energy = f.read('ions', step_selection=[0, steps])
#
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
# time *= s_to_fs
# field_energy = (ex + ey + ez + bx + by + bz) * J_to_kJ / dy
# e_energy = e_energy * J_to_kJ / dy
# i_energy = i_energy * J_to_kJ / dy
#
# fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8, 12), layout='constrained')
# ax00 = ax0.twinx()
#
# p1, = ax0.plot(time, e_data,  linestyle='-', label='e Temp')
# p2, = ax0.plot(time, i_data,  linestyle='-', label='i Temp')
# # ax0.set_ylim([0, 11])
# ax0.set_xlabel("Step")
# ax0.set_ylabel('Temp (eV)')
#
# p3, = ax00.plot(time, e_energy, 'b-', label='e Energy')
# p4, = ax00.plot(time, i_energy, 'r-', label='i Energy')
# ax00.set_xlabel('Time (s)')
# ax00.set_ylabel('Energy (eV)')
#
# ax0.legend(handles=[p1, p2, p3, p4])
#
# ax1.plot(time, field_energy)
#
# plt.show()
