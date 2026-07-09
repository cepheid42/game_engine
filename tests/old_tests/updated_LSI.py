#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from adios2 import FileReader
from datetime import datetime

from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
timestamp = datetime.today().strftime('%y%m%d%H%M')
sim_name = f'lsi_test_coulomb_{timestamp}'
project_path = '/'
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

save_interval = 75

# =====================
# ===== Particles =====
# =====================
px_range = (-5e-7, 5e-7) # meters
py_range = (ymin, ymax)
pz_range = (-1e-5, 1e-5)

ppc = (10, 1, 10)
density = 8.5e27
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
    particle_bcs=ParticleBCType.Outflow,
    push_type=ParticlePushType.Boris,
    interp_order=1,
    particle_data=(electrons, ions),
    collisions=(
        Collision(groups=(electrons, ions), channels=('coulomb',)),
        Collision(groups=(ions, ions), channels=('coulomb',), self_scatter=True),
        Collision(groups=(electrons, electrons), channels=('coulomb',), self_scatter=True),
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
        MetricType.ParticleDump,
        # MetricType.ParticleDiagnostics
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
# xs = np.linspace(xmin, xmax, shape[0])
# zs = np.linspace(zmin, zmax, shape[2])
#
# # for n in range(0, nt, save_interval):
# #     with FileReader(data_path + f'/fields_{n:010d}.bp') as f:
# #         print(n, f.read('Time'))
#
# # with FileReader(data_path + f'/ions_{3000:010d}.bp') as f:
# #     density = f.read("Density")[:, 0, :].T
# #
# # from matplotlib.cm import ScalarMappable
# # from matplotlib import colors, ticker
# # fig, ax = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
# #
# # norm = colors.LogNorm(vmin=1e26, vmax=1e28)
# # im = ax.contourf(density, levels=np.logspace(26, 28, 50), norm=norm,  cmap='jet')
# # fig.colorbar(ScalarMappable(norm=norm, cmap='jet'), ax=ax, shrink=0.82)
# # plt.show()
#
# # plt.savefig(data_path + f'/lsi_density_80fs.png')
# # plt.close(fig)
#
# # Jxs = []
# # Jys = []
# # Jzs = []
# # time = 0
# # with FileReader(data_path + f'/fields_{2500:010d}.bp') as f:
# #     # Jxs.append(f.read("Jx")[:, 0, :].T)
# #     # Jys.append(f.read("Jy")[:, 0, :].T)
# #     # Jzs.append(f.read("Jz")[:, 0, :].T)
# #
# #     # Jxs.append(f.read("Ex")[:, 0, :].T)
# #     # Jys.append(f.read("Ey")[:, 0, :].T)
# #     # Jzs.append(f.read("Ez")[:, 0, :].T)
# #
# #     Jxs.append(f.read("Hx")[:, 0, :].T * constants.mu_0)
# #     Jys.append(f.read("Hy")[:, 0, :].T * constants.mu_0)
# #     Jzs.append(f.read("Hz")[:, 0, :].T * constants.mu_0)
# #
# #     time = f.read("Time")
# #
# # fig, ax = plt.subplots(3, 1, figsize=(6, 12), layout='constrained')
# #
# # fig.suptitle(f't = {time * s_to_ns} ns')
# #
# # # a0 = ax[0].contourf(xs[:-1], zs, Jxs[-1], levels=70, cmap='jet', vmin=-5e12, vmax=5e12)
# # a0 = ax[0].contourf(xs, zs[:-1], Jxs[-1], levels=70, cmap='jet')
# # ax[0].set_xlabel('x')
# # ax[0].set_ylabel('z')
# # ax[0].set_title('Jx')
# #
# # # a1 = ax[1].contourf(xs, zs, Jys[-1], levels=70, cmap='jet')
# # a1 = ax[1].contourf(xs[:-1], zs[:-1], Jys[-1], levels=70, cmap='jet')
# # ax[1].set_xlabel('x')
# # ax[1].set_ylabel('z')
# # ax[1].set_title('Jy')
# #
# # # a2 = ax[2].contourf(xs, zs[:-1], Jzs[-1], levels=70, cmap='jet')
# # a2 = ax[2].contourf(xs[:-1], zs, Jzs[-1], levels=70, cmap='jet')
# #
# # ax[2].set_xlabel('x')
# # ax[2].set_ylabel('z')
# # ax[2].set_title('Jz')
# #
# # plt.colorbar(a0, ax=ax[0])
# # plt.colorbar(a1, ax=ax[1])
# # plt.colorbar(a2, ax=ax[2])
# #
# # plt.show()

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
# smith_field_data = np.genfromtxt('./data/smith_lsi_field_energy.csv', delimiter=',')
# smith_electron_data = np.genfromtxt('./data/smith_lsi_electron_energy.csv', delimiter=',')
# smith_proton_data = np.genfromtxt('./data/smith_lsi_proton_energy.csv', delimiter=',')
#
# y_labels = [
#     (r'Field (kJ m$^{-1}$)', smith_field_data, [0, 100]),
#     (r'Electron (kJ m$^{-1}$)', smith_electron_data, [0, 15]),
#     (r'Proton (kJ m$^{-1}$)', smith_proton_data, [0, 15])
# ]
#
# fig, ax = plt.subplots(3, 1, figsize=(6, 10), layout='constrained')
#
# for i, a in enumerate(ax):
#     label, _, _ = y_labels[i]
#     a.grid()
#     a.set_xlabel('Time (fs)')
#     a.set_ylabel(label)
#
# ax[0].plot(time, field_energy, 'b-', label='TriForce')
# ax[1].plot(time, e_energy, 'b-')
# ax[2].plot(time, i_energy, 'b-')
#
# ax[0].plot(smith_field_data[:, 0], smith_field_data[:, 1], 'r--', label='Smith')
# ax[1].plot(smith_electron_data[:, 0], smith_electron_data[:, 1], 'r--')
# ax[2].plot(smith_proton_data[:, 0], smith_proton_data[:, 1], 'r--')
#
# # ax[0].set_ylim([0, 100])
# # ax[1].set_ylim([0, 15])
# # ax[2].set_ylim([0, 15])
#
# ax[0].legend()
#
# # plt.savefig(data_path + f'/lsi_comp_normal.png')
# # plt.close(fig)
#
# plt.show()
