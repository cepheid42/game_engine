#!/usr/bin/env python3

# import sys
# sys.path.append('/home/cepheid/TriForce/game_engine')

import numpy as np
from scipy import constants

from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'lsi_smith_reduced_10eV'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'


shape = (1551, 2, 451)

xmin, xmax = -15.5e-6, 15.5e-6
ymin, ymax = 0.0, 0.01
zmin, zmax = -15.5e-6, 15.5e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 4.0e-17
t_end = 3.0e-13
nt = int(t_end / dt) + 1

# =====================
# ===== Particles =====
# =====================
px_range = (-5e-7, 5e-7) # meters
py_range = (ymin, ymax)
pz_range = (-1e-5, 1e-5)

ppc = (10, 1, 10)
density = 8.5e27 #m^-3
temp = tuple(3 * [10 / np.sqrt(3)])

protons = Particles(
    name='protons',
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

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
coll_interval = 1
production_mult = 1.0e8
particle_params = ParticleParams(
    particle_bcs=ParticleBCType.Outflow,
    push_type=ParticlePushType.Boris,
    interp_order=1,
    particle_data=(protons, electrons),
    # collisions=(
        # Collision(
        #     groups=(electrons, electrons),
        #     channels=('coulomb',),
        #     coulomb=CoulombParams(0.0, 1.0), # set LnLambda = 0 to calculate it on the fly
        #     self_scatter=True,
        #     step_interval=coll_interval
        # ),
        # Collision(
        #     groups=(electrons, protons),
        #     channels=('radiation',),
        #     coulomb=CoulombParams(0.0, 1.0), # set LnLambda = 0 to calculate it on the fly
        #     self_scatter=False,
        #     step_interval=coll_interval
        # ),
        # Collision(
        #     groups=(protons, protons),
        #     channels=('coulomb',),
        #     coulomb=CoulombParams(0.0, 1.0), # set LnLambda = 0 to calculate it on the fly
        #     self_scatter=True,
        #     step_interval=coll_interval
        # ),
    # )
)


# ==================================
# ===== Electromagnetic Params =====
# ==================================
em_params = EMParams(
    pml_depth=15,
    em_bcs=(1, 1, 2, 2, 1, 1),
    laser_spec=Laser(8.0e-7, -2.75e13, 2.5479e-6, 15.0e-6, 1.288)
)

# ==========================
# ===== Metrics Params =====
# ==========================
energy_interval = 1
dump_interval = 750
metric_params = Metrics(
    data_path,
    (
        (MetricType.ParticleEnergy, energy_interval),
        (MetricType.FieldEnergy, energy_interval),
        (MetricType.ParticleDump, dump_interval),
        (MetricType.ParticleDiagnostics, dump_interval),
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
    # collisions_enabled=True,
    push_enabled=True,
    jdep_enabled=True,
    velocity_backstep_enabled=True,
    collisions_enabled=False,
    # push_enabled=False,
    # jdep_enabled=False,
    # velocity_backstep_enabled=False,
    em_enabled=True
)

# ===========================
# ===== Compile and Run =====
# ===========================
run = True
# run = False

if run:
    print(f'Setting up "{sim_name}"')
    create_data_dir(data_path)
    create_particles(sim_params, (protons, electrons), data_path)
    update_header(sim_params, project_path=project_path, data_path=data_path)

    compile_project(build_path, output=True)
    run_project(build_path + '/game_engine', output=True)

# ===========================
# ===== Post Processing =====
# ===========================
data_path = project_path + f'/data/lsi_smith'
smith_data = '/home/cepheid/TriForce/game_engine/tests/smith_data'

# # Figure 2
# plot_energies(data_path, smith_data)
#
# # Figure 3
# plot_spectra(data_path, smith_data, 7500)
#
# # Figure 4
# # plot_density(['protons', 'electrons'], 7500, data_path, xs, zs)
#
# # # Figure 5
# # plot_average_density(data_path, smith_data, 7500)
#
# lsp_energies = np.genfromtxt(smith_data + '/Energy_Diagnostics/LSP_E/history.txt', comments='#', usecols=(1, 3, 5, 6))
# lsp_energies[:, 0] *= 1.0e6 # time ns -> fs
#
# # Read EM Field Energies
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
# field_energy = (ex + ey + ez + bx + by + bz)
#
# # Read Particle Energies
# with FileReader(data_path + '/particles_energy.bp') as f:
#     variables = f.available_variables()
#     steps = int(variables['Time']['AvailableStepsCount'])
#     time = f.read('Time', step_selection=[0, steps])
#     e_energy = f.read('electrons', step_selection=[0, steps])
#     i_energy = f.read('protons', step_selection=[0, steps])
#
# time *= s_to_fs
#
# y_labels = [
#     (r'Field (J cm$^{-1}$)', [0, 100]),
#     (r'Electron (J cm$^{-1}$)', [0, 15]),
#     (r'Proton (J cm$^{-1}$)', [0, 15])
# ]
#
# fig, ax = plt.subplots(3, 1, figsize=(6, 10), layout='constrained')
#
# # Plot EM Energies
# ax[0].plot(time, field_energy, 'b-', label='TriForce Reduced')
# ax[0].plot(lsp_energies[:, 0], lsp_energies[:, 2], 'r--', label='LSP(E)')
#
# # Plot Electron Energies
# ax[1].plot(time, e_energy, 'b-', label='TriForce Reduced')
# ax[1].plot(lsp_energies[:, 0], lsp_energies[:, 1], 'r--', label='LSP(E)')
#
# # Plot Ion Energies
# ax[2].plot(time, i_energy, 'b-', label='TriForce Reduced')
# ax[2].plot(lsp_energies[:, 0], lsp_energies[:, 3], 'r--', label='LSP(E)')
#
# # Read EM Field Energies
# data_path = project_path + f'/data/lsi_smith_reduced'
# with FileReader(data_path + '/fields_energy.bp') as f:
#     variables = f.available_variables()
#     steps = int(variables['Time']['AvailableStepsCount'])
#     ex = f.read('Ex Energy', step_selection=[0, steps])
#     ey = f.read('Ey Energy', step_selection=[0, steps])
#     ez = f.read('Ez Energy', step_selection=[0, steps])
#     bx = f.read('Bx Energy', step_selection=[0, steps])
#     by = f.read('By Energy', step_selection=[0, steps])
#     bz = f.read('Bz Energy', step_selection=[0, steps])
# field_energy = (ex + ey + ez + bx + by + bz)
# ax[0].plot(time, field_energy, 'g-.', label='TriForce Full')
#
# # Read Particle Energies
# with FileReader(data_path + '/particles_energy.bp') as f:
#     variables = f.available_variables()
#     steps = int(variables['Time']['AvailableStepsCount'])
#     # time = f.read('Time', step_selection=[0, steps])
#     e_energy = f.read('electrons', step_selection=[0, steps])
#     i_energy = f.read('protons', step_selection=[0, steps])
#
# ax[1].plot(time, e_energy, 'g-.', label='TriForce Full')
#
# # Plot Ion Energies
# ax[2].plot(time, i_energy, 'g-.', label='TriForce Full')
#
# for i, a in enumerate(ax):
#     label, ylim = y_labels[i]
#     a.grid()
#     a.set_xlabel('Time (fs)')
#     a.set_ylabel(label)
#     # a.set_ylim(ylim)
#     a.legend()
#
# plt.show()
