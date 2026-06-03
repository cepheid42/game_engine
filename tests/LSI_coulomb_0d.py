#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from adios2 import FileReader

from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = f'lsi_coulomb_0d'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

shape = (2, 2, 2)

xmin, xmax = 0.0, 2e-8
ymin, ymax = 0.0, 0.01
zmin, zmax = 0.0, 2e-8

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
px_range = (xmin, xmax) # meters
py_range = (ymin, ymax)
pz_range = (zmin, zmax)

ppc = (10, 1, 10)
density = 1.0e29 #8.5e27
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

coulombParam = CoulombParams(10, 1) # set LnLambda = 0 to calculate it on the fly

particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs=ParticleBCType.Periodic,
    push_type=ParticlePushType.Boris,
    interp_order=1,
    particle_data=(electrons, ions),
    collisions=(
        Collision(groups=(electrons, ions), channels=('coulomb',), step_interval=coll_interval),
        Collision(groups=(ions, ions), channels=('coulomb',), self_scatter=True, step_interval=coll_interval),
        Collision(groups=(electrons, electrons), channels=('coulomb',), self_scatter=True, step_interval=coll_interval),
    )
)

# # ==================================
# # ===== Electromagnetic Params =====
# # ==================================
# em_params = EMParams(
#     save_interval=save_interval,
#     pml_depth=15,
#     em_bcs=(1, 1, 2, 2, 1, 1),
#     laser_enabled=True,
# )

# ==========================
# ===== Metrics Params =====
# ==========================
metric_params = Metrics(
    data_path,
    (
        MetricType.ParticleEnergy,
        # MetricType.FieldEnergy,
        # MetricType.FieldDump,
        MetricType.ParticleDump,
        # MetricType.ParticleDiagnostics,
    )
)

# ============================
# ===== Simulation Class =====
# ============================
sim_params = Simulation(
    name=sim_name,
    shape=shape,
    nthreads=8,
    dt=dt,
    t_end=t_end,
    nt=nt,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    particle_params=particle_params,
    metric_params=metric_params,
    collisions_enabled=True,
    push_enabled=False,
    jdep_enabled=False,
    em_enabled=False,
    velocity_backstep_enabled=False,
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
def calculate_Temp(nn, m, group_name):
    filename = f'/{group_name}_dump_{nn:010d}.bp'
    with FileReader(data_path + filename) as f:
        weight = f.read('Weight')
        velocity = f.read('Velocity')

    ttl_weight = weight.sum()
    avg_velocity = (weight * velocity).sum() / ttl_weight # is multiply by weight, summing, and then dividing doing anything?
    dv = velocity - avg_velocity
    ttl_sum_dv2 = (weight * (dv**2).sum(axis=1)[:, None]).sum(axis=0)[0]
    avg_temp = ttl_sum_dv2 * m / (3.0 * constants.e * ttl_weight)
    return avg_temp

step = save_interval
stop = nt + step

e_data = []
for n in range(0, stop, step):
    e_data.append(calculate_Temp(n, constants.m_e, 'electrons'))
e_data = np.array(e_data)

i_data = []
for n in range(0, stop, step):
    i_data.append(calculate_Temp(n, constants.m_p, 'ions'))
i_data = np.array(i_data)

with FileReader(data_path + '/particles_energy.bp') as f:
    variables = f.available_variables()
    steps = int(variables['Time']['AvailableStepsCount'])
    time = f.read('Time', step_selection=[0, steps])
    e_energy = f.read('electrons', step_selection=[0, steps])
    i_energy = f.read('ions', step_selection=[0, steps])

fig, ax0 = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
ax1 = ax0.twinx()

p1, = ax0.plot(time, e_data,  linestyle='-', label='e Temp')
p2, = ax0.plot(time, i_data,  linestyle='-', label='i-Temp')
ax0.set_ylim([0, 11])
ax0.set_xlabel("Step")
ax0.set_ylabel('Temp (eV)')

p3, = ax1.plot(time, e_energy, 'b-', label='e Energy')
p4, = ax1.plot(time, i_energy, 'r-', label='i Energy')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Energy (eV)')

ax0.legend(handles=[p1, p2, p3, p4])

plt.show()

