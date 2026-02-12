#!/usr/bin/env python3

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from adios2 import FileReader
from scipy import constants

from scripts.particle_generation import create_particles
from scripts.domain_params import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'carbon_thermal_eq'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

shape = (2, 2, 2)
xmin, xmax = 0.0, 1e-6
ymin, ymax = 0.0, 1e-6
zmin, zmax = 0.0, 1e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

t_end = 5e-12
dt = 5e-16
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = nt // 100

# =====================
# ===== Particles =====
# =====================
px_range = (xmin, xmax)
py_range = (ymin, ymax)
pz_range = (zmin, zmax)
ppc = (20, 20, 20)
Z_carbon = 6
m_carbon = 1.9945e-26 # kg
n_carbon = 1e26 # m^3
T_hot = 1250 / np.sqrt(3)
T_cold = 250 / np.sqrt(3)

carbon1 = Particles(
    name='carbon1',
    mass=m_carbon,
    atomic_number=Z_carbon,
    charge=Z_carbon,
    temp=(T_hot, T_hot, T_hot), # eV
    density=n_carbon, # m^-3,
    ppc=ppc,
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

carbon2 = Particles(
    name='carbon2',
    mass=m_carbon,
    atomic_number=Z_carbon,
    charge=Z_carbon,
    temp=(T_cold, T_cold, T_cold), # eV
    density=n_carbon, # m^-3,
    ppc=ppc,
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
particle_params = ParticleParams(
    save_interval=save_interval,
    particle_bcs='periodic',
    interp_order=1,
    particle_data=(carbon1, carbon2),
    collisions=(
        Collision(groups=('carbon1', 'carbon2'), channels=('coulomb',)),
        Collision(groups=('carbon1', 'carbon1'), channels=('coulomb',), self_scatter=True),
        Collision(groups=('carbon2', 'carbon2'), channels=('coulomb',), self_scatter=True),
    )
)

# ==========================
# ===== Metrics Params =====
# ==========================
metric_params = Metrics(
    data_path,
    (MetricType.ParticleDump,)
)

# ============================
# ===== Simulation Class =====
# ============================
sim_params = Simulation(
    name='carbon_thermal_eq',
    shape=shape,
    nthreads=4,
    dt=dt,
    t_end=t_end,
    nt=nt,
    cfl=cfl,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    particle_params=particle_params,
    metric_params=metric_params,
    push_enabled=False,
    jdep_enabled=False,
    em_enabled=False
)

# ===========================
# ===== Compile and Run =====
# ===========================
print(f'Setting up "{sim_name}"')
create_particles(sim_params, carbon1, data_path)
create_particles(sim_params, carbon2, data_path)
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
def calculate_Temp(t, group_name):
    filename = f'/{group_name}_dump_{t:010d}.bp'
    with FileReader(data_path + filename) as f:
        weight = f.read('Weight')
        velocity = f.read('Velocity')

    ttl_weight = weight.sum()
    avg_velocity = (weight * velocity).sum() / ttl_weight # is multiply by weight, summing, and then dividing doing anything?
    dv = velocity - avg_velocity
    ttl_sum_dv2 = (weight * (dv**2).sum(axis=1)[:, None]).sum(axis=0)[0]
    avg_temp = ttl_sum_dv2 * m_carbon / (3.0 * constants.e * ttl_weight)
    return avg_temp

step = save_interval
stop = nt + step

carbon1_data = []
for n in range(0, stop, step):
    carbon1_data.append(calculate_Temp(n, 'carbon1'))
carbon1_data = np.array(carbon1_data)

carbon2_data = []
for n in range(0, stop, step):
    carbon2_data.append(calculate_Temp(n, 'carbon2'))
carbon2_data = np.array(carbon2_data)

time = np.linspace(0, t_end, stop // step) # picoseconds

plt.style.use('/home/cepheid/TriForce/game_engine/data/triforce.mplstyle')
good_colors = ['#db6d00', '#006ddb', '#920000', '#52a736', '#9B30FF']

fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')

eV_to_erg = 1.60218e-12
s_to_ps = 1.0e12
m_carbon_g = 1.9945e-23 # g
n_carbon_cm = 1e20 # cm^-3
coulomb_log = 10.0

coef = (8 / 3) * np.sqrt(np.pi)
Ts = 0.5 * (T_hot + T_cold) * eV_to_erg
mu = coef * (Z_carbon * 4.8e-10)**4 * coulomb_log * n_carbon_cm / (m_carbon_g**0.5 * Ts**1.5)

Tc1 = 0.5 * (T_hot + T_cold) + 0.5 * (T_hot - T_cold) * np.exp(-mu * time)
Tc2 = 0.5 * (T_hot + T_cold) + 0.5 * (T_cold - T_hot) * np.exp(-mu * time)

ax.plot(time * s_to_ps, Tc1, '-k', label='Fluid Model')
ax.plot(time * s_to_ps, Tc2, '-k', label='Fluid Model')

ax.plot(time * s_to_ps, carbon1_data,  linestyle='-', c=good_colors[0], label='Carbon1')
ax.plot(time * s_to_ps, carbon2_data,  linestyle='-', c=good_colors[1], label='Carbon2')

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Temperature (eV)')
# ax.set_xlim([0, 5])
# ax.set_ylim([200, 1300])
ax.legend()

# plt.savefig(data_path + f'/carbon_eq_test_temperature.png')
# plt.close(fig)
plt.show()
