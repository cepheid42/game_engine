#!/usr/bin/env python3

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from adios2 import FileReader, Stream
from scipy import constants

from scripts.particle_generation import create_particles
from scripts.domain_params import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'uniform_B_field'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

'''
Tests from https://iopscience.iop.org/article/10.3847/1538-4365/aab114
'''

shape = (16, 16, 16)

Bz_amp = constants.c #* 1e6
gamma_an = 1e6
Tc = 2 * np.pi * gamma_an / Bz_amp

xmin, xmax = -10.1, 10.1 # meters
ymin, ymax = -10.1, 10.1
zmin, zmax = -10.1, 10.1

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = Tc / 100 # seconds
t_end = 10000 * dt # seconds
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

save_interval = nt // 100

# =====================
# ===== Particles =====
# =====================
px_range = (0.0, 0.0) # meters
py_range = (0.0, 0.0)
pz_range = (0.0, 0.0)

mass = 1.0
charge = 1.0 / constants.e
v_perp = 0.00005 * constants.c * (1.0 - 5e-13)

singleton = Particles(
    name='singleton',
    mass=mass,
    charge=charge,
    atomic_number=0,
    # actually velocity for sp_uniformB distribution
    temp=(-v_perp, 0.0, 0.0),
    density=1.0, # m^-3,
    ppc=(1, 1, 1),
    distribution='sp_uniformB',
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
    bc_depth=0,
    interp_order=1,
    particle_data=(singleton,)
)

# ==================================
# ===== Electromagnetic Params =====
# ==================================
Bz_applied = np.full((shape[0] - 1, shape[1] - 1, shape[2]), Bz_amp)
with Stream(data_path + f'/{sim_name}_applied_fields.bp', 'w') as f:
    f.write('Bz', Bz_applied, Bz_applied.shape, (0, 0, 0), Bz_applied.shape)

em_params = EMParams(
    save_interval=save_interval,
    applied_fields=data_path + f'/{sim_name}_applied_fields.bp'
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
    name=sim_name,
    shape=shape,
    nthreads=1,
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
    em_enabled=False,
    jdep_enabled=False,
    collisions_enabled=False
)

# ===========================
# ===== Compile and Run =====
# ===========================
print(f'Setting up "{sim_name}"')

if not os.path.exists(data_path):
    print(f'Creating simulation data directory "{data_path}"...')
    os.makedirs(data_path)

create_particles(sim_params, singleton, data_path)
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
pos = []
for n in range(0, 2500, 100): #nt, save_interval):
    file = f'/singleton_dump_{n:010d}.bp'
    with FileReader(data_path + file) as f:
        pos.append(f.read("Position"))

pos = np.array(pos).reshape(-1, 3)

print(pos)
# Rc_an = gamma_an * 1.0e3 * v_perp / (1.0 * Bz_amp)
# p_Rc = plt.Circle((0, 0), Rc_an, fill=False)

fig, ax = plt.subplots(figsize=(6, 6), layout='constrained')
ax.set_aspect('equal')
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
# ax.add_patch(p_Rc)

ax.plot(pos[:, 0], pos[:, 1])

plt.show()