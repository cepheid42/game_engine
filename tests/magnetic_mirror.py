#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from adios2 import FileReader, Stream

from scripts.particle_generation import create_particles
from scripts.domain_params import *
from scripts.utilities import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'magnetic_mirror'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'

'''
Tests from https://iopscience.iop.org/article/10.3847/1538-4365/aab114
'''

L = 1.0e7
B0 = 1.0e6
v_perp = velocity_from_gamma(100.0000000001) / np.sqrt(2)
v = np.sqrt(v_perp * v_perp + v_perp * v_perp)
gamma = 1.0 / np.sqrt(1.0 - (v / constants.c)**2)
Rc = gamma * v_perp / B0

shape = (512, 512, 128)

xmin, xmax = -20.0 * L, 20 * L
ymin, ymax = -20.0 * L, 20 * L
zmin, zmax = -2.0e6 * L, 2.0e6 * L

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 1.0e-7 # seconds
t_end = np.pi / 10# seconds
nt = int(t_end / dt) + 1

save_interval = 500000

# =====================
# ===== Particles =====
# =====================
mass = 1.0
charge = 1.0 / constants.e

px_range = (-Rc, dx) # meters
py_range = (0.0, dy)
pz_range = (0.0, dz)

single_particle = Particles(
    name='sp',
    mass=mass,
    charge=charge,
    atomic_number=0,
    tracer=True,
    temp=(0.0, v_perp, v_perp), # eV
    density=1.0, # m^-3,
    ppc=(1, 1, 1),
    distribution='sp_magneticmirror',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
particle_params = ParticleParams(
    save_interval=save_interval,
    bc_depth=0,
    interp_order=1,
    particle_data=(single_particle,)
)

# ============================
# ===== Simulation Class =====
# ============================
pushers = [
    ParticlePushType.Boris,
    ParticlePushType.HC
]
sim_names = [sim_name + '_' + pusher.get_name() for pusher in pushers]

# Bx, By, Bz = gen_magnetic_mirror_fields(B0, L, np.linspace(zmin, zmax, shape[2]), shape)

xf, yf, zf = np.meshgrid(
    np.linspace(xmin, xmax, shape[0]),
    np.linspace(ymin, ymax, shape[1]),
    np.linspace(zmin, zmax, shape[2]),
    indexing='ij'
)

xh, yh, zh = np.meshgrid(
    np.linspace(xmin + dx/2, xmax + dx/2, shape[0]),
    np.linspace(ymin + dy/2, ymax + dy/2, shape[1]),
    np.linspace(zmin + dz/2, zmax + dz/2, shape[2]),
    indexing='ij'
)

Bx = -xf[:, :-1, :-1] * B0 * zh[:, :-1, :-1] / L**2
By = -yf[:-1, :, :-1] * B0 * zh[:-1, :, :-1] / L**2
Bz = B0 * zf[:-1, :-1, :]# (1.0 + (zf[:-1, :-1, :] / L)**2)

for pusher, name in zip(pushers, sim_names):
    data_path = project_path + f'/data/{name}'
    fields_path = data_path + f'/{name}_applied_fields.bp'

    with Stream(fields_path, 'w') as f:
        f.write('Bx', Bx, Bx.shape, (0, 0, 0), Bx.shape)
        f.write('By', By, By.shape, (0, 0, 0), By.shape)
        f.write('Bz', Bz, Bz.shape, (0, 0, 0), Bz.shape)

    em_params = EMParams(save_interval=save_interval, applied_fields=fields_path)
    metric_params = Metrics(data_path, (MetricType.ParticleDump,))
    particle_params.push_type = pusher

    sim_params = Simulation(
        name=name,
        shape=shape,
        nthreads=1,
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
        em_enabled=False,
        jdep_enabled=False,
        collisions_enabled=False,
        applied_fields_only=True
    )

    # ===========================
    # ===== Compile and Run =====
    # ===========================
    print(f'Setting up "{name}"')
    create_data_dir(data_path)
    create_particles(sim_params, single_particle, data_path)
    update_header(sim_params, project_path=project_path)

    compile_project(build_path, output=True)
    run_project(build_path + '/game_engine', output=True)

# ===========================
# ===== Post Processing =====
# ===========================
skip = 100

sims = dict()
for name in sim_names:
    sim_path = project_path + f'/data/{name}'
    # vel = []
    pos = []
    gs = []
    times = []
    for n in range(save_interval, nt, save_interval):
        file = f'/{single_particle.name}_tracer_{n:010d}.bp'
        with FileReader(sim_path + file) as f:
            # vel.append(f.read('Velocity')[::skip])
            pos.append(f.read('Position')[::skip])
            gs.append(f.read('Gamma')[::skip])
            times.append(f.read("Time")[::skip])

    sims[name] = ParticlePlotData(
        # velocities=np.array(vel).reshape(-1, 3),
        positions=np.array(pos).reshape(-1, 3),
        gammas=np.array(gs).flatten(),
        times=np.array(times).flatten()
    )


plot_params = [
    ('--', 'r', 'P', 8, 'full'), # Boris
    ('--', 'b', 'D', 8, 'none')  # HC
]

fig, ax = plt.subplots(1, 3, figsize=(14, 6), layout='constrained')
ax[0].set_xlabel('t')
ax[0].set_ylabel('x')
ax[1].set_xlabel('t')
ax[1].set_ylabel('y')
ax[2].set_xlabel('t')
ax[2].set_ylabel('z')
# # ax.set_xlim([0, np.pi])
# # ax.set_ylim([-1.5e7, 1.5e7])
#
for i, (name, data) in enumerate(sims.items()):
    name = name.split('_')[-1]
    ls, c, m, ms, fs = plot_params[i]

    mark_every = data.times.shape[0] // 20
    ax[0].plot(data.times, data.positions[:, 0], ls=ls, c=c, label=name)
    ax[1].plot(data.times, data.positions[:, 1], ls=ls, c=c, label=name)
    ax[2].plot(data.times, data.positions[:, 2], ls=ls, c=c, label=name)

ax[0].legend()
ax[1].legend()
ax[2].legend()
plt.show()
