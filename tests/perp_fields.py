#!/usr/bin/env python3

import matplotlib.pyplot as plt
from adios2 import FileReader, Stream

from scripts.particle_generation import create_particles
from scripts.domain_params import *
from scripts.utilities import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'perp_fields'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'

'''
Tests from https://iopscience.iop.org/article/10.3847/1538-4365/aab114
'''

shape = (16, 16, 16)

xmin, xmax = -1.0e13, 1.0e13 # meters
ymin, ymax = -2.1e16, 100.0
zmin, zmax = -1.0, 1.0

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 0.5 # seconds
t_end = 2.0 * np.pi * 1e7 # seconds
nt = int(t_end / dt) + 1

save_interval = 500000

# =====================
# ===== Particles =====
# =====================
mass = 1.0
charge = 1.0 / constants.e

Bz_amp = 1.0
Ex_amp = constants.c * (1.0 - 5.0e-5)
v_e = -Ex_amp
kappa = gamma_from_velocity(v_e)

px_range = (0.0, dx) # meters
py_range = (0.0, dy)
pz_range = (0.0, dz)

single_particle = Particles(
    name='sp',
    mass=mass,
    charge=charge,
    atomic_number=0,
    tracer=True,
    temp=(0.0, 0.0, 0.0), # eV
    density=1.0, # m^-3,
    ppc=(1, 1, 1),
    distribution='sp_perpfields',
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
Ex_applied = np.full((shape[0] - 1, shape[1], shape[2]), Ex_amp)
Bz_applied = np.full((shape[0] - 1, shape[1] - 1, shape[2]), Bz_amp)

for pusher, name in zip(pushers, sim_names):
    data_path = project_path + f'/data/{name}'
    fields_path = data_path + f'/{name}_applied_fields.bp'

    with Stream(fields_path, 'w') as f:
        f.write('Ex', Ex_applied, Ex_applied.shape, (0, 0, 0), Ex_applied.shape)
        f.write('Bz', Bz_applied, Bz_applied.shape, (0, 0, 0), Bz_applied.shape)

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
skip = 10000

sims = dict()
for name in sim_names:
    sim_path = project_path + f'/data/{name}'
    vel = []
    pos = []
    gs = []
    times = []
    for n in range(save_interval, nt, save_interval):
        file = f'/{single_particle.name}_tracer_{n:010d}.bp'
        with FileReader(sim_path + file) as f:
            vel.append(f.read('Velocity')[::skip])
            pos.append(f.read('Position')[::skip])
            gs.append(f.read('Gamma')[::skip])
            times.append(f.read("Time")[::skip])

    sims[name] = ParticlePlotData(
        velocities=np.array(vel).reshape(-1, 3),
        positions=np.array(pos).reshape(-1, 3),
        gammas=np.array(gs).flatten(),
        times=np.array(times).flatten()
    )


plot_params = [
    ('--', 'r', 'P', 8, 'full'), # Boris
    ('--', 'b', 'D', 8, 'none')  # HC
]

Bprime = Bz_amp / kappa
Rc_an = kappa * np.abs(v_e) / Bprime # because fuck you

fig, ax = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')

ax[0, 0].set_xlabel('x')
ax[0, 0].set_ylabel('y')
ax[0, 0].set_xlim([0, 8e12])
ax[0, 0].set_ylim([-2e16, 0])

ax[0, 1].set_xlabel('x\'')
ax[0, 1].set_ylabel('y\'')
ax[0, 1].set_aspect('equal')
ax[0, 1].set_xlim([-0.9e12, 6.9e12])
ax[0, 1].set_ylim([-3.1e12, 3.1e12])

ax[1, 0].set_xlim([0, 6.3e7])
ax[1, 0].set_ylim([1e-15, 1])
ax[1, 0].set_yscale('log')
ax[1, 0].set_xlabel('t')
ax[1, 0].set_ylabel(r'$|R_c - R_{c,an}|$ / $R_{c,an}$')

ax[1, 1].set_xlabel('t')
ax[1, 1].set_ylabel(r'$\gamma$')
ax[1, 1].set_ylim([0, 6.3e7])
ax[1, 1].set_ylim([1e1, 1e6])
ax[1, 1].set_yscale('log')

for i, (name, data) in enumerate(sims.items()):
    name = name.split('_')[-1]
    ls, c, m, ms, fs = plot_params[i]
    xp = data.positions[:, 0]
    yp = kappa * (data.positions[:, 1] - v_e * data.times)

    vx = data.velocities[:, 0]
    vy = data.velocities[:, 1]
    coef = 1.0 / (1.0 - (v_e * vy / constants.c**2))

    vxp = (vx / kappa) * coef
    vyp = ((vy / kappa) - v_e + (kappa * vy * v_e**2) / (constants.c**2 * (kappa + 1))) * coef

    gp = kappa * data.gammas * (1.0 - v_e * vy / constants.c**2)
    Rc = gp * np.sqrt(vxp**2 + vyp**2) / Bprime
    Rc_err = np.abs(Rc - Rc_an) / Rc_an

    mark_every = data.times.shape[0] // 20
    ax[0, 0].plot(data.positions[:, 0], data.positions[:, 1], c=c, label=name)
    ax[0, 1].plot(xp, yp, c=c,  label=name)
    ax[1, 0].plot(data.times, Rc_err, c=c, marker=m, ms=ms, markevery=mark_every, label=name)
    ax[1, 1].plot(data.times, data.gammas, c=c, label=name)
    ax[1, 1].plot(data.times, gp, ls=ls, c=c, label=name)

# ax[0, 0].legend()
# ax[0, 1].legend()
# ax[1, 0].legend()
# ax[1, 1].legend()
plt.show()
