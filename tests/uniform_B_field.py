#!/usr/bin/env python3

import matplotlib.pyplot as plt
from adios2 import FileReader, Stream

from scripts.particle_generation import create_particles
from scripts.domain_params import *
from scripts.utilities import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'uniform_B_field'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'

'''
Tests from https://iopscience.iop.org/article/10.3847/1538-4365/aab114
'''

shape = (16, 16, 16)

mass = 1.0
charge = 1.0 / constants.e # charge * e == 1
gamma_an = 1e6
v_perp = velocity_from_gamma(gamma_an)
Bz_amp = gamma_an * v_perp
omega_c = Bz_amp / gamma_an
Tc = 2 * np.pi / omega_c

xmin, xmax = -1.1, 1.1 # meters
ymin, ymax = -1.1, 1.1
zmin, zmax = -1.1, 1.1

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = Tc / 100 # seconds
t_end = 10000 * dt # seconds
nt = int(t_end / dt) + 1

save_interval = 2000

# =====================
# ===== Particles =====
# =====================
px_range = (-1.0, 0.0) # meters
py_range = (0.0, 0.0)
pz_range = (0.0, 0.0)

single_particle = Particles(
    name='sp',
    mass=mass,
    charge=charge,
    atomic_number=0,
    tracer=True,
    # actually velocity for sp_uniformB distribution
    temp=(0.0, -v_perp, 0.0),
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
Bz_applied = np.full((shape[0] - 1, shape[1] - 1, shape[2]), -Bz_amp)

for pusher, name in zip(pushers, sim_names):
    data_path = project_path + f'/data/{name}'
    fields_path = data_path + f'/{name}_applied_fields.bp'

    with Stream(fields_path, 'w') as f:
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
        collisions_enabled=False
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
            vel.append(f.read('Velocity'))
            pos.append(f.read('Position'))
            gs.append(f.read('Gamma'))
            times.append(f.read("Time"))

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

Rc_an = gamma_an * v_perp / Bz_amp
p_Rc = plt.Circle((0, 0), Rc_an, fill=False, label=r'$R_c$')

fig, ax = plt.subplots(1, 3, figsize=(16, 4), layout='constrained')

ax[0].set_xlabel('x')
ax[0].set_ylabel('y')
ax[0].set_aspect('equal')
ax[0].set_xlim([xmin, xmax])
ax[0].set_ylim([ymin, ymax])
ax[0].set_xlim([-1.1, 1.1])
ax[0].set_ylim([-1.1, 1.1])
ax[0].add_patch(p_Rc)

ax[1].set_xlabel(r't / $T_c$')
ax[1].set_ylabel(r'$|\gamma - \gamma_{an}|$ / $\gamma_{an}$')
ax[1].set_xlim([0, 100])
ax[1].set_ylim([1e-16, 1e-13])
ax[1].set_yscale('log')

ax[2].set_xlabel(r't / $T_c$')
ax[2].set_ylabel(r'$|\theta_c - \theta_{an}|$')
ax[2].set_xlim([0, 100])
ax[2].set_ylim([0, 0.25])

for i, (name, data) in enumerate(sims.items()):
    name = name.split('_')[-1]
    ls, c, m, ms, fs = plot_params[i]
    mark_every = data.times.shape[0] // 20
    gammas = np.abs(data.gammas - gamma_an) / gamma_an
    theta_an = -omega_c * data.times
    theta_c = np.pi - np.unwrap(np.arctan2(data.positions[:, 1], data.positions[:, 0]))
    thetas = np.abs(theta_c - theta_an)
    ax[0].plot(data.positions[:, 0], data.positions[:, 1], c=c, label=name)
    ax[1].plot(data.times / Tc, gammas, c=c, marker=m, ms=ms, markevery=mark_every, fillstyle=fs, label=name)
    ax[2].plot(data.times / Tc, thetas, c=c, marker=m, ms=ms, markevery=mark_every, fillstyle=fs, label=name)

ax[0].legend()
ax[1].legend()
ax[2].legend()
plt.show()