#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pyforce import *
from scipy import constants
import os, sys, h5py

# ========== Program Settings =============
# =========================================
prog_set = ProgramSettings(
    'binary_collisions_program',
    dims=2,
    debug_verbosity=3,
    output_verbosity=1,
    print_step=1,
    num_threads=1,
    dy=1.0E-6,
    halo=1
)

# ========== Simulation Params ============
# =========================================
sim_param = SimulationParameters(
    total_time=5.0E-12,
    dt=1.0E-15,
    cfl=1.0,
    enable_cfl=False
)

# =============    Mesh   =================
# =========================================
x_range = [0.0, 1.0e-6]
y_range = [0.0, 1.0e-6]
x = np.linspace(x_range[0], x_range[1], 2)
y = np.linspace(y_range[0], y_range[1], 2)

mesh_param = Mesh('simulation_mesh', xs=x, zs=y)

# ======= Materials and geometry ==========
# =========================================
margin = 1.0e-9
thebox = Polygon(
    'box',
    points=[[x_range[0] + margin, y_range[0] + margin],
            [x_range[1] - margin, y_range[0] + margin],
            [x_range[1] - margin, y_range[1] - margin],
            [x_range[0] + margin, y_range[1] - margin]])

# =========== IO parameters ===============
# =========================================
io_param = IO()
io_param.output_file_path = 'output'
io_param.timer_write_interval = Interval(value=1000, units='step_count')

io_param.checkpoint_enable = 1
io_param.checkpoint_file_path = 'checkpoints'
io_param.checkpoint_frequency = Interval(value=1.0, units='wall_time')
io_param.checkpoint_total_number = 1

# ================ EM =====================
# =========================================
em_param = Electromagnetics(solver_type='none')
mesh_param.finalize(prog_set, em_param)

# =========== Particle pusher =============
# =========================================
particle_pusher = ParticlePusher()
particle_pusher.enable_particle_pusher = False

# ============= Particles =================
# =========================================
# particle_group1 shows how to initialize group by list, file, and geometry
particle_group1 = ParticleGroup('carbon1')
particle_group1.mass = 12.011 * constants.m_p + 6 * constants.m_e
particle_group1.charge = 6
particle_group1.atomic_number = 6
particle_group1.max_count = int(1e5)
particle_group1.num_subcycling_groups = 1
particle_group1.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, 'step_count'), disable=True)
particle_group1.chaining_mesh_params = ChainingMeshParams(update_interval=10, sort_interval=100)
particle_group1.push_interval = Interval(value=1,  units="step_count")

particle_geometry1 = ParticleGeometry(
    thebox,
    ppc=[100, 100],
    density=1.0E26,
    temperature=[1.0E3, 1.0E3, 1.0E3],
    velocity=[0.0, 0.0, 0.0]
)
particle_group1.init_params = [particle_geometry1]

# particle_group2 shows how to initialize group without creating any particles
particle_group2 = ParticleGroup('carbon2')
particle_group2.mass = 12.011 * constants.m_p + 6 * constants.m_e
particle_group2.charge = 6
particle_group2.atomic_number = 6
particle_group2.max_count = int(1e5)
particle_group2.num_subcycling_groups = 1
particle_group2.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, 'step_count'), disable=True)
particle_group2.chaining_mesh_params = ChainingMeshParams(update_interval=10, sort_interval=100)
particle_group2.push_interval = Interval(value=1,  units="step_count")

# same geometry except with twice the temperature
particle_geometry2 = ParticleGeometry(
    thebox,
    ppc=[100, 100],
    density=1.0E26,
    temperature=[250.0, 250.0, 250.0],
    velocity=[0.0, 0.0, 0.0]
)
particle_group2.init_params = [particle_geometry2]

# ========= Timers and Metrics ===========
# =========================================
timer1 = Timer('timer_wall_time', Interval(value=1.0, units='wall_time'))
metric1 = IntegratedMetric('integrated', Interval(value=1, units='step_count'))
metric2 = PlasmaContour('plasma', Interval(value=1.0e-9, units='sim_time'), 'simulation_mesh',
                        ['density', 'temperature', 'momentum'])

# ============= Collisions ================
# =========================================
Coulomb_log = 10.0  # constant for all collisions

collision1 = Collision(group1_name='carbon1', interaction_type='binary')
collision1.channels = [Coulomb(Coulomb_log=Coulomb_log)]

collision2 = Collision(group1_name='carbon2', interaction_type='binary')
collision2.channels = [Coulomb(Coulomb_log=Coulomb_log)]

collision3 = Collision(group1_name='carbon1', group2_name='carbon2', interaction_type='binary')
collision3.channels = [Coulomb(Coulomb_log=Coulomb_log)]

# ================ Input ==================
# =========================================
global_input = Input(
    program_settings=prog_set,
    simulation_params=sim_param,
    meshes=mesh_param,
    io_params=io_param,
    electromagnetics=em_param,
    particle_groups=[particle_group1, particle_group2],
    particle_pusher=particle_pusher,
    collisions=[collision1, collision2, collision3],
    metrics=[metric1, metric2]
)

# ========= Compile and Run ===============
# =========================================
sim_dts = {0: [2.5e-15, r'$\nu\Delta t = 0.005$', '-'],
           1: [5.0e-15, r'$\nu\Delta t = 0.01$', '--'],
           2: [5.0e-14, r'$\nu\Delta t = 0.1$', ':'],
           3: [2.5e-13, r'$\nu\Delta t = 0.5$', '-.']}

run_ids = [2] #[0, 1, 2, 3]

for run_id in run_ids:

    io_param.output_file_path = 'output/Tequil_' + str(run_id)
    if not os.path.isdir(io_param.output_file_path):
        os.mkdir(io_param.output_file_path)

    sim_param.dt = sim_dts[run_id][0]

    global_input.compile('collisions.toml', relocate=True, bin_name='tflink')
    global_input.run('tflink', 'collisions.toml')

# ========= Post-processing ===============
# =========================================

proc_ids = [2] #[0, 1, 2, 3]

# theory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eV_to_erg = 1.60218e-12
s_to_ps = 1.0e12
m3_to_cm3 = 1e6
kg_to_g = 1.0e3
K_to_eV = constants.k / constants.e
e_sc = 4.8e-10  # statcoulombs


# def equilibrationRate(Z, m_g, n_cc, logL, TA_eV, TB_eV):
#     coef = 8 * (8 * np.pi)**0.5/3
#     return coef * (Z * e_sc)**4 * logL * n_cc / (m_g * (eV_to_erg * (TA_eV + TB_eV))**3)**0.5

def equilibrationRate(Z, m_g, n_cc, logL, TA_eV, TB_eV):
    coef = 8 * np.pi ** 0.5 / 3.0
    Ts = 0.5 * (TA_eV + TB_eV) * eV_to_erg
    return coef * (Z * e_sc) ** 4 * logL * n_cc / (m_g ** 0.5 * Ts ** 1.5)


T10 = particle_geometry1.temperature[0]
T20 = particle_geometry2.temperature[0]

mu = equilibrationRate(particle_group1.charge, particle_group1.mass * kg_to_g,
                       particle_geometry1.density / m3_to_cm3, Coulomb_log, T10, T20)
print(f'Equilibration time = 1/mu ={1 / mu:.2e}')

ntheory_steps = 50
time_theory = np.linspace(0, 5e-12, ntheory_steps)
temperature_C1_theory = np.zeros(ntheory_steps)
temperature_C2_theory = np.zeros(ntheory_steps)

for i in range(ntheory_steps):
    temperature_C1_theory[i] = 0.5 * (T10 + T20) + 0.5 * (T10 - T20) * np.exp(-mu * time_theory[i])
    temperature_C2_theory[i] = 0.5 * (T10 + T20) + 0.5 * (T20 - T10) * np.exp(-mu * time_theory[i])

time_theory *= s_to_ps

# plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plt.style.use('triforce.mplstyle')
good_colors = ['#db6d00', '#006ddb', '#920000', '#52a736', '#9B30FF']

fig, ax = plt.subplots(1, 2, figsize=(8, 4))

ax[0].set_ylabel('Temperature (eV)')
ax[0].set_ylim([600, 1000])
ax[1].set_ylim([250, 650])

for axis in ax:
    axis.set_xlabel('Time (ps)')
    axis.set_xlim([0, 5])
    axis.grid()

ax[0].plot(time_theory, temperature_C1_theory, '-k', label='Fluid model')
ax[1].plot(time_theory, temperature_C2_theory, '-k', label='Fluid model')

for proc_id in proc_ids:
    output_dir = 'output/Tequil_' + str(proc_id) + '/'
    data1 = np.loadtxt(output_dir + particle_group1.name + '_integrated.csv', skiprows=1)
    data2 = np.loadtxt(output_dir + particle_group2.name + '_integrated.csv', skiprows=1)

    # column order
    # Time Mass KE Tavg Tx Ty Tz px py pz vx vy vz

    time1 = data1[:, 0] * 1e12  # time in ps
    Temp1 = data1[:, 3]  # temperature in eV
    Temp2 = data2[:, 3]

    ax[0].plot(time1, Temp1, linestyle=sim_dts[proc_id][2], c=good_colors[proc_id], label=sim_dts[proc_id][1])
    ax[1].plot(time1, Temp2, linestyle=sim_dts[proc_id][2], c=good_colors[proc_id], label=sim_dts[proc_id][1])

ax[0].legend()
fig.tight_layout(pad=0.4, rect=[0, 0, 1, 1])
plt.show()

# plt.style.use('triforce.mplstyle')

# fig, ax = plt.subplots(1,2,figsize=(8,3))
#
# ax[0].set_ylabel('Temperature (eV)')
# ax[0].plot(time1, Temp1, '-b')
# ax[0].plot(time1, Temp2, '-r')
#
# ax[1].set_ylabel('Kinetic energy (J)')
# ax[1].plot(time1, KE1, '-b')
# ax[1].plot(time1, KE2, '-r')
# ax[1].plot(time1, KE1+KE2, '-k')
#
# for axis in ax:
#     axis.set_xlim([0,5])
#     axis.set_xlabel('Time (ps)')
#     axis.grid()
#
# fig.tight_layout(pad=0.5, rect=[0,0,1,1])
# plt.show()
