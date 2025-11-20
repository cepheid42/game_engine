#/usr/bin/env python3
"""
Created by M Lavell on 06/03/24

Simulate the ionization of sold aluminum by a 500 eV electron beam.

Test originally performed by Perez et al. 2012
"""

import numpy as np
from pyforce import *
from scipy import constants
import os, sys
import matplotlib.pyplot as plt

# load functions and mplstyle
home = os.getenv("TRIFORCE_ROOT")
plt.style.use(home + '/examples/triforce.mplstyle')
sys.path.append(home + '/examples/physics_tests')
from physics_utilities import *
from collision_diagnostics import *

# ========== Program Settings =============
# =========================================
prog_set = ProgramSettings(
    'burn1d_program',
    dims=2,
    debug_verbosity=3,
    output_verbosity=1,
    print_step=100,
    num_threads=1,
    dy=1.0E-2,
    halo=1
)

# ========== Simulation Params ============
# =========================================
sim_param = SimulationParameters(
    total_time=3.18e-15,
    dt=5.0e-18,
    cfl=1.0,
    enable_cfl=False,
    particle_buffer_size=int(10e6)
)

# =============    Mesh   =================
# =========================================
x_range = [0.0, 1.0e-6]
z_range = [0.0, 1.0e-6]
x = np.linspace(x_range[0], x_range[1], 2)
z = np.linspace(z_range[0], z_range[1], 2)
mesh_param = Mesh('simulation_mesh', xs=x, zs=z)

# ======= Materials and geometry ==========
# =========================================
sim_geometry = Polygon(
    "sim_geometry",
    points=[[x_range[0], z_range[0]],
            [x_range[0], z_range[1]],
            [x_range[1], z_range[1]],
            [x_range[1], z_range[0]]])

# =========== IO parameters ===============
# =========================================
io_param = IO()
io_param.output_file_path = "output"
io_param.timer_write_interval = Interval(value=1000, units='step_count')

io_param.checkpoint_enable = 1
io_param.checkpoint_file_path = "checkpoints"
io_param.checkpoint_frequency = Interval(value=1.0, units="wall_time")
io_param.checkpoint_total_number = 1

# ================ EM =====================
# =========================================
em_param = Electromagnetics(solver_type='none')
mesh_param.finalize(prog_set, em_param)

# =========== Particle pusher =============
# =========================================
particle_pusher = ParticlePusher(
    particle='particle2d3v',
    shape='quad',
    interpolator='energy_conserving',
    velocity_integrator='ballistic',
    trajectory_integrator='leapfrog1d')

# =========== Particle Groups =============
# =========================================

E_beam = 500  # eV
v_beam = calc_beam_velocity_EeV(E_beam, constants.m_e)

e_ppc_init = [100, 10]
i_ppc_init = [20, 20]
max_particles = int(1e6)

T_cold = 0.0

apm_target = [10, 10]
apm_interval = 100
apm_min_cnt = int(0.8 * e_ppc_init[0] * e_ppc_init[1])
apm_disabled = True

plasma_bcs = ParticleBoundaryConditions(x_min='periodic',
                                        x_max='periodic',
                                        z_min='periodic',
                                        z_max='periodic')
push_interval = Interval(value=1, units='step_count')

# define electron beam ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
electron_beam = ParticleGroup("electron_beam")
electron_beam.mass = constants.m_e
electron_beam.charge = -1
electron_beam.max_count = max_particles
electron_beam.num_subcycling_groups = 1
electron_beam.apm_params = ApmParams(target_ppc=apm_target, interval=Interval(apm_interval, "step_count"),
                                     min_count=apm_min_cnt, disable=apm_disabled)
electron_beam.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
electron_beam.boundary_conditions = plasma_bcs
electron_beam.push_interval = push_interval

electron_beam_geom = ParticleGeometry(
    sim_geometry,
    ppc=e_ppc_init,
    density=1.1e27,
    temperature=[T_cold, T_cold, T_cold],
    velocity=[0.0, 0.0, v_beam]
)
electron_beam.init_params = [electron_beam_geom]

# define electron products ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
electron_products = ParticleGroup("electron_products")
electron_products.mass = constants.m_e
electron_products.charge = -1
electron_products.max_count = max_particles
electron_products.num_subcycling_groups = 1
electron_products.apm_params = ApmParams(target_ppc=apm_target, interval=Interval(apm_interval, "step_count"),
                                         min_count=apm_min_cnt, disable=apm_disabled)
electron_products.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
electron_products.init_params = [ParticleNullInit('empty_initialization')]
electron_products.boundary_conditions = plasma_bcs
electron_products.push_interval = push_interval

# define neutral aluminum ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
neutral_aluminum = ParticleGroup("neutral_aluminum")
neutral_aluminum.mass = 26.9815 * constants.atomic_mass + 13.0 * constants.m_e
neutral_aluminum.charge = 0
neutral_aluminum.max_count = max_particles
neutral_aluminum.num_subcycling_groups = 1
neutral_aluminum.apm_params = ApmParams(target_ppc=apm_target, interval=Interval(apm_interval, "step_count"),
                                        min_count=apm_min_cnt, disable=apm_disabled)
neutral_aluminum.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
neutral_aluminum.boundary_conditions = plasma_bcs
neutral_aluminum.push_interval = push_interval

neutral_aluminum_geom = ParticleGeometry(
    sim_geometry,
    ppc=i_ppc_init,
    density=6.6e28,  # #/m^-3, solid density neutral Al
    temperature=[T_cold, T_cold, T_cold],
    velocity=[0.0, 0.0, 0.0]
)
neutral_aluminum.init_params = [neutral_aluminum_geom]

# define electron products ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aluminum_ion = ParticleGroup("aluminum_ion")
aluminum_ion.mass = 26.9815 * constants.atomic_mass + 12.0 * constants.m_e
aluminum_ion.charge = 1
aluminum_ion.max_count = max_particles
aluminum_ion.num_subcycling_groups = 1
aluminum_ion.apm_params = ApmParams(target_ppc=apm_target, interval=Interval(apm_interval, "step_count"),
                                    min_count=apm_min_cnt, disable=apm_disabled)
aluminum_ion.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
aluminum_ion.boundary_conditions = plasma_bcs
aluminum_ion.push_interval = push_interval

aluminum_ion.init_params = [ParticleNullInit('empty_initialization')]

# ============= Collisions ================
# =========================================

cs_dir = home + "/libs/common/data/cross_sections/"

eAl_collisions = BinaryCollision(group1_name="electron_beam", group2_name="neutral_aluminum")
eAl_collisions.channels = [Ionization(ionization_energy=5.9858,
                                      product_group_ids=[1, 3],  # first is electron product, second is ion
                                      production_multiplier=1.0,
                                      rejection_multiplier=1.0,
                                      cross_section_filename=cs_dir + "electron_Al0_impactIonz_crossSection.txt",
                                      scale_table_to_eV=1.0,
                                      scale_table_to_barns=1.0e28,
                                      emission_type="isotropic")]

# ========= Timers and Metrics ===========
# =========================================
metric_interval = 20
integrated_metric = IntegratedMetric("integrated", Interval(value=metric_interval, units="step_count"))
plasma_metric = PlasmaContour("plasma", Interval(value=metric_interval, units="step_count"), "simulation_mesh",
                              ["density", "temperature", "momentum"])
collision_metric = CollisionMetric("collisions", Interval(value=metric_interval, units="step_count"))

# ================ Input ==================
# =========================================
global_input = Input(
    program_settings=prog_set,
    simulation_params=sim_param,
    meshes=mesh_param,
    io_params=io_param,
    electromagnetics=em_param,
    particle_groups=[electron_beam, electron_products, neutral_aluminum, aluminum_ion],
    particle_pusher=particle_pusher,
    collisions=[eAl_collisions],
    metrics=[integrated_metric, plasma_metric, collision_metric]
)

# ========= Compile and Run ===============
# =========================================
# sim_id: production factor, labels
scan_params = {0: [1.0, r'$M_{\mathrm{p}}=10^0$', '-', 's'],
               1: [10.0, r'$M_{\mathrm{p}}=10^1$', '--', 'd'],
               2: [100.0, r'$M_{\mathrm{p}}=10^2$', '-.', 'o'],
               3: [1000.0, r'$M_{\mathrm{p}}=10^3$', '-.', '^']}

sim_ids = [0, 1, 2, 3]

for sim_id in sim_ids:
    # make directory for output data
    io_param.output_file_path = "output/ionization_rate_" + str(sim_id)
    os.makedirs(io_param.output_file_path, exist_ok=True)

    eAl_collisions.channels[0].production_multiplier = scan_params[sim_id][0]

    global_input.compile('ionization_rate.toml', relocate=True, bin_name='tflink')
    global_input.run('tflink', 'ionization_rate.toml')

# ========= Post-processing ===============
# =========================================

# plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig, ax = plt.subplots(1, 2, figsize=(8, 4))

ax[0].set_title('(a) Mean ion charge', fontsize=16, loc='left')
ax[1].set_title('(b) Particle count', fontsize=16, loc='left')
ax[0].set_ylim([0.0, 0.5])
# ax[1].set_ylim([])

for axis in ax:
    axis.set_xlabel('Time (ps)')
    axis.set_xlim([0, 6.0])
    axis.tick_params(axis='x', pad=8)
    axis.grid()

# Theory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
volume = (x_range[1] - x_range[0]) * (z_range[1] - z_range[1]) * prog_set.dy
omega_p_inv = 0.53e-15
sigma = 1.428e-20

nsteps_thry = 100
time_thry = np.linspace(0, 6, nsteps_thry) * omega_p_inv
thry_ion_charge = np.zeros(nsteps_thry)

for i in range(nsteps_thry):
    thry_ion_charge[i] = 1.0 - np.exp(-v_beam * electron_beam_geom.density * sigma * time_thry[i])

ax[0].plot(time_thry / omega_p_inv, thry_ion_charge, '-k', label="Theory")

for sim_id in sim_ids:
    output_dir = 'output/ionization_rate_' + str(sim_id) + '/'

    data_electron_beam = np.loadtxt(output_dir + electron_beam.name + '_integrated.csv', skiprows=1)
    data_electron_products = np.loadtxt(output_dir + electron_products.name + '_integrated.csv', skiprows=1)

    data_neutral_aluminum = np.loadtxt(output_dir + neutral_aluminum.name + '_integrated.csv', skiprows=1)
    data_aluminum_ion = np.loadtxt(output_dir + aluminum_ion.name + '_integrated.csv', skiprows=1)

    # column order
    # Time Mass KE Tavg Tx Ty Tz px py pz vx vy vz n_particles

    m0_Al0 = data_neutral_aluminum[0, 1] # / neutral_aluminum.mass  # initial number of aluminum atoms
    mean_ion_charge = data_aluminum_ion[:, 1] / m0_Al0

    N_ttl = data_electron_beam[:, 14] + data_electron_products[:, 14] + \
            data_neutral_aluminum[:, 14] + data_aluminum_ion[:, 14]

    time = data_aluminum_ion[:, 0] / omega_p_inv

    skip = 1

    ax[0].plot(time[::skip], mean_ion_charge[::skip], c=good_colors[sim_id], label=scan_params[sim_id][1],
               linestyle='none', marker=scan_params[sim_id][3], mec=good_colors[sim_id], mew=1.5, ms=5, fillstyle="none")
    ax[1].plot(time[::skip], N_ttl[::skip], c=good_colors[sim_id], label=scan_params[sim_id][1],
               linestyle='none', marker=scan_params[sim_id][3], mec=good_colors[sim_id], mew=1.5, ms=5, fillstyle="none")

ax[0].legend(loc=4)
fig.tight_layout(pad=0.4, rect=[0, 0, 1, 1])

plt.show()
