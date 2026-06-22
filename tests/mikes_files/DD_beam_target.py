#!/usr/bin/env python3
"""
Created by M Lavell on 07/19/24

Simulate DD and DT thermonuclear fusion.

Original test from Wu et al.
"""

import numpy as np
import matplotlib.pyplot as plt
from pyforce import *
from scipy import constants
import os, sys, h5py, math

# load functions and mplstyle
home = os.getenv("TRIFORCE_ROOT")
plt.style.use(home + '/examples/triforce.mplstyle')
sys.path.append(home + '/examples/physics_tests')
from physics_utilities import *
from collision_diagnostics import *

cross_section_dir = home + "/libs/common/data/cross_sections/"

# ========== Program Settings =============
# =========================================
prog_set = ProgramSettings(
    'beam_target_fusion',
    dims=2,
    debug_verbosity=3,
    output_verbosity=1,
    print_step=10,
    num_threads=1,
    dy=1.0E-6,
    halo=1
)

# ========== Simulation Params ============
# =========================================
sim_param = SimulationParameters(
    total_time=5.0E-15,
    dt=0.25E-15,
    cfl=1.0,
    enable_cfl=False,
    particle_buffer_size=int(1e6)
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
thebox = Polygon(
    'box',
    points=[[x_range[0], z_range[0]],
            [x_range[1], z_range[0]],
            [x_range[1], z_range[1]],
            [x_range[0], z_range[1]]])

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
particle_pusher = ParticlePusher(
    particle='particle2d3v',
    shape='quad',
    interpolator='energy_conserving',
    velocity_integrator='ballistic',
    trajectory_integrator='leapfrog1d')

# =========== Particle Groups =============
# =========================================
T_eV = 0.0257
ppc = [100, 100]
v_beam = 2.98304645e+07

plasma_bcs = ParticleBoundaryConditions(x_min='periodic',
                                        x_max='periodic',
                                        z_min='periodic',
                                        z_max='periodic')
push_interval = Interval(value=1, units='step_count')

# define cold ion group ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D_group = ParticleGroup('deuterium')
D_group.mass = 2.0014 * constants.atomic_mass
D_group.charge = 1
D_group.atomic_number = 1
D_group.max_count = int(1e6)
D_group.num_subcycling_groups = 1
D_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, 'step_count'), disable=True)
D_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
D_group.boundary_conditions = plasma_bcs
D_group.push_interval = push_interval

D_beam = ParticleGeometry(
    thebox,
    ppc=ppc,
    density=1.0E30,
    temperature=[T_eV, T_eV, T_eV],
    velocity=[0.0, 0.0, v_beam])

D_target = ParticleGeometry(
    thebox,
    ppc=ppc,
    density=1.0E30,
    temperature=[T_eV, T_eV, T_eV],
    velocity=[0.0, 0.0, 0.0])

D_group.init_params = [D_beam, D_target]

# define empty particle groups for fusion products
neutron_group = ParticleGroup('neutron')
neutron_group.mass = constants.m_n
neutron_group.charge = 0
neutron_group.max_count = int(1e6)
neutron_group.num_subcycling_groups = 1
neutron_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, 'step_count'), disable=True)
neutron_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
neutron_group.init_params = [ParticleNullInit('empty_initialization')]
neutron_group.boundary_conditions = plasma_bcs
neutron_group.push_interval = push_interval

helium3_group = ParticleGroup('helium3')
helium3_group.mass = 3.016029 * constants.atomic_mass
helium3_group.charge = 2
helium3_group.atomic_number = 2
helium3_group.max_count = int(1e6)
helium3_group.num_subcycling_groups = 1
helium3_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, 'step_count'), disable=True)
helium3_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
helium3_group.init_params = [ParticleNullInit('empty_initialization')]
helium3_group.boundary_conditions = plasma_bcs
helium3_group.push_interval = push_interval

# ========= Timers and Metrics ===========
# =========================================
metric1 = IntegratedMetric('integrated', Interval(value=1, units='step_count'))
DD_particle_metric = ParticleData("particle_data", Interval(value=1, units="step_count"),
                                  group_ids=[0, 1, 2], sample_fraction=1.0)

# ============= Collisions ================
# =========================================

DD_collisions = Collision(group1_name='deuterium', interaction_type='binary', disable=False)

DD_n3He_fusion = BinaryChannel('fusion', emission_type='isotropic')
DD_n3He_fusion.energy_gain = 3.269e6
DD_n3He_fusion.production_multiplier = 1.0e10
DD_n3He_fusion.product_group_ids = [1, 2]
DD_n3He_fusion.cross_section_filename = cross_section_dir + 'DD_n3He_BH_keV_mbarn.csv'
DD_n3He_fusion.scale_table_to_eV = 1.0e3
DD_n3He_fusion.scale_table_to_barns = 1.0e-3

DD_n3He_fusion.emission_type = 'fusion_polyfit'  # 'isotropic'
DD_n3He_fusion.differential_cs_filename = cross_section_dir + 'DD_n3He_anisotropy_poly_coefs.dat'

DD_collisions.channels = [DD_n3He_fusion]

# ================ Input ==================
# =========================================
DD_inputs = Input(
    program_settings=prog_set,
    simulation_params=sim_param,
    meshes=mesh_param,
    io_params=io_param,
    electromagnetics=em_param,
    particle_groups=[D_group, neutron_group, helium3_group],
    particle_pusher=particle_pusher,
    collisions=[DD_collisions],
    metrics=[metric1, DD_particle_metric])

# ========= Compile and Run ===============
# =========================================
run_sim = True

v_beam = [2.99777470e+06, 1.49709209e+07, 2.98304645e+07]
n_beams = len(v_beam)
labels = [r'(a) $u=0.01c$', r'(b) $u=0.05c$', r'(c) $u=0.1c$']

if run_sim:

    for i in range(n_beams):
        # change DD output directories and temperatures
        io_param.output_file_path = f'output/beam_target_{i}'
        os.makedirs(io_param.output_file_path, exist_ok=True)

        D_beam.velocity = [0.0, 0.0, v_beam[i]]

        DD_inputs.compile('DD_beam_target.toml', relocate=True, bin_name='tflink')
        DD_inputs.run('tflink', 'DD_beam_target.toml')

# ========= Post-processing ===============
# =========================================
print('\nPost-processing...\n')

post_process = True

n_theta_bins = 100
theta_bin_width = 180 / n_theta_bins
theta_bin_centers = np.arange(0.5 * theta_bin_width, 180 + 0.5 * theta_bin_width, theta_bin_width)

n_energy_bins = 100
energy_bin_ranges = [[2.0e6, 3.0e6],
                     [1.0e6, 6.0e6],
                     [0.1e6, 14.0e6]]
max_energy_plot_ranges = [4.0e6, 8.0e6, 16.0e6]


# Function for binning neutron energies
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


if post_process:

    ###########################################
    # start plot

    fig, ax = plt.subplots(3, 1, figsize=(8, 8))

    for i in range(len(labels)):
        ax[i].set_title(labels[i], fontsize=16, loc='left')

    for axis in ax:
        axis.set_ylabel(r'$\theta$ (Degree)')
        axis.set_ylim([0, 180])
        axis.tick_params(axis='x', pad=8)
        axis.grid()

    simulation_volume = (x_range[1] - x_range[0]) * (z_range[1] - z_range[0]) * 1.0e-6
    mn_c2 = neutron_group.mass * constants.c ** 2

    n_theory_pts = 100
    theta_theory = np.linspace(0.0, 180.0, n_theory_pts)
    En_theta_theory = np.zeros(n_theory_pts)

    for i in range(n_beams):

        ###########################################
        # Theory solution from Ball 1962
        v_cm = v_beam[i] / 2.0  # * D_group.mass / (D_group.mass + D_group.mass)
        E_beam = calc_beam_energy_v(v_beam[i], D_group.mass)

        for j in range(n_theory_pts):
            theta_radian = theta_theory[j] * np.pi / 180
            En_theta_theory[j] = calc_En_thetaBeam(theta_radian, v_cm,
                                                   neutron_group.mass, helium3_group.mass,
                                                   DD_n3He_fusion.energy_gain * constants.e,
                                                   E_beam * constants.e)

        ax[i].plot(En_theta_theory / constants.e, theta_theory, '-', c=good_colors[1], label="Ball 1962")

        ###########################################
        # process output data

        print(f'  -- post-processing v={v_beam[i]} m/s')
        data_path = f'output/beam_target_{i}/'
        neutron_particle_data = h5py.File(data_path + 'neutron_particles_00000020.h5', 'r')

        min_energy = energy_bin_ranges[i][0]
        max_energy = energy_bin_ranges[i][1]
        energy_bin_width = (max_energy - min_energy) / n_energy_bins
        energy_bin_centers = np.arange(min_energy + 0.5 * energy_bin_width,
                                       max_energy + 0.5 * energy_bin_width,
                                       energy_bin_width)

        avg_En_per_theta = np.zeros(n_theta_bins)
        theta_bins = np.zeros(n_theta_bins)
        dyde_bins = np.zeros(n_energy_bins)

        for j in range(neutron_particle_data['neutron'].attrs['size']):
            momentum = neutron_particle_data['neutron']['momentum'][j]
            weight = neutron_particle_data['neutron']['weight'][j]

            theta_degrees = math.acos(momentum[2] / np.linalg.norm(momentum)) * 180.0 / np.pi
            neutron_energy_eV = (neutron_particle_data['neutron']['gamma'][j] - 1.0) * mn_c2 / constants.e

            closest_theta_bin, theta_bin_id = find_nearest(theta_bin_centers, theta_degrees)
            avg_En_per_theta[theta_bin_id] += weight * neutron_energy_eV
            theta_bins[theta_bin_id] += weight

            closest_energy_bin, energy_bin_id = find_nearest(energy_bin_centers, neutron_energy_eV)
            dyde_bins[energy_bin_id] += weight

        # average energy per degree and normalized yield
        avg_En_per_theta /= theta_bins
        dyde_bins /= max(dyde_bins)

        ax[i].set_xlim([0.0, max_energy_plot_ranges[i]])
        ax[i].plot(avg_En_per_theta, theta_bin_centers, marker='o', linestyle='none', fillstyle='none',
                   ms=10, markeredgewidth=1.5, markeredgecolor=good_colors[0], label='simulated')

        ax_twin = ax[i].twinx()
        ax_twin.set_ylabel(r'$\mathrm{d}Y/\mathrm{d}E$ (arb.)', color=good_colors[2])
        ax_twin.set_ylim([0.0, 1.0])
        ax_twin.plot(energy_bin_centers, dyde_bins, '-', c=good_colors[2])
        ax_twin.tick_params(axis='y', labelcolor=good_colors[2])

    ax[0].legend(loc=4)
    fig.tight_layout(pad=0.1, rect=[0, 0, 1, 1])

    plt.show()
