#!/usr/bin/env python3
"""
Created by M Lavell on 02/10/24
Simulate DD and DT thermonuclear fusion.
"""

import numpy as np
import matplotlib.pyplot as plt
from pyforce import *
from scipy import constants
import os, sys, h5py

from collision_diagnostics import *

home = os.getenv('HOME')
tf_dir = home + '/TriForce/triforce/'
cs_dir = tf_dir + 'libs/common/data/cross_sections/'

# ========== Program Settings =============
# =========================================
prog_set = ProgramSettings(
    'thermonuclear_fusion',
    dims=2,
    debug_verbosity=3,
    output_verbosity=1,
    print_step=100,
    num_threads=1,
    dy=1.0E-6
)

# ========== Simulation Params ============
# =========================================
sim_param = SimulationParameters(
    total_time=1.0E-12,
    dt=50.0E-15,
    cfl=1.0,
    enable_cfl=False
)

# =============    Mesh   =================
# =========================================
x_range = [0.0, 1.0e-6]
z_range = [0.0, 1.0e-6]
x = np.linspace(x_range[0], x_range[1], 3)
z = np.linspace(z_range[0], z_range[1], 3)

mesh_param = Mesh('simulation_mesh', x, z)

# ======= Materials and geometry ==========
# =========================================
thebox = Polygon(
    'box',
    points=[[x_range[0], z_range[0]],
            [x_range[1], z_range[0]],
            [x_range[1], z_range[1]],
            [x_range[0], z_range[1]]]
)

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

# =========== Particle pusher =============
# =========================================
particle_pusher = ParticlePusher()
particle_pusher.enable_particle_pusher = False
particle_pusher.push_max_iterations = 50

# =========== Particle Groups =============
# =========================================
T_eV = 5.0e3

# define cold ion group ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D_group = ParticleGroup('deuterium')
D_group.mass = 2.0014 * constants.atomic_mass
D_group.charge = 1
D_group.atomic_number = 1
D_group.max_count = int(1e6)
D_group.num_subcycling_groups = 1
D_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, 'step_count'), disable=True)
D_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)

D_group_geom = ParticleGeometry(
    'box_geometry',
    thebox,
    ppc=[50, 50],
    density=1.0E26,
    temperature=[T_eV, T_eV, T_eV],
    velocity=[0.0, 0.0, 0.0]
)
D_group.init_params = [D_group_geom]

# define electron beam group ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
T_group = ParticleGroup('tritium')
T_group.mass = 3.01605 * constants.atomic_mass
T_group.charge = 1
T_group.atomic_number = 1
T_group.max_count = int(1e6)
T_group.num_subcycling_groups = 1
T_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, 'step_count'), disable=True)
T_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)

T_group_geom = ParticleGeometry(
    'box_geometry',
    thebox,
    ppc=[50, 50],
    density=1.0E26,
    temperature=[T_eV, T_eV, T_eV],
    velocity=[0.0, 0.0, 0.0]
)
T_group.init_params = [T_group_geom]

# define empty particle groups for fusion products
neutron_group = ParticleGroup('neutron')
neutron_group.mass = constants.m_n
neutron_group.charge = 0
neutron_group.max_count = int(1e6)
neutron_group.num_subcycling_groups = 1
neutron_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, 'step_count'), disable=True)
neutron_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
neutron_group.init_params = [ParticleNullInit('empty_initialization')]

helium3_group = ParticleGroup('helium3')
helium3_group.mass = 3.016029 * constants.atomic_mass
helium3_group.charge = 2
helium3_group.atomic_number = 2
helium3_group.max_count = int(1e6)
helium3_group.num_subcycling_groups = 1
helium3_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, 'step_count'), disable=True)
helium3_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
helium3_group.init_params = [ParticleNullInit('empty_initialization')]

helium4_group = ParticleGroup('helium4')
helium4_group.mass = 4.002603 * constants.atomic_mass
helium4_group.charge = 2
helium4_group.atomic_number = 2
helium4_group.max_count = int(1e6)
helium4_group.num_subcycling_groups = 1
helium4_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, 'step_count'), disable=True)
helium4_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
helium4_group.init_params = [ParticleNullInit('empty_initialization')]

# ========= Timers and Metrics ===========
# =========================================
timer1 = Timer('timer_wall_time', Interval(value=1.0, units='wall_time'))
metric1 = IntegratedMetric('integrated', Interval(value=1, units='step_count'))
DD_particle_metric = ParticleData("particle_data", Interval(value=1, units="step_count"),
                                    group_ids=[0, 1, 2], sample_fraction=1.0)

DT_particle_metric = ParticleData("particle_data", Interval(value=1, units="step_count"),
                                  group_ids=[0, 1, 2, 3], sample_fraction=1.0)
# ============= Collisions ================
# =========================================

# deuterium-deuterium fusion
DD_collisions = Collision(group1_name='deuterium', interaction_type='binary', disable=False)

DD_n3He_fusion = BinaryChannel('fusion', emission_type='isotropic')
DD_n3He_fusion.energy_gain = 3.269e6
DD_n3He_fusion.production_multiplier = 1.0e10
DD_n3He_fusion.product_group_ids = [1, 2]
DD_n3He_fusion.cross_section_filename = cs_dir + 'DD_n3He_BH_keV_mbarn.csv'
DD_n3He_fusion.scale_table_to_eV = 1.0e3
DD_n3He_fusion.scale_table_to_barns = 1.0e-3

# DD_pT_fusion = BinaryChannel('fusion', emission_type='isotropic')
# DD_pT_fusion.energy_gain = 4.03e6
# DD_pT_fusion.production_multiplier = 1.0e10
# DD_pT_fusion.product_group_ids = []
# DD_pT_fusion.cross_section_filename = '../../libs/common/data/cross_sections/DD_pT_BH_keV_mbarn.csv'
# DD_pT_fusion.scale_table_to_eV = 1.0e3
# DD_pT_fusion.scale_table_to_barns = 1.0e-31

DD_collisions.channels = [DD_n3He_fusion]

# deuterium-tritium fusion
DT_collisions = Collision(group1_name='deuterium', group2_name='tritium', interaction_type='binary', disable=False)
DT_fusion = BinaryChannel('fusion', emission_type='isotropic')
DT_fusion.energy_gain = 17.589e6
DT_fusion.production_multiplier = 1.0e10
DT_fusion.product_group_ids = [2, 3]
DT_fusion.cross_section_filename = cs_dir + 'DT_n4He_BH_keV_mbarn.csv'
DT_fusion.scale_table_to_eV = 1.0e3
DT_fusion.scale_table_to_barns = 1.0e-3
DT_collisions.channels = [DT_fusion]

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
    metrics=[metric1, DD_particle_metric]
)

DT_inputs = Input(
    program_settings=prog_set,
    simulation_params=sim_param,
    meshes=mesh_param,
    io_params=io_param,
    electromagnetics=em_param,
    particle_groups=[D_group, T_group, neutron_group, helium4_group],
    particle_pusher=particle_pusher,
    collisions=[DT_collisions],
    metrics=[metric1, DT_particle_metric]
)

# ========= Compile and Run ===============
# =========================================
run_DD_sim = True
run_DT_sim = True

# for each temperature value, run DD and DT simulations
temperatures_keV = [5, 10, 19]

for T_keV in temperatures_keV:
    T_eV = T_keV * 1e3

    if run_DD_sim:
        # change DD output directories and temperatures
        io_param.output_file_path = 'output/DD_' + str(T_keV) + 'keV'
        if not os.path.isdir(io_param.output_file_path):
            os.mkdir(io_param.output_file_path)

        D_group_geom.temperature = [T_eV, T_eV, T_eV]

        DD_inputs.compile('DD_thermonuclear.toml', relocate=True, bin_name='tflink-collisions')
        DD_inputs.run('tflink-collisions', 'DD_thermonuclear.toml')

    if run_DT_sim:
        # change DT output directories and temperatures
        io_param.output_file_path = 'output/DT_' + str(T_keV) + 'keV'
        if not os.path.isdir(io_param.output_file_path):
            os.mkdir(io_param.output_file_path)

        T_group_geom.temperature = [T_eV, T_eV, T_eV]

        DT_inputs.compile('DT_thermonuclear.toml', relocate=True, bin_name='tflink-collisions')
        DT_inputs.run('tflink-collisions', 'DT_thermonuclear.toml')

# ========= Post-processing ===============
# =========================================
print('\nPost-processing...\n')

DD_post_process = True
DT_post_process = True


###########################################
# functions give theory solutions (eq 12-14 from Higginson 2019)

# theory solution for reactivity from Bosch and Hale (1992) fitting parameters (eqs 12-14 and table 7)
def thermal_reactivity_BH1992(T, Bg, mrc2, Cparams):
    theta = T / (1.0 - T * (Cparams[1] + T * (Cparams[3] + T * Cparams[5])) /
                 (1.0 + T * (Cparams[2] + T * (Cparams[4] + T * Cparams[6]))))
    xsi = (Bg ** 2 / (4 * theta)) ** (1 / 3)
    reactivity = Cparams[0] * theta * np.sqrt(xsi / (mrc2 * T ** 3)) * np.exp(-3.0 * xsi)

    return reactivity


# theory solution for second moment of neutron spectrum (expressed as FWHM) is from
# Ballabio, Kallne, Gorini 1998, Relativistic calculation of fusion product spectra for thermonuclear plasmas
def calc_FWHM_Ballabio1998(T, alpha, w0):
    delta_w = alpha[0] * T ** (2 / 3) / (1.0 + alpha[1] * T ** alpha[2]) + alpha[3] * T
    whalf = w0 * (1.0 + delta_w) * np.sqrt(T)

    return whalf


if DD_post_process or DT_post_process:

    ###########################################
    # calculate theory solutions
    T_keV_theory = np.linspace(0.1, 20, 100)

    nD = D_group_geom.density
    nT = T_group_geom.density

    # fitting parameters for D+D -> n+He3 from table 7
    Bg = 31.3970  # sqrt(keV)
    mrc2 = 937814  # keV
    Cparams = [5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0]
    DD_reactivity_theory = thermal_reactivity_BH1992(T_keV_theory, Bg, mrc2, Cparams)  # units of cm3/s

    # fitting parameters for D+D -> n+He3 for T=[0,20], from table 3
    alpha = [1.7013e-3, 0.16888, 0.49, 7.9460e-4]
    omega0 = 82.542  # sqrt(keV), from table 1
    DD_FWHM_theory = calc_FWHM_Ballabio1998(T_keV_theory, alpha, omega0)

    # fitting parameters for D+T -> n+alpha from table 7
    Bg = 34.3827  # sqrt(keV)
    mrc2 = 1124656  # keV
    Cparams = [1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-6]
    DT_reactivity_theory = thermal_reactivity_BH1992(T_keV_theory, Bg, mrc2, Cparams)  # units of cm3/s

    # fitting parameters for D+T -> n+He4 for T=[0,20], from table 3
    alpha = [5.1068e-4, 7.6223e-3, 1.78, 8.7691e-5]
    omega0 = 177.259  # sqrt(keV), from table 1
    DT_FWHM_theory = calc_FWHM_Ballabio1998(T_keV_theory, alpha, omega0)

    ###########################################
    # start plot
    plt.style.use('triforce.mplstyle')
    good_colors = ['#db6d00', '#006ddb', '#920000', '#52a736', '#9B30FF']

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))

    ax[0].set_xlabel(r'$T_{\mathrm{i}}$ (keV)')
    ax[0].set_ylabel(r'$\langle \sigma v \rangle$ (cm$^3$/s)')
    ax[0].set_xlim([0, 20])
    ax[0].set_yscale('log')
    ax[0].set_ylim([1.0e-22, 1e-15])
    ax[0].grid()

    ax[0].plot(T_keV_theory, DT_reactivity_theory, '-k')
    ax[0].plot(T_keV_theory, DD_reactivity_theory, '--k')

    ax[1].set_xlabel(r'$T_{\mathrm{i}}$ (keV)')
    ax[1].set_ylabel(r'FWHM (keV)')
    ax[1].set_xlim([0, 20])
    ax[1].set_yscale('log')
    ax[1].set_ylim([99.9, 1e3])  # set min at 99.9 instead of 100 so intermediate ticks not printed
    ax[1].yaxis.set_label_coords(-0.1, 0.5)
    ax[1].grid(which='both')

    ax[1].plot(T_keV_theory, DT_FWHM_theory, '-k')
    ax[1].plot(T_keV_theory, DD_FWHM_theory, '--k')

    ###########################################
    # process output data

    n_temps = len(temperatures_keV)

    DD_simulated_reactivity = np.zeros(n_temps)
    DT_simulated_reactivity = np.zeros(n_temps)

    DD_fwhm = np.zeros(n_temps)
    DT_fwhm = np.zeros(n_temps)

    simulation_volume = (x_range[1] - x_range[0]) * (z_range[1] - z_range[0]) * 1.0e-6
    delta_ij_DD = 1.0  # variable for getting reactivity
    delta_ij_DT = 0.0
    mn_c2 = neutron_group.mass * constants.c ** 2

    for i in range(n_temps):
        print(f'  -- post-processing T={temperatures_keV[i]} keV')

        # get reactivities from integrated metrics and fwhm from particle dtaa
        DD_path = 'output/DD_' + str(temperatures_keV[i]) + 'keV/'
        DT_path = 'output/DT_' + str(temperatures_keV[i]) + 'keV/'

        if DD_post_process:
            # DD reactivity
            print(f'    -- DD reactivity')
            neutron_integrated_data = np.loadtxt(DD_path + 'neutron_integrated.csv', skiprows=1)
            final_time = neutron_integrated_data[-1, 0]
            ttl_neutron_mass = neutron_integrated_data[-1, 1]
            dYdt = (ttl_neutron_mass / neutron_group.mass) / final_time
            DD_simulated_reactivity[i] = dYdt * (1.0 + delta_ij_DD) / (D_group_geom.density ** 2 * simulation_volume)
            DD_simulated_reactivity[i] *= 1e6  # convert from m3/s to cm3/s

            # DD fwhm
            print(f'    -- DD FWHM')
            neutron_particle_data = h5py.File(DD_path + 'neutron_particles_00000020.h5', 'r')
            E_ttl = 0.0
            w_ttl = 0.0
            dE2_ttl = 0.0
            for j in range(neutron_particle_data['neutron'].attrs['size']):
                weight = neutron_particle_data['neutron']['weight'][j]
                E_ttl += weight * (neutron_particle_data['neutron']['gamma'][j] - 1.0) * mn_c2 / constants.e * 1e-3
                w_ttl += weight
            E_ttl /= w_ttl  # weight average total energy
            for j in range(neutron_particle_data['neutron'].attrs['size']):
                energy = (neutron_particle_data['neutron']['gamma'][j] - 1.0) * mn_c2 / constants.e * 1e-3
                dE2_ttl += (E_ttl - energy) ** 2 * neutron_particle_data['neutron']['weight'][j]

            DD_fwhm[i] = (8.0 * np.log(2) * dE2_ttl / w_ttl) ** 0.5

        if DT_post_process:
            # DT reactivity
            print(f'    -- DT reactivity')
            neutron_integrated_data = np.loadtxt(DT_path + 'neutron_integrated.csv', skiprows=1)
            final_time = neutron_integrated_data[-1, 0]
            ttl_neutron_mass = neutron_integrated_data[-1, 1]
            dYdt = (ttl_neutron_mass / neutron_group.mass) / final_time
            DT_simulated_reactivity[i] = dYdt * (1.0 + delta_ij_DT) / (T_group_geom.density ** 2 * simulation_volume)
            DT_simulated_reactivity[i] *= 1e6  # convert from m3/s to cm3/s

            # DT fwhm
            print(f'    -- DT FWHM')
            neutron_particle_data = h5py.File(DT_path + 'neutron_particles_00000020.h5', 'r')
            E_ttl = 0.0
            w_ttl = 0.0
            dE2_ttl = 0.0
            for j in range(neutron_particle_data['neutron'].attrs['size']):
                weight = neutron_particle_data['neutron']['weight'][j]
                E_ttl += weight * (neutron_particle_data['neutron']['gamma'][j] - 1.0) * mn_c2 / constants.e * 1e-3
                w_ttl += weight
            E_ttl /= w_ttl  # weight average total energy
            for j in range(neutron_particle_data['neutron'].attrs['size']):
                energy = (neutron_particle_data['neutron']['gamma'][j] - 1.0) * mn_c2 / constants.e * 1e-3
                dE2_ttl += (E_ttl - energy) ** 2 * neutron_particle_data['neutron']['weight'][j]

            DT_fwhm[i] = (8.0 * np.log(2) * dE2_ttl / w_ttl) ** 0.5

    ax[0].plot(temperatures_keV, DD_simulated_reactivity, marker='o', c=good_colors[3], linestyle='none',
               ms=10, markeredgewidth=1, markeredgecolor='k', label='D-D')
    ax[0].plot(temperatures_keV, DT_simulated_reactivity, marker='d', c=good_colors[2], linestyle='none',
               ms=10, markeredgewidth=1, markeredgecolor='k', label='D-T')

    ax[1].plot(temperatures_keV, DD_fwhm, marker='o', c=good_colors[3], linestyle='none',
               ms=10, markeredgewidth=1, markeredgecolor='k', label='D-D')
    ax[1].plot(temperatures_keV, DT_fwhm, marker='d', c=good_colors[2], linestyle='none',
               ms=10, markeredgewidth=1, markeredgecolor='k', label='D-T')

    ###########################################
    # finalize plot

    ax[0].legend(loc=4)
    ax[0].text(0.04, 0.9, '(a)', fontsize=18, weight='bold', transform=ax[0].transAxes)
    ax[1].text(0.04, 0.9, '(b)', fontsize=18, weight='bold', transform=ax[1].transAxes)

    fig.tight_layout(pad=0.1, rect=[0, 0, 1, 1])
    plt.show()
    # if (SAVE_FIGS): plt.savefig(SAVE_IMG_DIR+'thermonuclear_moments.pdf')
