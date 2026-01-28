#!/usr/bin/env python3
"""
Created by M Lavell on 02/13/24

Simulate Bremsstrahlung photons from 40 MeV electrons incident
on neutral copper atoms.

Test from Martinez et at. 2019, "High-energy radiation and pair
production by Coulomb processes in PIC simulations"

NOTE TEST PROBLEM IS WORK IN PROGRESS - need to revisit with new radiation model

"""

import numpy as np
import matplotlib.pyplot as plt
from pyforce import *
from scipy import constants
import os, sys, h5py

from collision_diagnostics import *

home = os.getenv("HOME")

# ========== Program Settings =============
# =========================================
prog_set = ProgramSettings(
    'copper_bremss',
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
    total_time=42.0E-15,
    dt=1.0E-15,
    cfl=1.0,
    enable_cfl=False
)

# =============    Mesh   =================
# =========================================
x_range = [0.0, 1.6e-6]
z_range = [0.0, 1.6e-6]
x = np.linspace(x_range[0], x_range[1], 2)
z = np.linspace(z_range[0], z_range[1], 2)

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
io_param.output_file_path = "output"
io_param.timer_write_interval = Interval(value=1000, units='step_count')

io_param.checkpoint_enable = 1
io_param.checkpoint_file_path = "checkpoints"
io_param.checkpoint_frequency = Interval(value=1.0, units="wall_time")
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

electron_beam_energy = 40.0e6 # 40 MeV
electron_velocity = calc_beam_velocity_EeV(electron_beam_energy, constants.m_e)

# define electron group
electron_group = ParticleGroup("electron")
electron_group.mass = constants.m_e
electron_group.charge = -1
electron_group.max_count = int(1e5)
electron_group.num_subcycling_groups = 1
electron_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, "step_count"), disable=True)
electron_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)

electron_temperature = 0.0

# same geometry except with twice the temperature
electron_group_geom = ParticleGeometry(
    'box_geometry',
    thebox,
    ppc=[100, 10],
    density=1.0e27,
    temperature=[electron_temperature, electron_temperature, electron_temperature],
    velocity=[0.0, 0.0, electron_velocity]
)
electron_group.init_params = [electron_group_geom]

# define boron group
copper_group = ParticleGroup("copper")
copper_group.mass = 63.55 * constants.atomic_mass
copper_group.charge = 0
copper_group.atomic_number = 29
copper_group.max_count = int(1e5)
copper_group.num_subcycling_groups = 1
copper_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, "step_count"), disable=True)
copper_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)

copper_temperature = 1.0

# same geometry except with twice the temperature
copper_group_geom = ParticleGeometry(
    'box_geometry',
    thebox,
    ppc=[100, 10],
    density=8.0e28,
    temperature=[copper_temperature, copper_temperature, copper_temperature],
    velocity=[0.0, 0.0, 0.0]
)
copper_group.init_params = [copper_group_geom]


# define photon group
photon_group = ParticleGroup("photon")
photon_group.mass = 0.0
photon_group.charge = 0
photon_group.max_count = int(10e6)
photon_group.num_subcycling_groups = 1
photon_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, "step_count"), disable=True)
photon_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
photon_group.init_params = [ParticleNullInit("empty_initialization")]

# ============= Collisions ================
# =========================================

cs_dir = home + "/TriForce/triforce/libs/common/data/cross_sections/"

collision_ei = Collision(group1_name="electron", group2_name="copper", interaction_type="binary")
collision_ei.channels = [Coulomb(disable=True),
                         Radiation(photon_group_id=2,
                                   production_multiplier=1.0e5,
                                   cross_section_filename=cs_dir+"SB_G4_Z29_kdsdk_barns.csv",
                                   scale_table_to_eV=1.0e6,
                                   scale_table_to_barns=1.0,
                                   reduce_electron_energy=False,
                                   disable=False)]

# ========= Timers and Metrics ===========
# =========================================
metric_interval = 1
integrated_metric = IntegratedMetric("integrated", Interval(value=metric_interval, units="step_count"))
collision_metric = CollisionMetric("collisions", Interval(value=metric_interval, units="step_count"))
photon_metric = PhotonSpectrum("photon_spectrum", Interval(value=metric_interval, units="step_count"),
                               group_id=2, n_bins=500)
particle_metric = ParticleData("particle_data", Interval(value=metric_interval, units="step_count"),
                               group_ids=[0,1,2], sample_fraction=1.0)


# ================ Input ==================
# =========================================
global_input = Input(
    program_settings=prog_set,
    simulation_params=sim_param,
    meshes=mesh_param,
    io_params=io_param,
    electromagnetics=em_param,
    particle_groups=[electron_group, copper_group, photon_group],
    particle_pusher=particle_pusher,
    collisions=[collision_ei],
    metrics=[integrated_metric, collision_metric, photon_metric, particle_metric]
)

# ========= Compile and Run ===============
# =========================================
run_sim = True

sim_names = ["wi2we","wiwe","2wiwe"]
electron_group_geom.ppc = [10, 10]
copper_ppcs = [[40,100], [80,100], [160,100]]
labels = [r"$w_{\mathrm{i}} = 2w_{\mathrm{e}}$",
          r"$w_{\mathrm{i}} = w_{\mathrm{e}}$",
          r"$2w_{\mathrm{i}} = w_{\mathrm{e}}$"]
markers = ["o", "s", "d"]

n_sims = len(sim_names)

for i in range(n_sims):

    if not os.path.isdir("output"):
        os.mkdir("output")

    io_param.output_file_path = "output/cuBrems_" + sim_names[i]
    if not os.path.isdir(io_param.output_file_path):
        os.mkdir(io_param.output_file_path)

    # change ion particle count
    copper_group_geom.ppc = copper_ppcs[i]

    if run_sim:
        global_input.compile('cupper_Brems.toml', relocate=True, bin_name='tflink-collisions')
        global_input.run('tflink-collisions', 'cupper_Brems.toml')
        # global_input.compile('cu_brems.toml', relocate=True, bin_name='tflink-1d-burn-wave')
        # global_input.run('tflink-1d-burn-wave', 'cu_brems.toml')



# # ========= Post-processing ===============
# # =========================================
run_postproc = True
if run_postproc:

    print("\nPost-processing...\n")

    plt.style.use("triforce.mplstyle")
    good_colors=["#db6d00","#006ddb","#920000","#52a736","#9B30FF"]

    # Start plot
    fig, ax = plt.subplots(1,2,figsize=(8,4))

    ax[0].set_xlabel("t (fs)")
    ax[0].set_ylabel(r"$N_{\gamma}$")
    ax[0].set_xlim([0, sim_param.total_time * 1e15])
    ax[0].set_ylim([0.0, 4.0e19])
    # ax[0].set_yscale('log')

    ax[1].set_xlabel(r"$k/(\gamma_e -1)$")
    ax[1].set_ylabel(r"$dN_{\gamma}/dk$")
    ax[1].set_xlim([1.0e-7, 1.0])
    ax[1].set_ylim([0.0, 3.0e18])
    ax[1].set_xscale('log')
    # ax[1].set_yscale('log')

    for axis in ax:
        axis.grid()

    ###########################################
    # theory solutions (eq 38 and 39 from Martinez 2019)

    N_theory_pts = 100
    # Ngamma_theory = np.zeros(N_theory_pts)
    # dNdk_theory = np.zeros(N_theory_pts)
    time_theory = np.linspace(0, sim_param.total_time, N_theory_pts)

    gamma_e = (1.0 - (electron_velocity / constants.c)**2)**(-0.5)
    sigmaB_mbarns = 3.54796e+01 # approximate total cross section at 40 MeV (from SB table)
    sigmaB = sigmaB_mbarns * 1.0e-28
    length = z_range[1]

    Ngamma_theory = electron_group_geom.density * copper_group_geom.density * electron_velocity \
                    * length * sigmaB * time_theory

    # ax[0].plot(time_theory * 1e15, Ngamma_theory, '-k')

    # Martinez_fig9b_data = np.loadtxt("temp_data/Martinez_fig9b_data.csv", delimiter=",")
    # ax[1].plot(Martinez_fig9b_data[:,0], Martinez_fig9b_data[:,1], '-k')


    ###########################################
    # process output data

    print("getting photon count")

    # photon diagnostics (KE,p) are normalized by me*c^2
    me_c2 = constants.electron_mass * constants.speed_of_light**2
    sim_volume = (x_range[1] - z_range[0]) * (z_range[1] - z_range[0]) * 1.0e-6
    area = (1.0e-6 * 1.6e-6) # since paper simulation is 1D

    for i in range(n_sims):
        print(f"processing sim: {i} {sim_names[i]}")

        photon_integrated_data = np.loadtxt(f"output/cuBrems_{sim_names[i]}/photon_integrated.csv")
        time_sim = photon_integrated_data[:,0] * 1.0e15
        photon_weight = photon_integrated_data[:,1]

        photon_density = photon_weight / area

        # total photons generated over time
        ax[0].plot(time_sim, photon_density, c=good_colors[i], mec=good_colors[i],
                   mew=1.5, ms=3.5, fillstyle="none", linestyle="none", marker=markers[i], label=labels[i])

        # next get photon energy spectrum at last step
        n_bins = 100
        n_cnt = np.zeros(n_bins)
        dndk = np.zeros(n_bins-1)
        k_over_gm1_bins = np.logspace(-7, 0, n_bins)

        i_step = 42
        photon_data = h5py.File(f"output/cuBrems_{sim_names[i]}/photon_particles_{i_step:08}.h5", "r")

        n_photons = photon_data["photon"].attrs["size"]
        photon_energies = np.zeros(n_photons)


        print(f"getting photon energy spectrum for n_photons = {n_photons}")

        for j in range(n_photons):
            k_over_gm1 = photon_data["photon"]["gamma"][j] / (gamma_e - 1.0)

            i_bin = 0
            while (i_bin < n_bins-1 and k_over_gm1 > k_over_gm1_bins[i_bin]):
                i_bin += 1
                if (i_bin == 100):
                    print(j,k_over_gm1,gamma_e-1.0)

            n_cnt[i_bin] += photon_data["photon"]["weight"][j] / area


        dk = np.zeros(n_bins-1)
        for j in range(n_bins-1):
            dk[j] = (k_over_gm1_bins[j+1] - k_over_gm1_bins[j]) * (gamma_e - 1.0)
            dndk[j] = n_cnt[j] / dk[j]

        ax[1].plot(k_over_gm1_bins[1:], dndk, c=good_colors[i], mec=good_colors[i],
                   mew=1.5, ms=3.5, fillstyle="none", linestyle="none", marker=markers[i], label=labels[i])

    ax[0].legend()
    fig.tight_layout(pad=0.5, rect=[0,0,1,1])
    plt.show()

###########################################################################################
# plot spectrum saved by code
plot_spectrum = False
if plot_spectrum:

    i_step = 42
    # filename = home + f"TriForce/triforce/examples/collisions/output/photon_spectrum_{i_step:08}.csv"
    filename = io_param.output_file_path + f"/photon_spectrum_{i_step:08}.csv"


    fig, ax = plt.subplots(1,1, figsize=(8,4))

    ax.set_ylabel("Intensity [arb.]")
    ax.set_xlabel(r"$\omega$")
    ax.grid()

    data = np.genfromtxt(filename, delimiter=',', dtype=float)
    bins = data[:,0]
    counts = data[:,1]
    cdf = data[:,2]

    ax.bar(bins, counts, width=(np.max(bins)-np.min(bins))/bins.shape[0])

    fig.tight_layout(pad=0.5, rect=[0,0,1,1])
    plt.show()
