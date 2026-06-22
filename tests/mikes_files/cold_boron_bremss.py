#!/usr/bin/env python3
"""
Created by M Lavell on 05/06/24

Simulate Bremsstrahlung photons from 1 MeV electrons incident on boron.

Similar to test from Martinez et at. 2019, "High-energy radiation and pair
production by Coulomb processes in PIC simulations"
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt
from pyforce import *
from scipy import constants
import os, sys, h5py

from tfd_functions import *

home = os.getenv("HOME")
project_dir = home + "/TriForce/proton_boron/"
cross_section_dir = project_dir + "simulations/cross_sections/"

plt.style.use(project_dir + "extras/triforce.mplstyle")
good_colors = ["#db6d00", "#006ddb", "#920000", "#52a736", "#9B30FF"]

# ========== Program Settings =============
# =========================================
prog_set = ProgramSettings(
    'boron_bremss',
    dims=2,
    debug_verbosity=3,
    output_verbosity=1,
    print_step=10,
    num_threads=1,
    dy=1.0E-6
)

# ========== Simulation Params ============
# =========================================
sim_param = SimulationParameters(
    total_time=40.0E-15,
    dt=1.0E-15,
    cfl=1.0,
    enable_cfl=False
)

# =============    Mesh   =================
# =========================================
x_range = [0.0, 1.0e-6]
z_range = [0.0, 1.0e-6]
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

boron_charge = 1
boron_temperature = 1.0  #30.0e3  # eV
boron_atomic_weight = 10.811
boron_density = 2.34  # g/cc
boron_ndensity = boron_density * constants.Avogadro / boron_atomic_weight * 1.0e6  # m^-3

electron_beam_energy = 1.0e6
electron_velocity = calc_beam_velocity_EeV(electron_beam_energy, constants.m_e)
electron_temperature = 0.0
electron_ndensity = boron_ndensity / 10.0  # * boron_charge

# define electron group
electron_group = ParticleGroup("electron")
electron_group.mass = constants.m_e
electron_group.charge = -1
electron_group.max_count = int(1e5)
electron_group.num_subcycling_groups = 1
electron_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, "step_count"), disable=True)
electron_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)

# same geometry except with twice the temperature
electron_group_geom = ParticleGeometry(
    'box_geometry',
    thebox,
    ppc=[20, 20],
    density=electron_ndensity,
    temperature=[electron_temperature, electron_temperature, electron_temperature],
    velocity=[0.0, 0.0, electron_velocity]
)
electron_group.init_params = [electron_group_geom]

# define boron group
boron_group = ParticleGroup("boron")
boron_group.mass = boron_atomic_weight * constants.atomic_mass
boron_group.charge = boron_charge
boron_group.atomic_number = 5
boron_group.max_count = int(1e5)
boron_group.num_subcycling_groups = 1
boron_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, "step_count"), disable=True)
boron_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)

boron_group_geom = ParticleGeometry(
    'box_geometry',
    thebox,
    ppc=[20, 20],
    density=boron_ndensity,
    temperature=[boron_temperature, boron_temperature, boron_temperature],
    velocity=[0.0, 0.0, 0.0],
    relativistic_init=True
)
boron_group.init_params = [boron_group_geom]

# define photon group
photon_group = ParticleGroup("photon")
photon_group.mass = 0.0
photon_group.charge = 0
photon_group.max_count = int(1e6)
photon_group.num_subcycling_groups = 1
photon_group.apm_params = ApmParams(target_ppc=[10, 10], interval=Interval(10, "step_count"), disable=True)
photon_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
photon_group.init_params = [ParticleNullInit("empty_initialization")]

# ========= Timers and Metrics ===========
# =========================================
metric_interval = 1
integrated_metric = IntegratedMetric("integrated", Interval(value=metric_interval, units="step_count"))
collision_metric = CollisionMetric("collisions", Interval(value=metric_interval, units="step_count"))
photon_metric = PhotonSpectrum("photon_spectrum", Interval(value=metric_interval, units="step_count"),
                               group_id=2, n_bins=500)
particle_metric = ParticleData("particle_data", Interval(value=metric_interval, units="step_count"),
                               group_ids=[0, 1, 2], sample_fraction=1.0)

# ============= Collisions ================
# =========================================

collision_ei = Collision(group1_name="electron", group2_name="boron", interaction_type="binary")
collision_ei.channels = [Coulomb(disable=True),
                         BremsstrahlungTFD(photon_group_id=2,
                                           production_multiplier=1.0e30,
                                           n_photon_energies=200,
                                           n_electron_energies=1,
                                           min_electron_energy=1.0e6,
                                           max_electron_energy=1.0e6,
                                           #enable_Debye_model=True,
                                           #update_at_delta_Debye=0.1,
                                           reduce_electron_energy=False,
                                           disable=False)]

# ================ Input ==================
# =========================================
global_input = Input(
    program_settings=prog_set,
    simulation_params=sim_param,
    meshes=mesh_param,
    io_params=io_param,
    electromagnetics=em_param,
    particle_groups=[electron_group, boron_group, photon_group],
    particle_pusher=particle_pusher,
    collisions=[collision_ei],
    metrics=[integrated_metric, collision_metric, photon_metric, particle_metric]
)

# ========= Compile and Run ===============
# =========================================
run_sim = True

# ni = 10*ne
# wi = 2*we, Ni = 20*Ne
# wi = we, Ni = 10*Ne
# 2*wi = we, Ni = 5*Ne

# sim_names = ["wi2we", "wiwe", "2wiwe"]
# electron_group_geom.ppc = [10, 100]
# ion_ppcs = [[200, 100], [100, 100], [50, 100]]
# labels = [r"$w_{\mathrm{i}} = 2w_{\mathrm{e}}$",
#           r"$w_{\mathrm{i}} = w_{\mathrm{e}}$",
#           r"$2w_{\mathrm{i}} = w_{\mathrm{e}}$"]

sim_names = ["wi5we", "wiwe", "5wiwe"]
electron_group_geom.ppc = [10, 100]
ion_ppcs = [[20, 100], [100, 100], [200, 100]]
labels = [r"$w_{\mathrm{i}} = 5w_{\mathrm{e}}$",
          r"$w_{\mathrm{i}} = w_{\mathrm{e}}$",
          r"$5w_{\mathrm{i}} = w_{\mathrm{e}}$"]

markers = ["o", "s", "d"]

n_sims = len(sim_names)

for i in range(n_sims):

    if not os.path.isdir("output"):
        os.mkdir("output")

    io_param.output_file_path = "output/boron_brems_" + sim_names[i]
    if not os.path.isdir(io_param.output_file_path):
        os.mkdir(io_param.output_file_path)

    # change ion particle count
    boron_group_geom.ppc = ion_ppcs[i]

    if run_sim:
        global_input.compile('boron_brems.toml', relocate=True, bin_name='tflink')
        global_input.run('tflink', 'boron_brems.toml')

# # ========= Post-processing ===============
# # =========================================
run_postproc = True
if run_postproc:

    print("\nPost-processing...\n")

    # Start plot
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))

    ax[0].set_xlabel("Time (fs)")
    ax[0].set_ylabel(r"$N_{\gamma}$")
    ax[0].set_xlim([0, sim_param.total_time * 1e15])
    ax[0].set_ylim([0.0, 1.1e19])
    # ax[0].set_yscale('log')

    ax[1].set_xlabel(r"$k/(\gamma_e -1)$")
    ax[1].set_ylabel(r"$k\mathrm{d}N_{\gamma}/\mathrm{d}k$")
    ax[1].set_xlim([1.0e-7, 1.0])
    ax[1].set_ylim([0.0, 8.0e17])
    ax[1].set_xscale('log')
    # ax[1].set_yscale('log')

    for axis in ax:
        axis.grid()
        axis.tick_params(axis='x', pad=8)

    ax[0].set_title(r"(a)", fontsize=16, loc='left', x=0.05, y=0.85)
    ax[1].set_title(r"(b)", fontsize=16, loc='right', x=0.95, y=0.85)

    ###########################################
    # theory solutions (eq 38 and 39 from Martinez 2019)

    gamma = 1.0 + (electron_beam_energy * constants.e) / me_c2
    gm1 = gamma - 1.0

    nks = 200
    k_over_gm1 = np.logspace(m.log10(1.0e-7), m.log10(1.0 - 1.0e-7), nks)
    ks = k_over_gm1 * gm1

    dsdk = np.zeros(nks)
    for i in range(nks):
        # dsdk[i] = f_diff_cs(ks[i], boron_group.atomic_number, boron_group.charge,
        #                     boron_temperature, boron_group_geom.density, gamma)
        dsdk[i] = f_diff_cs(ks[i], boron_group.atomic_number, 0.0, 0.0,
                            0.0, boron_group_geom.density, electron_group_geom.density, gamma)
    kdsdk = ks * dsdk

    # sigma_ttl = ttl_cs_gauss(boron_group.atomic_number, boron_group.charge,
    #                          boron_temperature, boron_group_geom.density, gamma)
    sigma_ttl = ttl_cs_gauss(boron_group.atomic_number, 0.0, 0.0,
                             0.0, boron_group_geom.density, electron_group_geom.density, gamma)
    print(f"sigma_ttl_py = {sigma_ttl}")
    # sigma_ttl = 5.59175e-28

    N_theory_pts = 100
    time_theory = np.linspace(0, sim_param.total_time, N_theory_pts)

    gamma_e = (1.0 - (electron_velocity / constants.c) ** 2) ** (-0.5)
    gm1 = gamma_e - 1.0
    length = z_range[1]
    Ngamma_theory = electron_group_geom.density * boron_group_geom.density * electron_velocity * \
                    length * sigma_ttl * time_theory
    ax[0].plot(time_theory * 1e15, Ngamma_theory, '-k', label="Theory")

    dNdk = electron_group_geom.density * boron_group_geom.density * electron_velocity * \
           length * time_theory[-1] * kdsdk
    ax[1].plot(k_over_gm1, dNdk, '-k')

    ###########################################
    # process output data

    print("getting photon count")

    # photon diagnostics (KE,p) are normalized by me*c^2
    me_c2 = constants.electron_mass * constants.speed_of_light ** 2
    sim_volume = (x_range[1] - x_range[0]) * (z_range[1] - z_range[0]) * prog_set.dy
    area = (z_range[1] - z_range[0]) * prog_set.dy  # since paper simulation is 1D

    for i in range(n_sims):
        print(f"processing sim: {i} {sim_names[i]}")

        photon_integrated_data = np.loadtxt(f"output/boron_brems_{sim_names[i]}/photon_integrated.csv")
        time_sim = photon_integrated_data[:, 0] * 1.0e15
        photon_weight = photon_integrated_data[:, 1]

        photon_density = photon_weight / area

        # total photons generated over time
        ax[0].plot(time_sim[::2], photon_density[::2], c=good_colors[i], mec=good_colors[i],
                   mew=1.5, ms=5, fillstyle="none", linestyle="none", marker=markers[i], label=labels[i])

        # next get photon energy spectrum at last step
        n_bins = 50
        n_cnt = np.zeros(n_bins)
        dndk = np.zeros(n_bins - 1)
        k_over_gm1_bins = np.logspace(-7, 0, n_bins)

        i_step = 40
        photon_data = h5py.File(f"output/boron_brems_{sim_names[i]}/photon_particles_{i_step:08}.h5", "r")

        n_photons = photon_data["photon"].attrs["size"]
        photon_energies = np.zeros(n_photons)

        print(f"getting photon energy spectrum for n_photons = {n_photons}")

        for j in range(n_photons):
            k_over_gm1 = photon_data["photon"]["gamma"][j] / (gamma_e - 1.0)

            i_bin = 0
            while i_bin < n_bins - 1 and k_over_gm1 > k_over_gm1_bins[i_bin]:
                i_bin += 1
                if i_bin == 100:
                    print(j, k_over_gm1, gamma_e - 1.0)

            n_cnt[i_bin] += photon_data["photon"]["weight"][j] / area

        for j in range(n_bins - 1):
            kj = k_over_gm1_bins[j] * gm1
            kjp1 = k_over_gm1_bins[j + 1] * gm1
            k = 0.5 * (kjp1 + kj)
            dk = kjp1 - kj

            dndk[j] = n_cnt[j] / dk * k

        ax[1].plot(k_over_gm1_bins[2:], dndk[1:], c=good_colors[i], mec=good_colors[i],
                   mew=1.5, ms=5, fillstyle="none", linestyle="none", marker=markers[i], label=labels[i])

    ax[0].legend(loc=4)
    fig.tight_layout(pad=0.2, rect=[0, 0, 1, 1])
    # plt.savefig("solid_cold_boron_bremss.pdf")
    plt.show()
