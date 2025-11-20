"""
Created by M Lavell on 02/16/2024

Simulate 50 MeV electrons slowing in bulk aluminum w/ and w/o radiation

Recreating Fig 5 from Wu 2018 atomic processes paper
"""

import numpy as np
import matplotlib.pyplot as plt
from pyforce import *
from scipy import constants
import os, sys, h5py

from collision_diagnostics import *

# ========== Program Settings =============
# =========================================
prog_set = ProgramSettings("binary_collisions_program")
prog_set.debug_verbosity = 3
prog_set.output_verbosity = 1
prog_set.n_threads = 1
prog_set.n_dims = 2
prog_set.gpu = False

# ========== Simulation Params ============
# =========================================
sim_param = SimulationParameters()
sim_param.total_time = 150.0e-12
sim_param.dt = 1.0e-12
sim_param.enable_cfl = False
sim_param.cfl = 1.0

# =============    Mesh   =================
# =========================================
x_range = [0.0, 1.0e-6]
z_range = [0.0, 1.0e-6]
x = np.linspace(x_range[0], x_range[1], 2)
z = np.linspace(z_range[0], z_range[1], 2)
gridlines = [x, z]
mesh_param = Mesh('simulation_mesh')
mesh_param.filename = "grid_1.csv"
mesh_param.grid_coordinates = gridlines

# ======= Materials and geometry ==========
# =========================================
mat = [Material("void", is_pec=True)]

internal_bcs1 = InternalBoundaryCondition("void")
internal_bcs1.geometry = "box"
internal_bcs1.material = "void"
internal_bcs1.particle_interaction = "transparent"
internal_bcs = [internal_bcs1]

geo = [Polygon("box", points=[
    [x_range[0], z_range[0]],
    [x_range[0], z_range[1]],
    [x_range[1], z_range[1]],
    [x_range[1], z_range[0]]
])]

# =========== IO parameters ===============
# =========================================
io_param = IO()
# io_param.output_file_path = "output"

io_param.checkpoint_enable = 1
io_param.checkpoint_file_path = "checkpoints"
io_param.checkpoint_frequency = Interval(value=1.0, units="wall_time")
io_param.checkpoint_total_number = 1

# ================ EM =====================
# =========================================
em_param = Electromagnetics()
em_param.solver_type = 'none'
em_param.boundary_conditions = {'x0': 'pml', 'x1': 'pml', 'y0': 'pml', 'y1': 'pml'}
em_param.pml_depth = 10
em_param.pml_damping_max = 0.35
em_param.initial_poisson = True
em_param.poisson_iterations = 100
em_param.poisson_tol = 1.0E-8
em_param.conformal_correction = 'corrected'

# =========== Particle pusher =============
# =========================================
particle_pusher = ParticlePusher()
particle_pusher.enable_particle_pusher = False
particle_pusher.push_max_iterations = 50

# ============ Particle Types =============
# =========================================

electron_type = ParticleType("electron_type")
electron_type.mass = constants.m_e
electron_type.charge = -1

aluminum_type = ParticleType("aluminum_type")
aluminum_type.mass = 63.55 * constants.atomic_mass
aluminum_type.charge = 13
aluminum_type.atomic_number = 13

photon_type = ParticleType("photon_type")
photon_type.mass = 0.0
photon_type.charge = 0

# =========== Particle Groups =============
# =========================================

electron_beam_energy = 5.0e6
electron_velocity = calc_beam_velocity_EeV(electron_beam_energy, electron_type.mass)

# define electron group
electron_group = ParticleGroup("electron")
electron_group.particle_type = "electron_type"
electron_group.max_count = int(1e5)
electron_group.num_subcycling_groups = 1
electron_group.apm_params = ApmParams(min_ppc_range=80, max_ppc_range=120, target_ppc=100,
                                      interval=10, activation_ppc=50, subcells=[5, 5])
electron_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)

electron_temperature = 0.0

# same geometry except with twice the temperature
electron_group_geom = ParticleGeometry("box_geometry")
electron_group_geom.geometry_name = "box"
electron_group_geom.ppc = [100, 100]
electron_group_geom.density = 7.76e29 #7.76e26
electron_group_geom.temperature = [electron_temperature, electron_temperature, electron_temperature]
electron_group_geom.velocity = [0.0, 0.0, electron_velocity]
electron_group_geom.relativistic_initialization = False

electron_group.init_params = [electron_group_geom]

# define boron group
aluminum_group = ParticleGroup("aluminum")
aluminum_group.particle_type = "aluminum_type"
aluminum_group.max_count = int(1e5)
aluminum_group.num_subcycling_groups = 1
aluminum_group.apm_params = ApmParams(min_ppc_range=80, max_ppc_range=120, target_ppc=100,
                                    interval=10, activation_ppc=50, subcells=[5, 5])
aluminum_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)

aluminum_temperature = 0.1

# same geometry except with twice the temperature
aluminum_group_geom = ParticleGeometry("box_geometry")
aluminum_group_geom.geometry_name = "box"
aluminum_group_geom.ppc = [100, 100]
aluminum_group_geom.density = 7.76e29
aluminum_group_geom.temperature = [aluminum_temperature, aluminum_temperature, aluminum_temperature]
aluminum_group_geom.velocity = [0.0, 0.0, 0.0]
aluminum_group_geom.relativistic_initialization = False

aluminum_group.init_params = [aluminum_group_geom]


# define photon group
photon_group = ParticleGroup("photon")
photon_group.particle_type = "photon_type"
photon_group.max_count = int(1e6)
photon_group.num_subcycling_groups = 1
photon_group.apm_params = ApmParams(min_ppc_range=80, max_ppc_range=120, target_ppc=100,
                                    interval=10, activation_ppc=50, subcells=[5, 5])
photon_group.chaining_mesh_params = ChainingMeshParams(update_interval=1, sort_interval=100)
photon_group.init_params = [ParticleNullInit("empty_initialization")]

# ========= Timers and Metrics ===========
# =========================================
timer1 = Timer("timer_wall_time", Interval(value=1.0, units="wall_time"))
metric1 = IntegratedMetric("integrated", Interval(value=1, units="step_count"))

# ============= Collisions ================
# =========================================

collision_ei = Collision(group1_name="electron", group2_name="aluminum", interaction_type="binary")
ei_Coulomb_channel = Coulomb(disable=False)
ei_Coulomb_channel.enable_radiation_correction = False
ei_Coulomb_channel.photon_group_id = 2
ei_Coulomb_channel.radiation_step_count_interval = 1
collision_ei.channels = [ei_Coulomb_channel]

# ================ Input ==================
# =========================================
global_input = Input()
global_input.program_settings = prog_set
global_input.simulation_parameters = sim_param
global_input.meshes = [mesh_param]
global_input.io_parameters = io_param
global_input.electromagnetics = em_param # void
global_input.materials = mat  # void
global_input.geometries = geo
global_input.internal_bcs = internal_bcs # void
global_input.particle_types = [electron_type, aluminum_type, photon_type]
global_input.particle_groups = [electron_group, aluminum_group, photon_group]
global_input.particle_pusher = particle_pusher # void
global_input.timers = [timer1]
global_input.metrics = [metric1]
global_input.collisions = [collision_ei]


# ========= Compile and Run ===============
# =========================================


# change DD output directories and temperatures
io_param.output_file_path = "output/aluminum_stopping"
if not os.path.isdir(io_param.output_file_path):
    os.mkdir(io_param.output_file_path)

global_input.compile('aluminum_stopping.toml', relocate=True, bin_name='tflink-collisions')
global_input.run('tflink-collisions', 'aluminum_stopping.toml')


# # ========= Post-processing ===============
# # =========================================
print("\nPost-processing...\n")

plt.style.use("triforce.mplstyle")
good_colors=["#db6d00","#006ddb","#920000","#52a736","#9B30FF"]

# Start plot
fig, ax = plt.subplots(2,1,figsize=(8,6))

ax[0].set_xlabel("t (ps)")
ax[0].set_ylabel(r"$E_k$ (MeV)")
ax[0].set_xlim([0, sim_param.total_time * 1e12])
# ax[0].set_ylim([0.0, 4.0e19])
# ax[0].set_yscale('log')

ax[1].set_xlabel(r"$E_k$ (MeV)")
ax[1].set_ylabel(r"$dN/dE_k$")
# ax[1].set_xlim([0.0, 55.0])
# ax[1].set_ylim([0.0, 3.0e18])
# ax[1].set_xscale('log')
# ax[1].set_yscale('log')

for axis in ax:
    axis.grid()

###########################################
# process output data

# Time Mass KE Tavg Tx Ty Tz px py pz vx vy vz vmag

electron_integrated_data = np.loadtxt("output/aluminum_stopping/electron_integrated.csv")
time = electron_integrated_data[:,0] * 1e12
v_z = electron_integrated_data[:,12]
E_beam1 = calc_beam_energy_v(v_z, electron_type.mass, abs(electron_type.charge)) * 1e-6

v_mag = electron_integrated_data[:,13]
E_beam2 = calc_beam_energy_v(v_mag, electron_type.mass, abs(electron_type.charge)) * 1e-6


ax[0].plot(time, E_beam1, '-k')
ax[0].plot(time, E_beam2, '-b')


me_c2 = constants.electron_mass * constants.c**2
n_bins = 100

# electron_data_t0 = h5py.File(f"output/aluminum_stopping/electron_particles_000000.hdf5", "r")
#
# n_electrons = electron_data_t0["electron"].attrs["size"]
# electron_energies = np.zeros(n_electrons)
# for i in range(n_electrons):
#     gamma = electron_data_t0["electron"]["gamma"][i]
#     electron_energies[i] = (gamma - 1.0) * me_c2 * 1e-6 / constants.e
#
# ax[1].hist(electron_energies, bins=n_bins)

electron_data_t150 = h5py.File(f"output/aluminum_stopping/electron_particles_000150.hdf5", "r")

n_electrons = electron_data_t150["electron"].attrs["size"]
electron_energies = np.zeros(n_electrons)
electron_energies2 = np.zeros(n_electrons)
for i in range(n_electrons):
    gamma = electron_data_t150["electron"]["gamma"][i]
    electron_energies[i] = (gamma - 1.0) * me_c2 * 1e-6 / constants.e

    momentum2 = sum(electron_data_t150["electron"]["momentum"][i]**2)
    electron_energies2[i] = momentum2 * 1e-6 / (constants.electron_mass * constants.e)
    # electron_energies[i] = 0.5 * constants.electron_mass * momentum[2]**2 * 1e-6 / (constants.e)

ax[1].hist(electron_energies, bins=n_bins)
# ax[1].hist(electron_energies2, bins=n_bins)



# def find_nearest(array,value):
#     array = np.asarray(array)
#     idx = (np.abs(array-value)).argmin()
#     return array[idx], idx




# photon diagnostics (KE,p) are normalized by me*c^2
# me_c2 = constants.electron_mass * constants.speed_of_light**2
#
# photon_integrated_data = np.loadtxt("output/cuBrems_0ppc/photon_integrated.csv")
# time_sim = photon_integrated_data[:,0] * 1e15
# photon_weight = photon_integrated_data[:,1]
#
# # total photons generated over time
# ax[0].plot(time_sim, photon_weight, 'd', c=good_colors[1], mec="k", mew=1.5, ms=6)

# # next get photon energy spectrum at last step
# n_kvalues = 100
# k_over_gamma_e_bins = np.zeros(n_kvalues)
# k_over_gamma_e_bin_centers = np.logspace(-7,0, n_kvalues)
#
# print("getting photon energy spectrum")
#
# # photon_data = h5py.File(f"output/cuBrems_0ppc/photon_particles_{i:06}.hdf5", "r")
# photon_data = h5py.File(f"output/cuBrems_0ppc/photon_particles_000042.hdf5", "r")
#
# n_photons = photon_data["photon"].attrs["size"]
# photon_energies = np.zeros(n_photons)
#
# def find_nearest(array,value):
#     array = np.asarray(array)
#     idx = (np.abs(array-value)).argmin()
#     return array[idx], idx
#
# for j in range(n_photons):
#
#     # save photon momentum as p_save = p_gamma / (me * c^2)
#     # p_gamma = E * c = omega * hbar * c
#     # E_photon_normalized = omega * hbar / (me * c^2) = p_save  * c
#
#     # get size of momentum vector
#     p_over_mec2 = sum(photon_data["photon"]["momentum"][j]**2)**0.5
#     E_photon_normalized = p_over_mec2 * constants.c
#     E_over_gamma_m1 = E_photon_normalized / (gamma_e - 1.0)
#
#     photon_energies[j] = E_over_gamma_m1
#
#     print(f"Enorm={E_over_gamma_m1}   gamma_e={gamma_e}  gm1={gamma_e-1}")
#
#     # find nearest bin and accumulate weight
#     closest_val, bin_id = find_nearest(k_over_gamma_e_bin_centers, E_over_gamma_m1)
#     k_over_gamma_e_bins[bin_id] += photon_data["photon"]["weight"][j]
#
#     # print(f"E_over_gamma_m1={E_over_gamma_m1}")
#
# ax[1].plot(k_over_gamma_e_bin_centers, k_over_gamma_e_bins, 'o', c=good_colors[2], mec="k", mew=1.5, ms=3)

# ax.legend()
fig.tight_layout(pad=0.5, rect=[0,0,1,1])
plt.show()