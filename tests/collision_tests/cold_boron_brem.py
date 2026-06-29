#!/usr/bin/env python3
import numpy as np

from scripts.pyforce import *

# =============================
# ===== Simulation Params =====
# =============================
sim_name = 'boron_brems'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
data_path = project_path + f'/data/{sim_name}'

shape = (2, 2, 2)

xmin, xmax = 0.0, 1.6e-6
ymin, ymax = 0.0, 1.0e-6
zmin, zmax = 0.0, 1.6e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 1.0e-15
t_end = 40.0e-15
nt = int(t_end / dt) + 1

# =====================
# ===== Particles =====
# =====================
px_range=(xmin, xmax)
py_range=(ymin, ymax)
pz_range=(zmin, zmax)

E_eV = 1.0e6
v_e = constants.c * np.sqrt(1 - (1 + constants.eV * E_eV / (constants.m_e * constants.c**2))**-2)
gamma_e = (1.0 - (v_e / constants.c)**2)**(-0.5)
B_den = 1.304e29 # m^-3, 2.34 g/cm^3
e_den = B_den / 10.0

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=(0.0, 0.0, E_eV), # eV
    density=e_den, # m^-3,
    ppc=(10, 10, 100),
    distribution='constant',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

boron = Particles(
    name='boron',
    mass=10.811 * constants.atomic_mass,
    charge=1,
    atomic_number=5,
    temp=(1.0, 1.0, 1.0), # eV
    density=1.304e29, # m^-3, 2.34 g/cm^3
    ppc=(20, 20, 100),
    distribution='relativistic',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

photons = Particles(
    name='photons',
    mass=0.0,
    charge=0.0,
    atomic_number=0,
    temp=(0.0, 0.0, 0.0), # eV
    density=0.0, # m^-3,
    ppc=(0, 0, 0),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ==========================================
# ===== Collisions and Particle Params =====
# ==========================================
radiation_params = RadiationParams(
    products=photons,
    reduce_electron_energy=False,
    production_multiplier=1.0e30,
    cross_section_file=project_path + '/tests/cross_section_data/SB_G4_Z5_kdsdk_MeV_barns.csv'
)

particle_params = ParticleParams(
    save_interval=1,
    particle_bcs=ParticleBCType.Periodic,
    interp_order=1,
    particle_data=(electrons, boron, photons),
    collisions=(
        Collision(
            groups=(electrons, boron),
            channels=('radiation',),
            radiation=radiation_params,
        ),
    )
)

# ==========================
# ===== Metrics Params =====
# ==========================
metric_params = Metrics(
    data_path,
    (MetricType.ParticleDump,)
)

# ============================
# ===== Simulation Class =====
# ============================
sim_params = Simulation(
    name=sim_name,
    shape=shape,
    nthreads=4,
    dt=dt,
    t_end=t_end,
    nt=nt,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    particle_params=particle_params,
    metric_params=metric_params,
    em_enabled=False,
    push_enabled=False,
    jdep_enabled=False,
    velocity_backstep_enabled=False,
    collisions_enabled=True
)

# # ===========================
# # ===== Compile and Run =====
# # ===========================
# print(f'Setting up "{sim_name}"')
# create_data_dir(data_path)
# create_particles(sim_params, [electrons, boron, photons], data_path)
# update_header(sim_params, project_path=project_path, data_path=data_path)
#
# compile_project(build_path, output=True)
# run_project(build_path + '/game_engine', output=True)

# ===========================
# ===== Post Processing =====
# ===========================
fig, ax = plt.subplots(1, 2, figsize=(9, 4))

gamma = 1.0 + (E_eV * constants.e) / me_c2
gm1 = gamma - 1.0

nks = 200
k_over_gm1 = np.logspace(np.log10(1.0e-7), np.log10(1.0 - 1.0e-7), nks)
ks = k_over_gm1 * gm1

dsdk = np.zeros(nks)
for i in range(nks):
    dsdk[i] = f_diff_cs(ks[i], 5, 0, 0, 0, B_den, e_den, gamma)

kdsdk = ks * dsdk

sigma_ttl = ttl_cs_gauss(5, 0, 0, 0, B_den, e_den, gamma)

N_theory_pts = 100
time_theory = np.linspace(0, t_end, N_theory_pts)

gamma_e = (1.0 - (v_e / constants.c) ** 2) ** (-0.5)
gm1 = gamma_e - 1.0
length = zmax
Ngamma_theory = e_den * B_den * v_e * length * sigma_ttl * time_theory
ax[0].plot(time_theory * 1e15, Ngamma_theory, '-k', label="Theory")

dNdk_thry = e_den * B_den * v_e * length * time_theory[-1] * kdsdk
ax[1].plot(k_over_gm1, dNdk_thry, '-k')
ax[1].set_xscale('log')

area = dx * dy
photon_weights = np.zeros(nt)
for i in range(1, nt):
    file_name = f'/data/{sim_name}/photons_dump_{i:010d}.bp'
    with FileReader(project_path + file_name) as f:
        photon_weights[i] = f.read('Weight').sum()

photon_density = photon_weights / area

ax[0].plot(photon_density)


i_step = 40
file_name = f'/data/{sim_name}/photons_dump_{i_step:010d}.bp'
with FileReader(project_path + file_name) as f:
    weights = f.read('Weight').flatten()
    gammas = f.read('Gamma').flatten()

n_bins = 100
k_over_gm1 = gammas / (gamma_e - 1.0)
k_over_gm1_bins = np.logspace(-7, 0, n_bins)

n_count, bins = np.histogram(k_over_gm1, bins=k_over_gm1_bins, weights=weights / area)

kgm1 = k_over_gm1_bins * (gamma_e - 1.0)
ks = 0.5 * (kgm1[1:] + kgm1[:-1])
dk = kgm1[1:] - kgm1[:-1]
dNdk = ks * n_count / dk

ax[1].scatter(k_over_gm1_bins[1:], dNdk)

plt.show()











