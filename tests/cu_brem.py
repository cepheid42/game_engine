#!/usr/bin/env python3

import numpy as np
import subprocess
import matplotlib.pyplot as plt
from scipy import constants
from adios2 import FileReader

from scripts.domain_params import *
from scripts.particle_generation import create_particles

sim_name = 'cu_brems'
project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
particle_data = project_path + '/data/'

shape = (2, 2, 2)

xmin, xmax = 0.0, 1.6e-6
ymin, ymax = 0.0, 1.0e-6
zmin, zmax = 0.0, 1.6e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 1.0e-15
t_end = 4.2e-14
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

# ===== Particles =====
px_range=(xmin, xmax)
py_range=(ymin, ymax)
pz_range=(zmin, zmax)
e_eV = 40.0e6
v_e = constants.c * np.sqrt(1 - (1 + constants.eV * e_eV / (constants.m_e * constants.c**2))**-2)
e_den = 1.0e27 # m^3
cu_den = 8.0e28 # m^3

electrons = Particles(
    name='electrons',
    mass=constants.m_e,
    charge=-1,
    atomic_number=0,
    temp=(0.0, 0.0, e_eV), # eV
    density=e_den, # m^-3,
    ppc=(10, 1, 10),
    distribution='constant',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

cu_temp = 1.0e-50
copper = Particles(
    name='copper',
    mass=63.55 * constants.atomic_mass,
    charge=0,
    atomic_number=29,
    temp=(cu_temp, cu_temp, cu_temp), # eV
    # temp=(0.0, 0.0, 0.0), # eV
    density=cu_den, # m^-3,
    ppc=(100, 10, 10),
    distribution='thermal',
    # distribution='constant',
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

# ===== Collisions and Particle Params =====
radiation_params = RadiationParams(
    products='photons',
    reduce_electron_energy=False,
    production_multiplier=1.0e5,
    cross_section_file='/tests/cross_section_data/SB_G4_Z29_kdsdk_MeV_barns.csv'
)

particle_params = ParticleParams(
    save_interval=1,
    particle_bcs='periodic',
    interp_order=1,
    particle_data=(electrons, copper, photons),
    collisions=(
        Collision(
            groups=('electrons', 'copper'),
            channels=('radiation',),
            radiation=radiation_params,
        ),
    )
)

sim_params = Simulation(
    name=sim_name,
    shape=shape,
    nthreads=4,
    dt=dt,
    t_end=t_end,
    nt=nt,
    cfl=cfl,
    x_range=(xmin, xmax),
    y_range=(ymin, ymax),
    z_range=(zmin, zmax),
    deltas=(dx, dy, dz),
    particle_params=particle_params,
    em_enabled=False,
    push_enabled=False,
    jdep_enabled=False
)

sims = [
    ["wi2we", ( 40, 1, 100), r"$w_{\mathrm{i}} = 2w_{\mathrm{e}}$", 'o', "#db6d00"],
    # ["wiwe",  ( 80, 1, 100), r"$w_{\mathrm{i}} = w_{\mathrm{e}}$", 's', "#006ddb"],
    # ["2wiwe", (160, 1, 100), r"$2w_{\mathrm{i}} = w_{\mathrm{e}}$", 'd', "#920000"]
]


for name, cppc, _, _, _ in sims:
    sim_params.name = name
    copper.ppc = cppc
    print(f'Setting up "{name}"')

    create_particles(sim_params, electrons, particle_data)
    create_particles(sim_params, copper, particle_data)
    update_header(sim_params, project_path=project_path)

    subprocess.run(
        ['meson', 'compile', '-C', build_path, '-j4'],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    ).check_returncode()

    subprocess.run(build_path + '/game_engine').check_returncode()


plt.style.use('/home/cepheid/TriForce/game_engine/data/triforce.mplstyle')
good_colors=["#db6d00","#006ddb","#920000","#52a736","#9B30FF"]

# Start plot
fig, ax = plt.subplots(1,2,figsize=(8,4))

ax[0].set_xlabel("t (fs)")
ax[0].set_ylabel(r"$N_{\gamma}$")
ax[0].set_xlim([0, t_end * 1e15])
# ax[0].set_ylim([0.0, 1.0e19])
# ax[0].set_yscale('log')

ax[1].set_xlabel(r"$k/(\gamma_e -1)$")
ax[1].set_ylabel(r"$dN_{\gamma}/dk$")
# ax[1].set_xlim([1.0e-7, 1.0])
# ax[1].set_ylim([0.0, 3.0e18])
ax[1].set_xscale('log')
ax[1].set_yscale('log')

for axis in ax:
    axis.grid()

N_theory_pts = 100
time_theory = np.linspace(0, t_end, N_theory_pts)

gamma_e = (1.0 - (v_e / constants.c)**2)**(-0.5)
sigmaB_mbarns = 3.518027889e+01 # approximate total cross section at 40 MeV (from SB table)
sigmaB = sigmaB_mbarns * 1.0e-28
length = zmax

Ngamma_theory = e_den * cu_den * v_e * length * sigmaB * time_theory

ax[0].plot(time_theory * 1e15, Ngamma_theory, '-k', label='Theory')

Martinez_fig9b_data = np.loadtxt("/home/cepheid/TriForce/game_engine/tests/mikes_files/Martinez_fig9b_data.csv", delimiter=",")
ax[1].plot(Martinez_fig9b_data[:,0], Martinez_fig9b_data[:,1], '-k')

me_c2 = constants.electron_mass * constants.speed_of_light**2
# sim_volume = dx * dy * dz
area = (1.0e-6 * 1.6e-6) # since paper simulation is 1D

start = 1
stop = nt
step = 1
for name, cppc, label, marker, color in sims:
    times = np.zeros(nt)
    photon_weights = np.zeros(nt)
    for i in range(start, stop, step):
        file_name = f'{name}/photons_dump_{i:010d}.bp'
        with FileReader(particle_data + file_name) as f:
            times[i] = f.read('Time')[0] * 1.0e15
            photon_weights[i] = f.read('Weight').sum()

    photon_density = photon_weights / area

    ax[0].plot(times, photon_density, c=color, mec=color,
               mew=1.5, ms=3.5, fillstyle="none", linestyle="none", marker=marker, label=label)

    n_bins = 100
    n_count = np.zeros(n_bins)
    dndk = np.zeros(n_bins - 1)
    k_over_gm1_bins = np.logspace(-7, 0, n_bins)

    i_step = 42
    file_name = f'{name}/photons_dump_{i_step:010d}.bp'
    with FileReader(particle_data + file_name) as f:
        weights = f.read('Weight')
        gammas = f.read('Gamma')

    n_photons = gammas.shape[0]

    for j in range(n_photons):
        k_over_gm1 = gammas[j] / (gamma_e - 1.0)
        i_bin = 0
        while i_bin < n_bins - 1 and k_over_gm1 > k_over_gm1_bins[i_bin]:
            i_bin += 1
        n_count[i_bin] += weights[j, 0] / area

    dk = np.zeros(n_bins - 1)
    for j in range(n_bins - 1):
        dk[j] = (k_over_gm1_bins[j + 1] - k_over_gm1_bins[j]) * (gamma_e - 1.0)
        dndk[j] = n_count[j] / dk[j]

    ax[1].plot(k_over_gm1_bins[1:], dndk, c=color, mec=color,
               mew=1.5, ms=3.5, fillstyle="none", linestyle="none", marker=marker, label=label)

ax[0].legend()
fig.tight_layout(pad=0.5, rect=[0,0,1,1])
plt.show()
