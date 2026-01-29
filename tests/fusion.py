#!/usr/bin/env python3

import numpy as np
import subprocess
import matplotlib.pyplot as plt
from scipy import constants
from adios2 import FileReader

from scripts.domain_params import *
from scripts.particle_generation import create_particles

project_path = '/home/cepheid/TriForce/game_engine'
build_path = project_path + '/buildDir'
particle_data = project_path + '/data'

shape = (2, 2, 2)

xmin, xmax = 0.0, 1.0e-6
ymin, ymax = 0.0, 1.0e-6
zmin, zmax = 0.0, 1.0e-6

dx = (xmax - xmin) / (shape[0] - 1)
dy = (ymax - ymin) / (shape[1] - 1)
dz = (zmax - zmin) / (shape[2] - 1)

dt = 5.0e-14
t_end = 1.0e-12
nt = int(t_end / dt) + 1
cfl = constants.c * dt * np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2)

# ===== Particles =====
px_range=(xmin, xmax)
py_range=(ymin, ymax)
pz_range=(zmin, zmax)
T_eV = 5.0e3

deuterium = Particles(
    name='deuterium',
    mass=2.0014 * constants.atomic_mass,
    charge=1,
    atomic_number=1,
    temp=(T_eV, T_eV, T_eV), # eV
    density=1.0e26, # m^-3,
    ppc=(25, 20, 20),
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

tritium = Particles(
    name='tritium',
    mass=3.01605 * constants.atomic_mass,
    charge=1,
    atomic_number=1,
    temp=(T_eV, T_eV, T_eV), # eV
    density=1.0e26, # m^-3,
    ppc=(25, 20, 20),
    distribution='thermal',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

neutrons = Particles(
    name='neutrons',
    mass=constants.m_n,
    charge=0,
    atomic_number=0,
    temp=(0, 0, 0), # eV
    density=0.0, # m^-3,
    ppc=(0, 0, 0),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

helium3 = Particles(
    name='helium3',
    mass=3.016029 * constants.atomic_mass,
    charge=2,
    atomic_number=2,
    temp=(0, 0, 0), # eV
    density=0.0, # m^-3,
    ppc=(0, 0, 0),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

helium4 = Particles(
    name='helium4',
    mass=4.002603 * constants.atomic_mass,
    charge=2,
    atomic_number=2,
    temp=(0, 0, 0), # eV
    density=0.0, # m^-3,
    ppc=(0, 0, 0),
    distribution='none',
    px_range=px_range,
    py_range=py_range,
    pz_range=pz_range
)

# ===== Collisions and Particle Params =====
DD_params = FusionParams(
    products=('neutrons', 'helium3'),
    energy_gain=3.269e6,
    production_multiplier=1.0e10,
    cross_section_file='/tests/cross_section_data/DD_nHe3_BH_eV_m2.txt'
)

DT_params = FusionParams(
    products=('neutrons', 'helium4'),
    energy_gain=17.589e6,
    production_multiplier=1.0e10,
    cross_section_file='/tests/cross_section_data/DT_nHe4_BH_eV_m2.txt'
)

collision_types = ['DD', 'DT']
collision_temps = [5, 10, 19]

has_DD = 'DD' in collision_types
has_DT = 'DT' in collision_types
plot_DD = True
plot_DT = True

for c in collision_types:
    if c == 'DD':
        particles = (deuterium, neutrons, helium3)
        collision = Collision(
            groups=('deuterium', 'deuterium'),
            channels=('fusion',),
            step_interval=1,
            self_scatter=True,
            fusion=DD_params
        )
    else: # DT
        particles = (tritium, deuterium, neutrons, helium4)
        collision = Collision(
            groups=('tritium', 'deuterium'),
            channels=('fusion',),
            step_interval=1,
            fusion=DT_params
        )

    for t in collision_temps:
        sim_name = f'{c}_fusion_{t}keV'
        print(f'Setting up "{sim_name}"')
        TeV = float(t) * 1.0e3
        deuterium.temp = (TeV, TeV, TeV)
        tritium.temp = (TeV, TeV, TeV)

        particle_params = ParticleParams(
            particle_bcs='periodic',
            particle_data=particles,
            collisions=(collision,)
        )

        sim_params = Simulation(
            name=sim_name,
            shape=shape,
            nthreads=1,
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

        create_particles(sim_params, deuterium, particle_data)
        create_particles(sim_params, tritium, particle_data)
        update_header(sim_params, project_path=project_path)

        subprocess.run(
            ['meson', 'compile', '-C', build_path, '-j4'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        ).check_returncode()

        subprocess.run(build_path + '/game_engine').check_returncode()
        print()

if not (plot_DD or plot_DT):
    exit()

def thermal_reactivity_BH1992(T, bg, mrc2, c):
    T2 = T**2
    T3 = T**3
    theta = T / (1.0 - ((T * c[1] + T2 * c[3] + T3 * c[5]) / (1.0 + T * c[2] + T2 * c[4] + T3 * c[6])))
    xsi = (bg**2 / (4.0 * theta))**(1.0/3.0)
    return c[0] * theta * np.sqrt(xsi / (mrc2 * T3)) * np.exp(-3.0 * xsi)

# theory solution for second moment of neutron spectrum (expressed as FWHM) is from
# Ballabio, Kallne, Gorini 1998, Relativistic calculation of fusion product spectra for thermonuclear plasmas
def calc_FWHM_Ballabio1998(T, alpha, w0):
    delta_w = alpha[0] * T ** (2 / 3) / (1.0 + alpha[1] * T ** alpha[2]) + (alpha[3] * T)
    whalf = w0 * (1.0 + delta_w) * np.sqrt(T)

    return whalf

def calc_sigma_fwhm_dyde(weight, gamma, t_final, delta_ij):
    ttl_weight = weight.sum()
    dYdt = ttl_weight / t_final

    sim_volume = 1.0e-6 ** 3
    density = 1e26
    sigma = 1e6 * dYdt * (1.0 + delta_ij) / (density**2 * sim_volume) # m^3 -> cm^3

    E_ttl = (weight * (gamma - 1.0) * constants.m_n * constants.c ** 2 * 1e-3 / constants.e).sum() / ttl_weight
    energy = ((gamma - 1.0) * constants.m_n * constants.c**2 * 1e-3 / constants.e)
    dE2_ttl = ((E_ttl - energy) ** 2 * weight).sum()
    fwhm = np.sqrt(8.0 * np.log(2) * dE2_ttl / ttl_weight)

    return sigma, fwhm


# start plot
plt.style.use('/home/cepheid/TriForce/game_engine/data/triforce.mplstyle')
good_colors = ['#db6d00', '#006ddb', '#920000', '#52a736', '#9B30FF']

fig, ax = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')

ax[0].set_xlabel(r'$T_{\mathrm{i}}$ (keV)')
ax[0].set_ylabel(r'$\langle \sigma v \rangle$ (cm$^3$/s)')
ax[0].set_xlim([0, 20])
ax[0].set_yscale('log')
ax[0].set_ylim([1.0e-22, 1e-15])
ax[0].grid()

ax[1].set_xlabel(r'$T_{\mathrm{i}}$ (keV)')
ax[1].set_ylabel(r'FWHM (keV)')
ax[1].set_xlim([0, 20])
ax[1].set_yscale('log')
ax[1].set_ylim([99.9, 1.5e3])  # set min at 99.9 instead of 100 so intermediate ticks not printed
ax[1].yaxis.set_label_coords(-0.1, 0.5)
ax[1].grid(which='both')

T_keV_theory = np.linspace(0.1, 20, 100)

if has_DD:
    dd_reactivity_theory = thermal_reactivity_BH1992(
        T_keV_theory,
        31.3970,  # Bg, sqrt(keV)
        937814,    # mrc2, keV
        [5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0]
    )

    dd_fwhm_theory = calc_FWHM_Ballabio1998(
        T_keV_theory,
        [1.7013e-3, 0.16888, 0.49, 7.9460e-4],
        82.542 # omega0, sqrt(keV)
    )
    ax[0].plot(T_keV_theory, dd_reactivity_theory, '--k', label='DD Theory')
    ax[1].plot(T_keV_theory, dd_fwhm_theory, '--k', label='DD Theory')

if has_DT:
    dt_reactivity_theory = thermal_reactivity_BH1992(
        T_keV_theory,
        34.3827,  # Bg, sqrt(keV)
        1124656,    # mrc2, keV
        [1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-6]
    )

    dt_fwhm_theory = calc_FWHM_Ballabio1998(
        T_keV_theory,
        [5.1068e-4, 7.6223e-3, 1.78, 8.7691e-5],
        177.259 # omega0, sqrt(keV)
    )

    ax[0].plot(T_keV_theory, dt_reactivity_theory, '-k', label='DT Theory')
    ax[1].plot(T_keV_theory, dt_fwhm_theory, '-k', label='DT Theory')

DD_reactivity_sim = np.zeros(3)
DD_fwhm = np.zeros(3)
DT_reactivity_sim = np.zeros(3)
DT_fwhm = np.zeros(3)

for i, t in enumerate(collision_temps):
    if has_DD:
        # plot DD
        dd_file = f'/DD_fusion_{t}keV/neutrons_dump_{nt-1:010d}.bp'
        with FileReader(particle_data + dd_file) as f:
            final_time = f.read("Time")[0]
            weights = f.read('Weight')
            gammas = f.read('Gamma')

        sigmaDD, fwhmDD = calc_sigma_fwhm_dyde(weights, gammas, final_time, 1.0)
        DD_reactivity_sim[i] = sigmaDD
        DD_fwhm[i] = fwhmDD

    if has_DT:
        # plot DT
        dt_file = f'/DT_fusion_{t}keV/neutrons_dump_{nt-1:010d}.bp'

        with FileReader(particle_data + dt_file) as f:
            final_time = f.read("Time")[0]
            weights = f.read('Weight')
            gammas = f.read('Gamma')

        sigmaDT, fwhmDT = calc_sigma_fwhm_dyde(weights, gammas, final_time, 0.0)
        DT_reactivity_sim[i] = sigmaDT
        DT_fwhm[i] = fwhmDT

if has_DD:
    ax[0].plot(collision_temps, DD_reactivity_sim, marker='o', c=good_colors[3], linestyle='none', ms=10, markeredgewidth=1, markeredgecolor='k', label='D-D')
    ax[1].plot(collision_temps, DD_fwhm, marker='o', c=good_colors[3], linestyle='none', ms=10, markeredgewidth=1, markeredgecolor='k', label='D-D')

if has_DT:
    ax[0].plot(collision_temps, DT_reactivity_sim, marker='o', c=good_colors[2], linestyle='none', ms=10, markeredgewidth=1, markeredgecolor='k', label='D-T')
    ax[1].plot(collision_temps, DT_fwhm, marker='o', c=good_colors[2], linestyle='none', ms=10, markeredgewidth=1, markeredgecolor='k', label='D-T')

ax[0].legend(loc=4)
ax[0].text(0.04, 0.9, '(a)', fontsize=18, weight='bold', transform=ax[0].transAxes)
ax[1].text(0.04, 0.9, '(b)', fontsize=18, weight='bold', transform=ax[1].transAxes)
plt.savefig(particle_data + f'/fusion_test_reactivity_fwhm.png')
plt.close(fig)
# plt.show()
