#!/usr/bin/env python3

from adios2 import FileReader
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Circle
from scipy import constants


data_dir = '/home/cepheid/TriForce/game_engine/data'

EPS0 = 8.8541878188E-12 # F/m
MU0 = 1.25663706127E-6 # N/A^2

# dt = 1.829541541469147e-11
dt = 4.5738543e-12

def load_field(n, name, file_dir):
    filename = f'/fields_{n:010d}.bp'
    with FileReader(data_dir + file_dir + filename) as f:
        return f.read(name)

def plot3d(n, name, step, file_dir):
    print(f'Processing file {n}')
    vmin, vmax = -0.1, 0.1
    # norm = Normalize(vmin=vmin, vmax=vmax)

    frame = load_field(n, name, file_dir)

    fig, ax = plt.subplots(figsize=(10, 10))

    nx = frame.shape[0] // 2
    ny = frame.shape[1] // 2
    nz = frame.shape[2] // 2

    im = ax.contourf(frame[:, :, nz].T, levels=100, vmin=vmin, vmax=vmax)
    #im = ax.pcolormesh(frame[:, :, nz].T, vmin=vmin, vmax=vmax)

    fig.colorbar(ScalarMappable(norm=im.norm), ax=ax)
    # fig.colorbar(im)

    ax.set_xlabel('x (arb)')
    ax.set_ylabel('y (arb)')
    ax.set_title(f'{name}[:, ny/2, :]')

    plt.savefig(data_dir + f'/pngs/{name}_{n // step:08}.png')
    plt.clf()
    plt.close(fig)

def total_field_energy(n, file_dir):
    Ex = load_field(n, 'Ex', file_dir)[:, :-1, :-1]
    Ey = load_field(n, 'Ey', file_dir)[:-1, :, :-1]
    Ez = load_field(n, 'Ez', file_dir)[:-1, :-1, :]
    Hx = load_field(n, 'Hx', file_dir)[:-1, :, :]
    Hy = load_field(n, 'Hy', file_dir)[:, :-1, :]
    Hz = load_field(n, 'Hz', file_dir)[:, :, :-1]

    E = 0.5 * EPS0 * (Ex**2 + Ey**2 + Ez**2)
    H = (0.5 * MU0) * (Hx**2 + Hy**2 + Hz**2)

    total = E + H
    return E.sum() * (0.01**3), H.sum() * (0.01**3), total.sum() * (0.01**3)

def plot_total_energy(start, stop, step):
    def plot_regular(arr1, arr2, log=False):
        times = np.linspace(0, stop * dt, len(arr1))
        fig, ax = plt.subplots(figsize=(8,8))

        # ax.plot(times, arr1[:, 0], label='E-Field float')
        # ax.plot(times, arr1[:, 1], label='B-Field float')
        # ax.plot(times, arr2[:, 0], label='E-Field double')
        # ax.plot(times, arr2[:, 1], label='B-Field double')
        ax.plot(times, arr1[:, 2], label='Total single')
        ax.plot(times, arr2[:, 2], label='Total double')

        if log:
            ax.set_yscale('log')

        ax.grid(True, which='major', axis='both', linestyle='--')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Energy (J)')
        ax.set_title(f'Total Field Energy')
        ax.legend()
        plt.savefig(data_dir + f'/pngs/total_energy.png')
        plt.close(fig)

    targs1 = [(n, '/single_test') for n in range(start, stop + step, step)]
    targs2 = [(n, '/double_test') for n in range(start, stop + step, step)]

    with mp.Pool(16) as p:
        result1 = p.starmap(total_field_energy, targs1)
        result2 = p.starmap(total_field_energy, targs2)

    plot_regular(np.asarray(result1), np.asarray(result2), True)


def compare_fields(n, name, step):
    print(f'Processing file {n}')
    # vmin, vmax = -0.1, 0.1

    single = load_field(n, name, '/single_test')
    double = load_field(n, name, '/double_test')

    diff = double - single

    fig, ax = plt.subplots(figsize=(10, 10))

    nx = diff.shape[0] // 2
    ny = diff.shape[1] // 2
    nz = diff.shape[2] // 2

    im = ax.contourf(diff[:, :, nz].T, levels=100)#, vmin=vmin, vmax=vmax)
    #im = ax.pcolormesh(frame[:, :, nz].T, vmin=vmin, vmax=vmax)

    # fig.colorbar(ScalarMappable(norm=im.norm), ax=ax)
    fig.colorbar(im)

    ax.set_xlabel('x (arb)')
    ax.set_ylabel('y (arb)')
    ax.set_title(f'{name}[:, ny/2, :]')

    plt.savefig(data_dir + f'/pngs/diff_{name}_{n // step:08}.png')
    plt.clf()
    plt.close(fig)


def plot_particles(start, stop, step, group_name, file_dir):
    def v_from_KE(KE = 2, mass = constants.m_e, return_rel_v = False):
        rest_mass = mass * constants.c**2
        gamma = gamma_from_KE(KE,mass)
        rel_p = (((gamma*rest_mass)**2-rest_mass**2)**0.5)/constants.c
        if return_rel_v:
            return rel_p/mass
        return rel_p/mass/gamma


    def gamma_from_KE(KE = 2, mass = constants.m_e):
        return KE*constants.e/mass/constants.c**2 + 1

    def load_particles():
        frames = {
            'location': [],
            'velocity': [],
            'weight': []
        }
        for n in range(start, stop, step):
            file = f'/{group_name}_{n:010d}.bp'
            # print(file)
            with FileReader(data_dir + file_dir + file) as f:
                frames['location'].append(f.read('Position'))
                frames['velocity'].append(f.read('Velocity'))
                frames['weight'].append(f.read('Weight'))

        frames['location'] = np.asarray(frames['location'])
        frames['velocity'] = np.asarray(frames['velocity'])
        frames['weight'] = np.asarray(frames['weight'])
        return frames

    total_time = stop * dt

    particles = load_particles()

    particle_E = 1000
    Bz_mag = 1600e-4

    # analytic relativistic gyroradius and gyroperiod
    gyroradius = v_from_KE(particle_E, return_rel_v = True) * constants.m_e / constants.e / Bz_mag
    gyroperiod = 1 / (constants.e * Bz_mag / constants.m_e / 2 / np.pi) * gamma_from_KE(particle_E)
    steps_per_gyroperiod = gyroperiod / dt
    complete_gyroperiods = int(total_time / gyroperiod)

    # obtain the centroid of the orbit, looking at only complete gyro orbits
    end_frame = int(complete_gyroperiods * steps_per_gyroperiod / step) + 1

    center = particles["location"][:end_frame,0,:].mean(axis=0)

    print(f"center: {center}, gyroradius: {gyroradius}")

    simulated_gyroradii=(((particles["location"][:,0,:]-center)**2).sum(axis=1))**0.5

    fig,ax = plt.subplots(figsize=(10,10))
    # gyrocircle = Circle(center*100,gyroradius*100, facecolor=(0,0,0,0), lw=2, edgecolor="tab:orange", label="analytic")
    # plt.gca().add_artist(gyrocircle)
    plt.plot(particles["location"][:end_frame,0,0]*100,particles["location"][:end_frame,0,1]*100,marker='+', markersize = 10,lw=0, label="simulated")
    plt.xlabel("x (cm)")
    plt.ylabel("y (cm)")
    plt.title("simulated vs analytic gyro-orbit")
    plt.legend(title="gyro orbit")
    plt.grid(ls="--")
    filename = data_dir + f"/pngs/{group_name}_gyroorbit"
    # plt.savefig(f"{filename}.pdf", bbox_inches="tight", pad_inches=0.05)
    plt.savefig(f"{filename}.png", bbox_inches="tight", pad_inches=0.05)
    plt.close()

    # fig,ax = plt.subplots(figsize=(12,10))
    # plt.plot((simulated_gyroradii-gyroradius)/gyroradius)
    # plt.title("fractional error in gyroradius vs step")
    # plt.xlabel("step")
    # plt.ylabel("fractional error")
    # plt.grid(ls="--")
    # filename = data_dir + f"/pngs/{group_name}_gyroradius_error"
    # # plt.savefig(f"{filename}.pdf", bbox_inches="tight", pad_inches=0.05)
    # plt.savefig(f"{filename}.png", bbox_inches="tight", pad_inches=0.05)

    print(f"mean fractional error in gyro orbit: {((simulated_gyroradii-gyroradius)/gyroradius).mean()}")



def main():
    step = 50
    start = 0
    stop = 1000

    file_dir = '/particle_test'

    plot_particles(start, stop, step, 'electrons', file_dir)

    # targs = [(n, 'electrons', 'density', step) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(plot_metric, targs)

    # field = 'Ez'
    # targs = [(n, field, step, file_dir) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(plot3d, targs)

    # field = 'Ez'
    # targs = [(n, field, step) for n in range(start, stop + step, step)]
    # with mp.Pool(16) as p:
    #     p.starmap(compare_fields, targs)

    # plot_total_energy(start, stop, step)


if __name__ == '__main__':
    main()