import numpy as np
from scipy import constants
from dataclasses import dataclass

@dataclass
class Particles:
    name: str
    mass: float
    charge: float
    temp: float
    density: float
    ppc: tuple
    distribution: str = 'relativistic'
    px_range: tuple = ()
    py_range: tuple = ()
    pz_range: tuple = ()


def thermal_distribution(mass, T_M, num_particles, velocity=0.0):
    rng = np.random.default_rng()
    v_thermal = np.sqrt(constants.elementary_charge * T_M / mass)
    velocities = rng.normal(velocity, v_thermal, (num_particles, 3))
    return velocities


def maxwell_juttner_distribution(mass, t_M, num_particles):
    rng = np.random.default_rng()
    t_norm = constants.elementary_charge * t_M / (mass * constants.c**2)
    a, b, r0 = 0.56, 0.35, 0.95
    root2 = np.sqrt(2)
    w3 = np.sqrt(np.pi)
    w4 = a * np.sqrt(2 * t_norm)
    w5 = 1.5 * np.sqrt(np.pi) * b * t_norm
    w6 = (2 * t_norm)**1.5
    s_sum = w3 + w4 + w5 + w6
    pi3 = w3 / s_sum
    pi4 = w4 / s_sum
    pi5 = w5 / s_sum
    def sample():
        def R(xx):
            return (1 + xx) * np.sqrt(xx + 2) / (root2 + (a * np.sqrt(xx)) + (b * root2 * xx) + xx ** (3 / 2))
        while True:
            X1 = rng.uniform()
            X2 = rng.uniform()
            i = 6
            if X1 < pi3:
                i = 3
            elif X1 < pi3 + pi4:
                i = 4
            elif X1 < pi3 + pi4 + pi5:
                i = 5
            x = rng.gamma(i / 2, t_norm)
            if X2 < r0 or X2 < R(x):
                break

        x3 = rng.uniform()
        x4 = rng.uniform()
        u_mag = constants.c * np.sqrt(x * (x + 2))
        u = u_mag * np.asarray([(2 * x3 - 1),
                                2.0 * np.sqrt(x3 * (1.0 - x3)) * np.cos(2.0 * np.pi * x4),
                                2.0 * np.sqrt(x3 * (1.0 - x3)) * np.sin(2.0 * np.pi * x4)])
        return u / np.sqrt(1 + (u @ u) / constants.c**2)

    velocities = np.zeros((num_particles, 3))
    for p in range(num_particles):
        velocities[p, :] = sample()
    return velocities


def create_particles(domain, particles, file_dir, plot=False):
    print(f'Creating {particles.name}...', end=' ')
    xmin, xmax = domain.x_range
    ymin, ymax = domain.y_range
    zmin, zmax = domain.z_range
    dx, dy, dz = domain.deltas

    ppc_x, ppc_y, ppc_z = particles.ppc
    px_min, px_max = particles.px_range
    py_min, py_max = particles.py_range
    pz_min, pz_max = particles.pz_range
    distribution = particles.distribution
    density = particles.density
    charge = particles.charge
    mass = particles.mass
    temp = particles.temp
    name = particles.name

    pnx = int((px_max - px_min) / dx)
    pny = int((py_max - py_min) / dy)
    pnz = int((pz_max - pz_min) / dz)

    pdx = dx / ppc_x
    pdy = dy / ppc_y
    pdz = dz / ppc_z

    xs = np.arange(xmin, xmax, dx)
    ys = np.arange(ymin, ymax, dy)
    zs = np.arange(zmin, zmax, dz)

    if xs.shape[0] > 2:
        px0 = xs[np.searchsorted(xs, px_min)]
        px1 = xs[np.searchsorted(xs, px_max)]
        pxs = np.arange(px0 + pdx / 2, px1, pdx)
    else:
        pxs = np.array([dx / 2])

    if ys.shape[0] > 2:
        py0 = ys[np.searchsorted(ys, py_min)]
        py1 = ys[np.searchsorted(ys, py_max)]
        pys = np.arange(py0 + pdy / 2, py1, pdy)
    else:
        pys = np.array([dy / 2])

    if zs.shape[0] > 2:
        pz0 = zs[np.searchsorted(zs, pz_min)]
        pz1 = zs[np.searchsorted(zs, pz_max)]
        pzs = np.arange(pz0 + pdz / 2, pz1, pdz)
    else:
        pzs = np.array([dz / 2])

    xc, yc, zc = np.meshgrid(pxs, pys, pzs)

    ppc_total = ppc_x * ppc_y * ppc_z
    num_particles = ppc_total * pnx * pny * pnz

    weights = np.full((num_particles, 1), density * dx * dy * dz / ppc_total)
    positions = np.vstack((xc.flatten(), yc.flatten(), zc.flatten())).T
    
    if distribution == 'thermal':
        velocities = thermal_distribution(mass, temp, num_particles)
    else:
        velocities = maxwell_juttner_distribution(mass, temp, num_particles)

    if plot:
        import matplotlib.pyplot as plt
        plt.scatter(positions[:, 0], positions[:, 2], s=2)
        plt.xticks(xs)
        plt.yticks(zs)
        plt.grid()
        plt.show()
    else:
        output = np.hstack((positions, velocities, weights))
        header = f'{name} {mass} {charge}'
        np.savetxt(file_dir + f'/{name}.dat', output, delimiter=' ', header=header)
    print('Done')
