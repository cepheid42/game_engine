import numpy as np
from scipy import constants
from dataclasses import dataclass
from adios2 import Stream

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

# ===== Position Sampling Utilities =====
# =======================================
def generate_normalized_positions_1d(ppc):
    return np.linspace(1.0 / ppc, 1.0, ppc) - 0.5 / ppc

def generate_normalized_positions_3d(ppc):
    return [points.flatten() for points in np.meshgrid(generate_normalized_positions_1d(ppc[0]),
                                                       generate_normalized_positions_1d(ppc[1]),
                                                       generate_normalized_positions_1d(ppc[2]),
                                                       indexing='ij')]

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

def maxwell_juttner_distribution_alt(mass, temperature, count):
    rng = np.random.default_rng()
    a_MJ, b_MJ, R0_MJ = 0.56, 0.35, 0.95
    T_norm = constants.elementary_charge * temperature / (mass * constants.c**2)
    root2 = np.sqrt(2)

    w = np.array([np.sqrt(constants.pi),
                  a_MJ * np.sqrt(2.0 * T_norm),
                  3.0 * b_MJ * np.sqrt(constants.pi) * T_norm / 2.0,
                  (2.0 * T_norm)**1.5])
    s_sum = w.sum()
    pi = w / s_sum

    def mc_sample(num):
        def R(xx):
            return (1.0 + xx) * np.sqrt(xx + 2.0) / (root2 + (a_MJ * np.sqrt(xx)) + (b_MJ * root2 * xx) + xx ** (3.0 / 2.0))

        x1 = rng.uniform(size=num)
        x2 = rng.uniform(size=num)

        igamma_distribution = np.full(num, 6.0)
        igamma_distribution[x1 < (pi[0] + pi[1] + pi[2])] = 5.0
        igamma_distribution[x1 < (pi[0] + pi[1])] = 4.0
        igamma_distribution[x1 < (pi[0])] = 3.0

        x_variate = rng.gamma(igamma_distribution / 2.0, T_norm)
        acceptance_mask = np.logical_or((x2 < R0_MJ), (x2 < R(x_variate)))

        return x_variate[acceptance_mask]

    x = np.array([])
    remaining_samples = count
    while remaining_samples > 0:
        x = np.concatenate((x, mc_sample(int(np.ceil(1.10 * remaining_samples)))))
        remaining_samples = count - x.shape[0]
    x = x[:count]

    x3 = rng.uniform(size=count)
    x4 = rng.uniform(size=count)
    u_mag = constants.c * np.sqrt(x * (x + 2))
    u = u_mag * np.asarray([(2.0 * x3 - 1.0),
                            2.0 * np.sqrt(x3 * (1.0 - x3)) * np.cos(2.0 * constants.pi * x4),
                            2.0 * np.sqrt(x3 * (1.0 - x3)) * np.sin(2.0 * constants.pi * x4)])

    u = u.T / np.sqrt(1.0 + np.einsum('ij,ij->i', u.T, u.T) / constants.c**2)[:, None]
    velocities = np.empty((count, 3)) # Because reasons
    velocities[:, :] = u
    return velocities

def create_particles(domain, particles, file_dir):
    print(f'Creating {particles.name}...', end=' ', flush=True)
    # todo: add buffer flush

    # ----- Temporary mesh generation -----
    # <AJK> - This will be replaced with an input mesh when we get there
    nx, ny, nz = domain.shape
    dx, dy, dz = domain.deltas

    xmin, xmax = domain.x_range
    x_coords = np.linspace(xmin, xmax, nx, endpoint=True)
    x_center = 0.5 * (x_coords[:-1] + x_coords[1:])

    ymin, ymax = domain.y_range
    y_coords = np.linspace(ymin, ymax, ny, endpoint=True)
    y_center = 0.5 * (y_coords[:-1] + y_coords[1:])

    zmin, zmax = domain.z_range
    z_coords = np.linspace(zmin, zmax, nz, endpoint=True)
    z_center = 0.5 * (z_coords[:-1] + z_coords[1:])

    # ----- Unpack particle params -----
    name = particles.name
    mass = particles.mass
    charge = particles.charge
    temp = particles.temp       # Eventually, these will have to be arrays, or maybe functions?
    density = particles.density # ^
    ppc_x, ppc_y, ppc_z = particles.ppc
    distribution = particles.distribution
    px_min, px_max = particles.px_range
    py_min, py_max = particles.py_range
    pz_min, pz_max = particles.pz_range

    # ----- Generate particle positions -----
    # "Fake" geometry generation (will be replaced when full geometry support is added)
    xx, yy, zz = np.meshgrid(x_center, y_center, z_center, indexing='ij')

    fake_geom = (xx >= px_min)
    fake_geom = np.logical_and(fake_geom, (xx <= px_max))
    fake_geom = np.logical_and(fake_geom, (yy >= py_min))
    fake_geom = np.logical_and(fake_geom, (yy <= py_max))
    fake_geom = np.logical_and(fake_geom, (zz >= pz_min))
    fake_geom = np.logical_and(fake_geom, (zz <= pz_max))
    fake_geom = np.argwhere(fake_geom)

    # Get cell positions
    xc = x_coords[fake_geom[:, 0]]
    xcp1 = x_coords[fake_geom[:, 0] + 1]

    yc = y_coords[fake_geom[:, 1]]
    ycp1 = y_coords[fake_geom[:, 1] + 1]

    zc = z_coords[fake_geom[:, 2]]
    zcp1 = z_coords[fake_geom[:, 2] + 1]

    x_norm, y_norm, z_norm = generate_normalized_positions_3d(particles.ppc) # Get particle positions in each cell

    # Generate all particle positions
    xp = (np.matmul(xc[:, None], (1.0 - x_norm)[:, None].T) + np.matmul(xcp1[:, None], x_norm[:, None].T)).flatten()
    yp = (np.matmul(yc[:, None], (1.0 - y_norm)[:, None].T) + np.matmul(ycp1[:, None], y_norm[:, None].T)).flatten()
    zp = (np.matmul(zc[:, None], (1.0 - z_norm)[:, None].T) + np.matmul(zcp1[:, None], z_norm[:, None].T)).flatten()

    positions = np.vstack((xp, yp, zp)).T
    # print(positions.shape)

    # Generate particle weights
    ppc_total = ppc_x * ppc_y * ppc_z
    num_particles = ppc_total * fake_geom.shape[0]

    weights = np.full((num_particles, 1), density * dx * dy * dz / ppc_total)
    # print(weights.shape)

    # Sample particle velocities
    if distribution == 'thermal':
        velocities = thermal_distribution(mass, temp, num_particles)
    else:
        velocities = maxwell_juttner_distribution_alt(mass, temp, num_particles)
        # velocities = maxwell_juttner_distribution(mass, temp, num_particles)

    # print(velocities.shape)
    output = np.hstack((positions, velocities, weights))
    # header = f'{name} {mass} {charge}'
    # np.savetxt(file_dir + f'/{name}.dat', output, delimiter=' ', header=header)
    with Stream(file_dir + f'/{name}.bp', 'w') as f:
        f.write_attribute('name', name)
        f.write_attribute('mass', mass)
        f.write_attribute('charge', charge)
        f.write_attribute('num_particles', num_particles)
        f.write('particles', output, output.shape, [0, 0], output.shape)
    print('Done')
