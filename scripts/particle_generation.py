from collections.abc import Callable, Iterable

import numpy as np
from scipy import constants
from adios2 import Stream


# ===== Position Sampling Utilities =====
# =======================================
def generate_normalized_positions_1d(ppc):
    return np.linspace(1.0 / ppc, 1.0, ppc) - 0.5 / ppc

def generate_normalized_positions_3d(ppc):
    return [points.flatten() for points in np.meshgrid(generate_normalized_positions_1d(ppc[0]),
                                                       generate_normalized_positions_1d(ppc[1]),
                                                       generate_normalized_positions_1d(ppc[2]),
                                                       indexing='ij')]


def generate_density_mapping(density : float | Callable | Iterable, xcenter, ycenter, zcenter):
    assert(xcenter.shape == ycenter.shape == zcenter.shape) # check shapes for self-consistency

    # if centers aren't already meshed, mesh them
    if xcenter.ndim == 1:
        xmesh, ymesh, zmesh = np.meshgrid(xcenter, ycenter, zcenter, sparse=True)
    elif xcenter.ndim == 3:
        xmesh, ymesh, zmesh = xcenter, ycenter, zcenter
    else:
        raise ValueError

    if isinstance(density, Iterable):
        return density
    elif isinstance(density, Callable):
        return density(xmesh, ymesh, zmesh)
    else:
        return np.full(xmesh.shape, density)


def constant_distribution(mass, temp, num_particles):
    temp = np.asarray(temp)
    vx, vy, vz = constants.c * np.sqrt(1 - (1 + constants.eV * temp / (mass * constants.c**2))**-2)
    velocities = np.zeros((num_particles, 3))
    velocities[:, 0] = vx
    velocities[:, 1] = vy
    velocities[:, 2] = vz
    return velocities


def thermal_distribution(mass, T_M, num_particles, velocity=0.0):
    rng = np.random.default_rng()
    v_thermal = np.sqrt(constants.elementary_charge * T_M / mass)
    velocities = rng.normal(velocity, v_thermal, (num_particles, 3))
    return velocities


def maxwell_juttner_distribution(mass, temperature, count):
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

def write_particle_file(name, file_dir, mass, charge, positions, velocities, weights, gammas, atomic_number = 0, tracer = False, sourcer = False):
    with Stream(file_dir + f'/{name.lower().replace(" ","_")}.bp', 'w') as f:
        f.write_attribute('Name', name)
        f.write_attribute('Mass', mass)
        f.write_attribute("Mass/Unit", "kg")
        f.write_attribute('Charge', charge)
        f.write_attribute("Charge/Unit", "C")
        f.write_attribute("Atomic Number", np.array([atomic_number],np.uint64)[0])
        # f.write_attribute('Tracer', np.array([tracer],np.uint64)[0])
        # f.write_attribute('Sourcer', np.array([sourcer],np.uint64)[0])
        f.write("Position", positions.copy(), positions.shape, [0, 0], positions.shape)
        f.write_attribute("Unit", "m", "Position")
        f.write("Velocity", velocities.copy(), velocities.shape, [0, 0], velocities.shape)
        f.write_attribute("Unit", "m/s", "Velocity")
        f.write("Weight", weights, weights.shape, [0], weights.shape)
        f.write("Gamma", gammas, gammas.shape, [0], gammas.shape)



def create_particles(domain, particles, file_dir):
    print(f'Creating {particles.name}...', end=' ', flush=True)

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
    atomic_number = particles.atomic_number
    temp = particles.temp       # Eventually, these will have to be arrays, or maybe functions?
    density = particles.density # ^
    ppc_x, ppc_y, ppc_z = particles.ppc
    distribution = particles.distribution
    px_min, px_max = particles.px_range
    py_min, py_max = particles.py_range
    pz_min, pz_max = particles.pz_range

    if particles.distribution == 'none':
        return

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

    if not fake_geom.shape[0]:
        print(f'!!! Warning: Particle Geometry Empty !!!', end=' ', flush=True)
        print(f'Done.', end='\n', flush=True)
        return

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

    # Generate particle weights
    ppc_total = ppc_x * ppc_y * ppc_z
    num_particles = ppc_total * fake_geom.shape[0]

    density_mesh = generate_density_mapping(density,
                                            x_center[fake_geom[:, 0]],
                                            y_center[fake_geom[:, 1]],
                                            z_center[fake_geom[:, 2]])
    weights = (dx * dy * dz / ppc_total) * np.matmul(density_mesh.flatten()[:, None], np.ones((ppc_total,))[:, None].T).flatten()

    # Sample particle velocities
    if distribution == 'thermal':
        velocities = thermal_distribution(mass, temp[0], num_particles)
    elif distribution == 'constant':
        velocities = constant_distribution(mass, temp, num_particles)
    else:
        velocities = maxwell_juttner_distribution(mass, temp[0], num_particles)

    gammas = 1.0 / np.sqrt(1 - ((velocities/constants.c)**2).sum(axis=1))
    write_particle_file(name, file_dir, mass, charge, positions, velocities, weights, gammas, atomic_number, False, False)
    print('Done')
