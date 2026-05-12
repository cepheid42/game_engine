from collections.abc import Callable, Iterable

import numpy as np
from scipy import constants
from adios2 import Stream

from .particle import ParticleDistributionType, ParticleGroup

# ===== Position Sampling Utilities =====
# =======================================
def generate_normalized_positions_1d(ppc):
    return np.linspace(1.0 / ppc, 1.0, ppc) - 0.5 / ppc


def generate_normalized_positions_3d(ppc):
    return [points.flatten() for points in np.meshgrid(generate_normalized_positions_1d(ppc[0]),
                                                       generate_normalized_positions_1d(ppc[1]),
                                                       generate_normalized_positions_1d(ppc[2]),
                                                       indexing='ij')]


def generate_scalar_mapping(scalar : float | Callable | Iterable, xcenter, ycenter, zcenter):
    assert(xcenter.shape == ycenter.shape == zcenter.shape) # check shapes for self-consistency

    # if centers aren't already meshed, mesh them
    if xcenter.ndim == 1:
        xmesh, ymesh, zmesh = np.meshgrid(xcenter, ycenter, zcenter, sparse=True)
    elif xcenter.ndim == 3:
        xmesh, ymesh, zmesh = xcenter, ycenter, zcenter
    else:
        raise ValueError

    if isinstance(scalar, Iterable):
        return scalar
    elif isinstance(scalar, Callable):
        return scalar(xmesh, ymesh, zmesh)
    else:
        return np.full(xmesh.shape, scalar)


def sample_maxwell_distribution(mass, temperature, count, v_drift=np.zeros(3)):
    rng = np.random.default_rng()
    v_thermal = np.sqrt(constants.elementary_charge * temperature / mass)

    velocities = np.matmul(np.ones((count,))[:, None], v_drift[None, :])

    if v_thermal != 0.0:
        vrand = rng.normal(0.0, v_thermal, (count, 3))
        v_avg = vrand.mean(axis=0)
        v2_avg = np.square(vrand).mean(axis=0)
        denom = np.sqrt(np.abs(v2_avg) - np.square(v_avg))

        velocities += (v_thermal / denom * (vrand - v_avg))

    return np.clip(velocities, -0.9999 * constants.c, 0.9999 * constants.c)


def sample_maxwell_juttner_distribution(mass, temperature, count, v_drift=np.zeros(3)):
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

    velocities += np.matmul(np.ones((count,))[:, None], v_drift[None, :])
    return velocities


def generate_particle_positions(geometry : Iterable, x_coords, y_coords, z_coords, ppc):
    # Get cell positions
    xc = x_coords[geometry[:, 0]]
    xcp1 = x_coords[geometry[:, 0] + 1]

    yc = y_coords[geometry[:, 1]]
    ycp1 = y_coords[geometry[:, 1] + 1]

    zc = z_coords[geometry[:, 2]]
    zcp1 = z_coords[geometry[:, 2] + 1]

    # Get cell-normalized positions
    x_norm, y_norm, z_norm = generate_normalized_positions_3d(ppc) # Get particle positions in each cell

    # Generate all particle positions
    xp = (np.matmul(xc[:, None], (1.0 - x_norm)[None, :]) + np.matmul(xcp1[:, None], x_norm[None, :])).flatten()
    yp = (np.matmul(yc[:, None], (1.0 - y_norm)[None, :]) + np.matmul(ycp1[:, None], y_norm[None, :])).flatten()
    zp = (np.matmul(zc[:, None], (1.0 - z_norm)[None, :]) + np.matmul(zcp1[:, None], z_norm[None, :])).flatten()

    return np.vstack((xp, yp, zp)).T


def generate_particle_weights(density, x_center, y_center, z_center, deltas, ppc_total):
    dx, dy, dz = deltas

    density_mesh = generate_scalar_mapping(density, x_center, y_center, z_center)
    
    return (dx * dy * dz / ppc_total) * np.matmul(density_mesh.flatten()[:, None], np.ones((ppc_total,))[None, :]).flatten()


def generate_particle_velocities(mass, temp, x_center, y_center, z_center, ppc_total, v_drift=np.zeros(3), relativistic=False):
    temp_mesh = generate_scalar_mapping(temp, x_center, y_center, z_center)

    if relativistic:
        velocities = np.array([sample_maxwell_juttner_distribution(mass, T, ppc_total, v_drift) for T in temp_mesh.flatten()]).reshape((-1, 3))
    else:
        velocities = np.array([sample_maxwell_distribution(mass, T, ppc_total, v_drift) for T in temp_mesh.flatten()]).reshape((-1, 3))

    return velocities

def write_particle_file(name, filename, mass, charge, positions, velocities, weights, gammas, atomic_number = 0, tracer = False, sourcer = False):
    with Stream(filename, 'w') as f:
        f.write_attribute('Name', name)
        f.write_attribute('Mass', mass)
        f.write_attribute("Mass/Unit", "kg")
        f.write_attribute('Charge', charge)
        f.write_attribute("Charge/Unit", "C")
        f.write_attribute("Atomic Number", np.array([atomic_number],np.uint64)[0])
        f.write_attribute('Tracer', np.array([tracer],np.uint64)[0])
        f.write_attribute('Sourcer', np.array([sourcer],np.uint64)[0])
        f.write("Position", positions.copy(), positions.shape, [0,0], positions.shape)
        f.write_attribute("Unit", "m", "Position")
        f.write("Velocity", velocities.copy(), velocities.shape, [0,0], velocities.shape)
        f.write_attribute("Unit", "m/s", "Velocity")
        f.write("Weight", weights, weights.shape, [0], weights.shape)
        f.write("Gamma", gammas, gammas.shape, [0], gammas.shape)

def create_particles(simulation):
    extra_groups = []

    xdim = simulation.grid_params.x
    x_coords = xdim.coords if xdim is not None else np.array((0.0,))
    x_center = 0.5 * (x_coords[:-1] + x_coords[1:]) if xdim is not None else np.array((0.0,))

    ydim = simulation.grid_params.y
    y_coords = ydim.coords if ydim is not None else np.array((0.0,))
    y_center = 0.5 * (y_coords[:-1] + y_coords[1:]) if ydim is not None else np.array((0.0,))

    zdim = simulation.grid_params.z
    z_coords = zdim.coords
    z_center = 0.5 * (z_coords[:-1] + z_coords[1:]) if zdim is not None else np.array((0.0,))

    for group in simulation.particle_params.groups:

        grid = simulation.grid_params.get_mesh("geom", staggered_fields=True)
        group.geometry.voxelize(grid)
        cell_map = np.argwhere(group.geometry.boolean_map).reshape(-1,3)

        positions = generate_particle_positions(cell_map, x_coords, y_coords, z_coords, group.ppc)

        # ----- Generate particle weight -----
        ppc_total = np.prod(group.ppc)

        weights = generate_particle_weights(
            group.density,
            x_center[cell_map[:, 0]],
            y_center[cell_map[:, 1]],
            z_center[cell_map[:, 2]],
            simulation.grid_params.delta,
            ppc_total
        )

        # ----- Generate particle velocities -----
        velocities = generate_particle_velocities(
            group.mass,
            group.temp,
            x_center[cell_map[:, 0]],
            y_center[cell_map[:, 1]],
            z_center[cell_map[:, 2]],
            ppc_total,
            relativistic=(group.distribution is ParticleDistributionType.Relativistic)
        )

        # print(velocities.shape)
        gammas = 1.0 / np.sqrt(1 - ((velocities/constants.c)**2).sum(axis=1))

        num_particles = weights.shape[0]
        write_particle_file(
            group.name,
            group.filename,
            group.mass,
            group.charge,
            positions,
            velocities,
            weights,
            gammas,
            group.atomic_number,
            False,
            False
        )

        # check if tracer particles should be generated for this particle group
        if group.tracer_fraction > 0.0:
            tracer_group = ParticleGroup(
                name = f"{group.name} Tracers",
                filename = group.filename.rpartition("/")[0],
                mass = group.mass,
                charge = group.charge,
                atomic_number = group.atomic_number,
                temp = group.temp,
                density = group.density,
                ppc = group.ppc,
                tracer = True
            )

            num_tracers = int(num_particles * group.tracer_fraction)
            rng = np.random.default_rng()
            pids = np.arange(0,num_particles)
            chosen_pids = rng.choice(pids, num_tracers, replace=False)

            write_particle_file(tracer_group.name, tracer_group.filename, tracer_group.mass, tracer_group.charge, positions[chosen_pids,:], velocities[chosen_pids,:], weights[chosen_pids], gammas[chosen_pids], tracer_group.atomic_number, True, False)

            extra_groups.append(tracer_group)

    # add any extra tracer particle groups to the list of groups to load
    simulation.particle_params.groups=tuple([*simulation.particle_params.groups, *extra_groups])
