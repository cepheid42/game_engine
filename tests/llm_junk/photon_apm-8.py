"""
Batched exact-moment adaptive particle manager for photon macroparticles.

This standalone example initializes a 1 keV blackbody photon population,
reduces the photon count from 200 to 100, and preserves the following
quantities to floating-point roundoff by changing only photon weights:

    1. Total macro-weight.
    2. Total photon energy.
    3. All three components of total photon momentum.
    4. The full symmetric radiation pressure tensor.
    5. The number-weighted x, y, and z spatial centroids.
    6. The second and third photon-energy spectral moments.

The accelerated implementation removes many photons in each pass by solving
several independent small linear systems in one batched NumPy call.  It uses
only NumPy and Matplotlib; no SciPy, KD-tree, Numba, or other special library
is required.
"""

# Import the high-resolution timer used to report the APM execution time.
import time
# Import NumPy for arrays, linear algebra, random numbers, and statistics.
import numpy as np
# Import Matplotlib for the same before-and-after plots used previously.
import matplotlib.pyplot as plt

# -----------------------------
# Physical constants (SI units)
# -----------------------------

# Store the Boltzmann constant in joules per kelvin.
k_B = 1.380649e-23
# Store Planck's constant in joule-seconds.
h = 6.62607015e-34
# Store the speed of light in meters per second.
c = 299792458.0
# Store the joule conversion factor for one electron-volt.
eV_to_J = 1.602176634e-19
# Store the Stefan-Boltzmann constant in SI units.
sigma_sb = 5.670374419e-8

# -----------------------------
# Problem setup
# -----------------------------

# Set the initial number of photon macroparticles in this example cell.
N_particles = 200
# Set the target number of photon macroparticles after the APM reduction.
N_target = 100
# Set the initial blackbody thermal scale kT in keV.
kT_keV = 1.0
# Define the cell size in the code-coordinate system.
cell_size = np.array([1.0, 1.0, 1.0])
# Define the physical cell volume used to normalize the radiation energy.
cell_volume_m3 = 1.0e-6
# Set the random-number seed so this example is exactly reproducible.
random_seed = 0
# Set the maximum number of random regrouping attempts per batched pass.
maximum_batch_retries = 50
# Set the relative threshold used to classify a weight as numerically zero.
relative_zero_tolerance = 1.0e-13
# Set the relative residual allowed for a computed local null vector.
relative_null_tolerance = 1.0e-10

# Convert kT from keV to joules.
kT_J = kT_keV * 1.0e3 * eV_to_J
# Convert the photon thermal-energy scale to temperature in kelvin.
T_K = kT_J / k_B
# Compute the radiation constant a = 4 sigma_SB / c.
a_rad = 4.0 * sigma_sb / c
# Compute the blackbody energy density u = a T^4.
u_rad = a_rad * T_K**4
# Compute the target radiation energy represented by this cell.
E_total_cell = u_rad * cell_volume_m3

# -----------------------------
# Photon initialization
# -----------------------------

# Construct a tabulated cumulative distribution for Planck photon energies.
def make_planck_sampler(x_max=20.0, n_grid=5000):
    # Build a dimensionless energy grid x = E/(kT) without including x = 0.
    x_grid = np.linspace(1.0e-6, x_max, n_grid)
    # Evaluate the unnormalized Planck photon-number density x^2/(exp(x)-1).
    probability_density = x_grid**2 / np.expm1(x_grid)
    # Approximate the cumulative distribution with a cumulative sum.
    cumulative_distribution = np.cumsum(probability_density)
    # Normalize the cumulative distribution so that it ends at one.
    cumulative_distribution /= cumulative_distribution[-1]
    # Return the tabulated energy grid and cumulative distribution.
    return x_grid, cumulative_distribution


# Sample photon energies from a tabulated Planck cumulative distribution.
def sample_planck_energies_keV(
    number_of_photons,
    thermal_energy_keV,
    x_grid,
    cumulative_distribution,
    rng,
):
    # Draw one uniform random variate for each photon.
    uniform_samples = rng.random(number_of_photons)
    # Invert the cumulative distribution by linear interpolation.
    dimensionless_energies = np.interp(
        uniform_samples,
        cumulative_distribution,
        x_grid,
    )
    # Convert the dimensionless energies to physical energies in keV.
    energy_keV = dimensionless_energies * thermal_energy_keV
    # Return the sampled photon energies.
    return energy_keV


# Sample isotropic unit vectors uniformly over the sphere.
def sample_isotropic_directions(number_of_directions, rng):
    # Draw random numbers used to sample the polar angle.
    polar_samples = rng.random(number_of_directions)
    # Draw random numbers used to sample the azimuthal angle.
    azimuth_samples = rng.random(number_of_directions)
    # Sample cos(theta) uniformly on the interval [-1, 1].
    cosine_theta = 2.0 * polar_samples - 1.0
    # Recover sin(theta) while protecting against tiny negative roundoff.
    sine_theta = np.sqrt(np.maximum(0.0, 1.0 - cosine_theta**2))
    # Sample the azimuthal angle uniformly on [0, 2 pi).
    azimuth = 2.0 * np.pi * azimuth_samples
    # Compute the x component of every unit vector.
    direction_x = sine_theta * np.cos(azimuth)
    # Compute the y component of every unit vector.
    direction_y = sine_theta * np.sin(azimuth)
    # Use cos(theta) as the z component.
    direction_z = cosine_theta
    # Combine the three components into an N-by-3 array.
    directions = np.column_stack((direction_x, direction_y, direction_z))
    # Return the isotropically sampled directions.
    return directions


# Initialize positions, energies, momenta, and weights for one photon cell.
def initialize_photons_in_cell(number_of_photons, thermal_energy_keV, rng):
    # Sample photon positions uniformly within the cell.
    positions = rng.random((number_of_photons, 3)) * cell_size
    # Build the inverse-CDF sampler for the Planck spectrum.
    x_grid, cumulative_distribution = make_planck_sampler()
    # Sample photon energies from the Planck distribution.
    energy_keV = sample_planck_energies_keV(
        number_of_photons,
        thermal_energy_keV,
        x_grid,
        cumulative_distribution,
        rng,
    )
    # Convert photon energies from keV to joules.
    energy_J = energy_keV * 1.0e3 * eV_to_J
    # Convert energy to frequency through E = h f.
    frequency_Hz = energy_J / h
    # Convert frequency to wavelength through lambda = c/f.
    wavelength_m = c / frequency_Hz
    # Compute the geometric center of the cell.
    center = 0.5 * cell_size
    # Compute radial vectors from the cell center to every photon.
    radial_vectors = positions - center
    # Compute the magnitude of every radial vector.
    radial_norms = np.linalg.norm(radial_vectors, axis=1)
    # Allocate the photon propagation-direction array.
    directions = np.zeros_like(radial_vectors)
    # Identify photons that are not exactly at the cell center.
    nonzero_radius = radial_norms > 0.0
    # Normalize all nonzero radial vectors.
    directions[nonzero_radius] = (
        radial_vectors[nonzero_radius]
        / radial_norms[nonzero_radius, None]
    )
    # Handle the extremely unlikely case of a photon at the exact cell center.
    if np.any(~nonzero_radius):
        # Give each center-position photon an isotropic propagation direction.
        directions[~nonzero_radius] = sample_isotropic_directions(
            np.count_nonzero(~nonzero_radius),
            rng,
        )
    # Represent photon momentum in keV/c, for which |p| numerically equals E_keV.
    momentum_keVc = directions * energy_keV[:, None]
    # Sum the unweighted photon energy sampled in this cell.
    sampled_energy_J = np.sum(energy_J)
    # Choose one initial macro-weight so the represented energy is exact.
    initial_weight = E_total_cell / sampled_energy_J
    # Assign the same initial weight to every photon.
    weights = np.full(number_of_photons, initial_weight)
    # Return all initialized photon properties.
    return (
        positions,
        energy_keV,
        energy_J,
        frequency_Hz,
        wavelength_m,
        weights,
        momentum_keVc,
    )

# -----------------------------
# Exact moment constraints
# -----------------------------

# Construct the scaled linear constraint matrix used by the exact APM.
def build_exact_constraint_matrix(positions, energy_keV, momentum_keVc):
    # Recover photon unit directions from p/E for on-shell photons.
    directions = momentum_keVc / energy_keV[:, None]
    # Compute a characteristic energy used to nondimensionalize E^2 and E^3.
    mean_energy_keV = np.mean(energy_keV)
    # Create the photon-number row whose weighted sum is total macro-weight.
    number_row = np.ones_like(energy_keV)
    # Create the energy row whose weighted sum is total radiation energy.
    energy_row = energy_keV
    # Use each momentum component as an independent linear constraint row.
    momentum_x_row = momentum_keVc[:, 0]
    # Use the y momentum component as another conserved row.
    momentum_y_row = momentum_keVc[:, 1]
    # Use the z momentum component as another conserved row.
    momentum_z_row = momentum_keVc[:, 2]
    # Preserve the xx component of sum(w E n n).
    pressure_xx_row = energy_keV * directions[:, 0] ** 2
    # Preserve the yy component of sum(w E n n).
    pressure_yy_row = energy_keV * directions[:, 1] ** 2
    # Preserve the xy component of sum(w E n n).
    pressure_xy_row = energy_keV * directions[:, 0] * directions[:, 1]
    # Preserve the xz component of sum(w E n n).
    pressure_xz_row = energy_keV * directions[:, 0] * directions[:, 2]
    # Preserve the yz component of sum(w E n n).
    pressure_yz_row = energy_keV * directions[:, 1] * directions[:, 2]
    # Preserve the number-weighted x spatial first moment.
    centroid_x_row = positions[:, 0]
    # Preserve the number-weighted y spatial first moment.
    centroid_y_row = positions[:, 1]
    # Preserve the number-weighted z spatial first moment.
    centroid_z_row = positions[:, 2]
    # Preserve the second spectral moment in dimensionless form.
    spectral_E2_row = (energy_keV / mean_energy_keV) ** 2
    # Preserve the third spectral moment in dimensionless form.
    spectral_E3_row = (energy_keV / mean_energy_keV) ** 3
    # Stack all fifteen independent moment rows into one matrix.
    raw_matrix = np.vstack(
        (
            number_row,
            energy_row,
            momentum_x_row,
            momentum_y_row,
            momentum_z_row,
            pressure_xx_row,
            pressure_yy_row,
            pressure_xy_row,
            pressure_xz_row,
            pressure_yz_row,
            centroid_x_row,
            centroid_y_row,
            centroid_z_row,
            spectral_E2_row,
            spectral_E3_row,
        )
    )
    # Compute an RMS magnitude for every row to improve conditioning.
    row_scales = np.sqrt(np.mean(raw_matrix * raw_matrix, axis=1))
    # Replace any zero scale by one to avoid division by zero.
    row_scales[row_scales == 0.0] = 1.0
    # Scale every row without changing its null space.
    scaled_matrix = raw_matrix / row_scales[:, None]
    # Return the scaled matrix and row scales for diagnostics if needed.
    return scaled_matrix, row_scales


# Compute candidate weights at both positivity boundaries for one batch.
def compute_batched_weight_candidates(group_weights, null_vectors):
    # Mark positive components of every local null vector.
    positive_components = null_vectors > 1.0e-14
    # Mark negative components of every local null vector.
    negative_components = null_vectors < -1.0e-14
    # Compute positive-direction limiting ratios only where null components are negative.
    positive_ratios = np.where(
        negative_components,
        group_weights / (-null_vectors),
        np.inf,
    )
    # Find the largest allowed positive step for every group.
    positive_steps = np.min(positive_ratios, axis=1)
    # Record which local weight reaches zero on each positive step.
    positive_zero_indices = np.argmin(positive_ratios, axis=1)
    # Compute negative-direction limiting ratios only where null components are positive.
    negative_ratios = np.where(
        positive_components,
        group_weights / null_vectors,
        np.inf,
    )
    # Find the magnitude of the largest allowed negative step for every group.
    negative_steps = np.min(negative_ratios, axis=1)
    # Record which local weight reaches zero on each negative step.
    negative_zero_indices = np.argmin(negative_ratios, axis=1)
    # Advance every group to its positive positivity boundary.
    candidate_plus = group_weights + positive_steps[:, None] * null_vectors
    # Advance every group to its negative positivity boundary.
    candidate_minus = group_weights - negative_steps[:, None] * null_vectors
    # Set the analytically limiting positive-boundary weights to exactly zero.
    candidate_plus[
        np.arange(group_weights.shape[0]),
        positive_zero_indices,
    ] = 0.0
    # Set the analytically limiting negative-boundary weights to exactly zero.
    candidate_minus[
        np.arange(group_weights.shape[0]),
        negative_zero_indices,
    ] = 0.0
    # Compute the squared-weight objective for the positive candidates.
    positive_scores = np.sum(candidate_plus * candidate_plus, axis=1)
    # Compute the squared-weight objective for the negative candidates.
    negative_scores = np.sum(candidate_minus * candidate_minus, axis=1)
    # Prefer the candidate with the smaller sum of squared weights.
    choose_positive = positive_scores <= negative_scores
    # Select one candidate vector for every independent group.
    selected_candidates = np.where(
        choose_positive[:, None],
        candidate_plus,
        candidate_minus,
    )
    # Remove tiny negative roundoff while leaving significant negatives detectable.
    selected_candidates[
        (selected_candidates < 0.0)
        & (selected_candidates > -1.0e-12 * np.max(group_weights))
    ] = 0.0
    # Return the selected nonnegative weights for every group.
    return selected_candidates


# Reduce photons by batched local null-space weight pruning.
def batched_exact_moment_photon_apm(
    positions,
    energy_keV,
    energy_J,
    momentum_keVc,
    weights,
    target_count,
    rng,
):
    """
    Reduce the photon count while preserving fifteen selected moments exactly.

    The algorithm changes only macro-weights.  Surviving photons therefore retain
    their original positions, energies, directions, and on-shell photon momenta.
    For m independent constraints, each local group contains m+1 photons and has
    a one-dimensional null space.  Several disjoint groups are solved together,
    so many photons are removed in one batched NumPy linear-algebra operation.
    """
    # Validate that the requested target count can support all constraints.
    if target_count < 15:
        # Explain why fewer than fifteen particles is not supported here.
        raise ValueError("N_target must be at least the number of constraints.")
    # Validate that the target does not exceed the initial particle count.
    if target_count > weights.size:
        # Reject an invalid request that would require particle creation.
        raise ValueError("N_target cannot exceed the initial particle count.")
    # Store the initial weighted total energy for reporting.
    old_total_energy_J = np.sum(weights * energy_J)
    # Build the scaled linear moment matrix once for the fixed photon states.
    constraint_matrix, _ = build_exact_constraint_matrix(
        positions,
        energy_keV,
        momentum_keVc,
    )
    # Copy the input weights because the APM updates them in place internally.
    current_weights = weights.copy()
    # Mark every original photon as initially active.
    active_indices = np.arange(weights.size, dtype=np.int64)
    # Read the number of independent conserved constraints.
    constraint_count = constraint_matrix.shape[0]
    # Use one more photon than constraints to obtain a one-dimensional null space.
    group_size = constraint_count + 1
    # Continue until exactly the requested number of positive-weight photons remains.
    while active_indices.size > target_count:
        # Compute how many additional photons must be removed.
        removals_needed = active_indices.size - target_count
        # Compute how many disjoint local groups fit in the active population.
        maximum_groups = active_indices.size // group_size
        # Do not process more groups than the number of removals still required.
        groups_this_pass = min(maximum_groups, removals_needed)
        # Protect against an impossible late-stage configuration.
        if groups_this_pass <= 0:
            # Report an informative failure rather than entering an infinite loop.
            raise RuntimeError("Too few active photons remain for another exact group.")
        # Track whether a numerically acceptable batch was found.
        successful_batch = False
        # Retry with a different random grouping if a block is singular or ill-conditioned.
        for _ in range(maximum_batch_retries):
            # Randomize the active photons to avoid persistent grouping bias.
            randomized_active = rng.permutation(active_indices)
            # Take exactly enough indices to fill all groups in this pass.
            grouped_indices = randomized_active[
                : groups_this_pass * group_size
            ].reshape(groups_this_pass, group_size)
            # Extract each local 15-by-16 constraint block.
            local_blocks = np.transpose(
                constraint_matrix[:, grouped_indices],
                (1, 0, 2),
            )
            # Use the first fifteen columns as a square batched coefficient matrix.
            square_blocks = local_blocks[:, :, :constraint_count]
            # Move the last column to the right-hand side of A c = -a_last.
            right_hand_sides = -local_blocks[:, :, constraint_count]
            # Attempt to solve all independent 15-by-15 systems in one NumPy call.
            try:
                # Compute the first fifteen components of every null vector.
                dependent_coefficients = np.linalg.solve(
                    square_blocks,
                    right_hand_sides[:, :, None],
                )[:, :, 0]
            # Catch a singular block and retry with different random groups.
            except np.linalg.LinAlgError:
                # Skip directly to the next random grouping attempt.
                continue
            # Append a final component of one to complete each null vector.
            null_vectors = np.concatenate(
                (
                    dependent_coefficients,
                    np.ones((groups_this_pass, 1)),
                ),
                axis=1,
            )
            # Compute the maximum absolute component of every null vector.
            null_scales = np.max(np.abs(null_vectors), axis=1)
            # Reject nonfinite or identically zero null vectors.
            if np.any(~np.isfinite(null_scales)) or np.any(null_scales == 0.0):
                # Retry with a fresh random grouping.
                continue
            # Normalize every null vector to improve numerical consistency.
            null_vectors /= null_scales[:, None]
            # Compute all local null residuals in one batched matrix product.
            null_residuals = np.linalg.norm(
                np.einsum("gij,gj->gi", local_blocks, null_vectors),
                axis=1,
            )
            # Compute a scale for each local constraint block.
            local_norms = np.linalg.norm(local_blocks, axis=(1, 2))
            # Reject the batch if any local null vector is insufficiently accurate.
            if np.any(
                null_residuals
                > relative_null_tolerance * np.maximum(local_norms, 1.0)
            ):
                # Retry all groups with a different random assignment.
                continue
            # Extract the current weights belonging to every group.
            group_weights = current_weights[grouped_indices]
            # Require every null vector to contain both positive and negative entries.
            has_positive = np.any(null_vectors > 1.0e-14, axis=1)
            # Check separately for at least one negative null component.
            has_negative = np.any(null_vectors < -1.0e-14, axis=1)
            # Reject a pathological batch that cannot reach a positivity boundary.
            if not np.all(has_positive & has_negative):
                # Retry with a new grouping.
                continue
            # Compute the lower-weight-variance boundary candidate for every group.
            selected_group_weights = compute_batched_weight_candidates(
                group_weights,
                null_vectors,
            )
            # Reject significant negative weights caused by numerical failure.
            if np.any(selected_group_weights < -1.0e-11 * np.max(group_weights)):
                # Retry instead of silently violating nonnegativity.
                continue
            # Clip any remaining tiny negative roundoff to zero.
            selected_group_weights[selected_group_weights < 0.0] = 0.0
            # Compute a scale-dependent zero threshold for the current pass.
            zero_tolerance = (
                relative_zero_tolerance
                * max(np.max(current_weights[active_indices]), 1.0)
            )
            # Count how many photons each group would remove at this tolerance.
            removed_per_group = np.count_nonzero(
                selected_group_weights <= zero_tolerance,
                axis=1,
            )
            # Prefer the normal case of exactly one removed photon per group.
            if not np.all(removed_per_group == 1):
                # Retry because multiple simultaneous zeros could overshoot the target.
                continue
            # Store all selected group weights back into the global weight array.
            current_weights[grouped_indices] = selected_group_weights
            # Rebuild the active-index list from strictly positive weights.
            active_indices = active_indices[
                current_weights[active_indices] > zero_tolerance
            ]
            # Confirm that the batched pass did not cross below the target count.
            if active_indices.size < target_count:
                # This should be prevented by the one-zero-per-group test above.
                raise RuntimeError("Batched pruning overshot the target count.")
            # Mark this batched pass as successful.
            successful_batch = True
            # Leave the retry loop and begin the next reduction pass.
            break
        # Stop with a clear message if no acceptable random grouping was found.
        if not successful_batch:
            # Include the active population to aid debugging difficult distributions.
            raise RuntimeError(
                "Batched exact-moment pruning failed with "
                f"{active_indices.size} active photons."
            )
    # Copy the surviving photon positions without altering their coordinates.
    final_positions = positions[active_indices].copy()
    # Copy the surviving photon energies in keV without spectral averaging.
    final_energy_keV = energy_keV[active_indices].copy()
    # Copy the surviving photon energies in joules.
    final_energy_J = energy_J[active_indices].copy()
    # Copy the original on-shell momenta of the surviving photons.
    final_momentum_keVc = momentum_keVc[active_indices].copy()
    # Copy the final nonuniform weights of the surviving photons.
    final_weights = current_weights[active_indices].copy()
    # Compute the final weighted total radiation energy in joules.
    new_total_energy_J = np.sum(final_weights * final_energy_J)
    # Return the reduced photons and before-and-after total energies.
    return (
        final_positions,
        final_energy_keV,
        final_energy_J,
        final_momentum_keVc,
        final_weights,
        old_total_energy_J,
        new_total_energy_J,
    )

# -----------------------------
# Plotting utilities
# -----------------------------

# Plot photon positions in the x-y plane and color them by photon energy.
def plot_positions(positions, energy_keV, weights, title_suffix):
    # Create one figure for the x-y position distribution.
    plt.figure(figsize=(6, 5))
    # Scale marker area weakly with macro-weight for the nonuniform final weights.
    marker_sizes = 12.0 * weights / np.mean(weights)
    # Draw the photon positions and color every marker by energy.
    scatter = plt.scatter(
        positions[:, 0],
        positions[:, 1],
        s=marker_sizes,
        c=energy_keV,
        alpha=0.75,
    )
    # Add a colorbar identifying the photon energy.
    plt.colorbar(scatter, label="Photon energy [keV]")
    # Label the horizontal coordinate.
    plt.xlabel("x (cell units)")
    # Label the vertical coordinate.
    plt.ylabel("y (cell units)")
    # Add the selected before-or-after description to the title.
    plt.title(f"Photon positions in cell {title_suffix}")
    # Reduce excess whitespace around the plot.
    plt.tight_layout()


# Plot the weighted photon-energy distribution and the target Planck curve.
def plot_weighted_energy_histogram(
    energy_keV,
    weights,
    thermal_energy_keV,
    title_suffix,
):
    # Create one figure for the weighted spectral distribution.
    plt.figure(figsize=(7, 5))
    # Use fifty linearly spaced histogram bins as in the earlier scripts.
    number_of_bins = 50
    # Set the lower plotted energy bound to zero.
    minimum_energy = 0.0
    # Extend the upper plotted bound beyond the largest sampled photon energy.
    maximum_energy = 1.1 * np.max(energy_keV)
    # Construct the histogram bin edges.
    bins = np.linspace(minimum_energy, maximum_energy, number_of_bins)
    # Normalize macro-weights so the histogram integrates to one.
    normalized_weights = weights / np.sum(weights)
    # Plot the weighted photon-number energy distribution.
    plt.hist(
        energy_keV,
        bins=bins,
        weights=normalized_weights,
        density=True,
        alpha=0.6,
        label="Weighted photon distribution",
    )
    # Construct a smooth energy grid for the target Planck curve.
    plotted_energy = np.linspace(bins[0] + 1.0e-6, bins[-1], 500)
    # Convert energy to the dimensionless variable x = E/(kT).
    dimensionless_energy = plotted_energy / thermal_energy_keV
    # Evaluate the unnormalized Planck photon-number probability density.
    target_pdf = dimensionless_energy**2 / np.expm1(dimensionless_energy)
    # Integrate the target density numerically with the trapezoidal rule.
    target_normalization = np.sum(
        0.5 * (target_pdf[1:] + target_pdf[:-1])
        * np.diff(plotted_energy)
    )
    # Normalize the target probability density over the displayed range.
    target_pdf /= target_normalization
    # Overlay the target Planck photon-number distribution.
    plt.plot(
        plotted_energy,
        target_pdf,
        linewidth=2,
        label="Target Planck PDF",
    )
    # Label the horizontal energy axis.
    plt.xlabel("Photon energy [keV]")
    # Label the vertical probability-density axis.
    plt.ylabel("Probability density")
    # Add the before-or-after method description to the title.
    plt.title(f"Photon energy distribution {title_suffix}")
    # Display the histogram and target-curve labels.
    plt.legend()
    # Reduce excess whitespace around the plot.
    plt.tight_layout()


# Plot all nine position-momentum correlations in a 3-by-3 grid.
def plot_position_momentum_correlations(
    positions,
    momentum_keVc,
    weights,
    title_suffix,
):
    # Separate the three spatial coordinates for convenient plotting.
    coordinates = [positions[:, 0], positions[:, 1], positions[:, 2]]
    # Separate the three photon momentum components.
    momenta = [
        momentum_keVc[:, 0],
        momentum_keVc[:, 1],
        momentum_keVc[:, 2],
    ]
    # Define readable labels for the spatial coordinates.
    coordinate_labels = ["x", "y", "z"]
    # Define readable labels for the momentum components.
    momentum_labels = ["p_x", "p_y", "p_z"]
    # Create a 3-by-3 array of plotting axes.
    figure, axes = plt.subplots(3, 3, figsize=(12, 10))
    # Add a title that spans all nine panels.
    figure.suptitle(
        f"Position-momentum correlations {title_suffix} (p in keV/c)",
        fontsize=14,
    )
    # Scale marker area weakly with the nonuniform macro-weights.
    marker_sizes = 6.0 * weights / np.mean(weights)
    # Loop over spatial coordinates for the subplot rows.
    for coordinate_index in range(3):
        # Loop over momentum components for the subplot columns.
        for momentum_index in range(3):
            # Select the current subplot axis.
            axis = axes[coordinate_index, momentum_index]
            # Draw the weighted photon phase-space samples.
            axis.scatter(
                coordinates[coordinate_index],
                momenta[momentum_index],
                s=marker_sizes,
                alpha=0.35,
            )
            # Label the spatial-coordinate axis.
            axis.set_xlabel(
                f"{coordinate_labels[coordinate_index]} (cell units)"
            )
            # Label the momentum-component axis.
            axis.set_ylabel(
                f"{momentum_labels[momentum_index]} [keV/c]"
            )
    # Reserve space for the spanning title while tightening the panels.
    plt.tight_layout(rect=[0.0, 0.03, 1.0, 0.95])

# -----------------------------
# Diagnostic utilities
# -----------------------------

# Compute a weighted mean and population standard deviation.
def weighted_mean_and_std(values, weights):
    # Sum the macro-weights.
    total_weight = np.sum(weights)
    # Compute the weighted arithmetic mean.
    mean = np.sum(weights * values) / total_weight
    # Compute the weighted population variance.
    variance = np.sum(weights * (values - mean) ** 2) / total_weight
    # Return the weighted mean and standard deviation.
    return mean, np.sqrt(variance)


# Compute weighted L1 and L2 distances between two energy histograms.
def weighted_histogram_distance(
    old_energy_keV,
    old_weights,
    new_energy_keV,
    new_weights,
    number_of_bins=60,
):
    # Set a common upper energy bound for both histograms.
    maximum_energy = 1.05 * max(
        np.max(old_energy_keV),
        np.max(new_energy_keV),
    )
    # Create common histogram bin edges.
    bins = np.linspace(0.0, maximum_energy, number_of_bins + 1)
    # Form the old weighted histogram.
    old_histogram = np.histogram(
        old_energy_keV,
        bins=bins,
        weights=old_weights,
    )[0]
    # Form the new weighted histogram.
    new_histogram = np.histogram(
        new_energy_keV,
        bins=bins,
        weights=new_weights,
    )[0]
    # Normalize the old histogram to unit total weight.
    old_histogram /= np.sum(old_histogram)
    # Normalize the new histogram to unit total weight.
    new_histogram /= np.sum(new_histogram)
    # Compute the bin-by-bin probability difference.
    difference = new_histogram - old_histogram
    # Compute the discrete L1 histogram distance.
    L1_distance = np.sum(np.abs(difference))
    # Compute the discrete L2 histogram distance.
    L2_distance = np.sqrt(np.sum(difference**2))
    # Return both spectral-shape diagnostics.
    return L1_distance, L2_distance


# Compute the photon radiation pressure tensor sum(w E n_a n_b).
def photon_pressure_tensor(energy_keV, momentum_keVc, weights):
    # Recover the photon propagation directions from p/E.
    directions = momentum_keVc / energy_keV[:, None]
    # Contract the weighted direction dyads into a 3-by-3 tensor.
    pressure_tensor = np.einsum(
        "i,ij,ik->jk",
        weights * energy_keV,
        directions,
        directions,
    )
    # Return the pressure tensor in consistent code units.
    return pressure_tensor


# Compute a normwise relative error with a safe denominator.
def relative_norm_error(new_value, old_value):
    # Compute the Euclidean or Frobenius norm of the reference quantity.
    denominator = np.linalg.norm(old_value)
    # Compute and return the normalized difference.
    return np.linalg.norm(new_value - old_value) / max(denominator, 1.0e-300)


# Print the same basic before-and-after metrics plus exact-moment diagnostics.
def report_error_metrics(
    positions_old,
    positions_new,
    momentum_old,
    momentum_new,
    energy_old_J,
    energy_new_J,
    energy_old_keV,
    energy_new_keV,
    weights_old,
    weights_new,
    total_energy_old,
    total_energy_new,
    elapsed_time_ms,
):
    # Print a descriptive header for the APM diagnostics.
    print("\n=== Batched Exact-Moment Photon APM Error Metrics ===")
    # Print the initial photon count.
    print(f"Initial particle count: {positions_old.shape[0]}")
    # Print the final photon count.
    print(f"Final   particle count: {positions_new.shape[0]}")
    # Print the measured APM execution time for this call.
    print(f"APM execution time: {elapsed_time_ms:.8f} ms")
    # Compute the relative total-energy error.
    relative_energy_error = (
        total_energy_new - total_energy_old
    ) / total_energy_old
    # Print the initial weighted energy.
    print(f"Total energy old: {total_energy_old:.8e} J")
    # Print the final weighted energy.
    print(f"Total energy new: {total_energy_new:.8e} J")
    # Print the relative energy-conservation error.
    print(f"Relative total-energy error: {relative_energy_error:.8e}")
    # Compute the old weighted mean and standard deviation of photon energy.
    mean_energy_old, std_energy_old = weighted_mean_and_std(
        energy_old_J,
        weights_old,
    )
    # Compute the new weighted mean and standard deviation of photon energy.
    mean_energy_new, std_energy_new = weighted_mean_and_std(
        energy_new_J,
        weights_new,
    )
    # Print the old mean energy in keV.
    print(
        "Weighted mean energy old: "
        f"{mean_energy_old / eV_to_J / 1.0e3:.8f} keV"
    )
    # Print the new mean energy in keV.
    print(
        "Weighted mean energy new: "
        f"{mean_energy_new / eV_to_J / 1.0e3:.8f} keV"
    )
    # Print the old weighted spectral width in keV.
    print(
        "Weighted std  energy old: "
        f"{std_energy_old / eV_to_J / 1.0e3:.8f} keV"
    )
    # Print the new weighted spectral width in keV.
    print(
        "Weighted std  energy new: "
        f"{std_energy_new / eV_to_J / 1.0e3:.8f} keV"
    )
    # Loop over all three spatial dimensions.
    for component, label in enumerate(("x", "y", "z")):
        # Compute the old weighted spatial mean and standard deviation.
        mean_old, std_old = weighted_mean_and_std(
            positions_old[:, component],
            weights_old,
        )
        # Compute the new weighted spatial mean and standard deviation.
        mean_new, std_new = weighted_mean_and_std(
            positions_new[:, component],
            weights_new,
        )
        # Print the old and new weighted spatial means.
        print(f"{label}-pos mean old/new: {mean_old:.8f}, {mean_new:.8f}")
        # Print the old and new weighted spatial standard deviations.
        print(f"{label}-pos std  old/new: {std_old:.8f}, {std_new:.8f}")
    # Loop over all three momentum components.
    for component, label in enumerate(("p_x", "p_y", "p_z")):
        # Compute the old weighted momentum mean and standard deviation.
        mean_old, std_old = weighted_mean_and_std(
            momentum_old[:, component],
            weights_old,
        )
        # Compute the new weighted momentum mean and standard deviation.
        mean_new, std_new = weighted_mean_and_std(
            momentum_new[:, component],
            weights_new,
        )
        # Print the old and new weighted momentum means.
        print(
            f"{label} mean old/new: {mean_old:.8f}, "
            f"{mean_new:.8f} [keV/c]"
        )
        # Print the old and new weighted momentum standard deviations.
        print(
            f"{label} std  old/new: {std_old:.8f}, "
            f"{std_new:.8f} [keV/c]"
        )
    # Compute weighted energy-histogram distances.
    histogram_L1, histogram_L2 = weighted_histogram_distance(
        energy_old_keV,
        weights_old,
        energy_new_keV,
        weights_new,
    )
    # Print the weighted L1 spectral-shape error.
    print(f"Weighted energy histogram L1 distance: {histogram_L1:.8e}")
    # Print the weighted L2 spectral-shape error.
    print(f"Weighted energy histogram L2 distance: {histogram_L2:.8e}")
    # Compute the old total momentum vector.
    total_momentum_old = np.sum(
        weights_old[:, None] * momentum_old,
        axis=0,
    )
    # Compute the new total momentum vector.
    total_momentum_new = np.sum(
        weights_new[:, None] * momentum_new,
        axis=0,
    )
    # Compute the old photon pressure tensor.
    pressure_old = photon_pressure_tensor(
        energy_old_keV,
        momentum_old,
        weights_old,
    )
    # Compute the new photon pressure tensor.
    pressure_new = photon_pressure_tensor(
        energy_new_keV,
        momentum_new,
        weights_new,
    )
    # Compute the old number-weighted spatial centroid.
    centroid_old = np.sum(
        weights_old[:, None] * positions_old,
        axis=0,
    ) / np.sum(weights_old)
    # Compute the new number-weighted spatial centroid.
    centroid_new = np.sum(
        weights_new[:, None] * positions_new,
        axis=0,
    ) / np.sum(weights_new)
    # Compute the relative total-weight error.
    total_weight_error = abs(
        np.sum(weights_new) - np.sum(weights_old)
    ) / np.sum(weights_old)
    # Compute the relative total-momentum error.
    momentum_error = relative_norm_error(
        total_momentum_new,
        total_momentum_old,
    )
    # Compute the relative pressure-tensor error.
    pressure_error = relative_norm_error(pressure_new, pressure_old)
    # Compute the relative centroid error.
    centroid_error = relative_norm_error(centroid_new, centroid_old)
    # Compute the old second weighted spectral moment.
    spectral_E2_old = np.sum(weights_old * energy_old_keV**2)
    # Compute the new second weighted spectral moment.
    spectral_E2_new = np.sum(weights_new * energy_new_keV**2)
    # Compute the old third weighted spectral moment.
    spectral_E3_old = np.sum(weights_old * energy_old_keV**3)
    # Compute the new third weighted spectral moment.
    spectral_E3_new = np.sum(weights_new * energy_new_keV**3)
    # Compute the relative E^2-moment error.
    spectral_E2_error = abs(
        spectral_E2_new - spectral_E2_old
    ) / spectral_E2_old
    # Compute the relative E^3-moment error.
    spectral_E3_error = abs(
        spectral_E3_new - spectral_E3_old
    ) / spectral_E3_old
    # Compute the effective sample size implied by nonuniform weights.
    effective_sample_size = (
        np.sum(weights_new) ** 2 / np.sum(weights_new**2)
    )
    # Compute the ratio of the largest output weight to the mean output weight.
    maximum_to_mean_weight = np.max(weights_new) / np.mean(weights_new)
    # Print the total-weight conservation error.
    print(f"Relative total-weight error: {total_weight_error:.8e}")
    # Print the total-vector-momentum conservation error.
    print(f"Relative total-momentum error: {momentum_error:.8e}")
    # Print the pressure-tensor conservation error.
    print(f"Relative pressure-tensor error: {pressure_error:.8e}")
    # Print the spatial-centroid conservation error.
    print(f"Relative spatial-centroid error: {centroid_error:.8e}")
    # Print the second spectral-moment conservation error.
    print(f"Relative E^2 spectral-moment error: {spectral_E2_error:.8e}")
    # Print the third spectral-moment conservation error.
    print(f"Relative E^3 spectral-moment error: {spectral_E3_error:.8e}")
    # Print the maximum-to-mean weight ratio.
    print(f"Maximum/mean output weight: {maximum_to_mean_weight:.8f}")
    # Print the effective sample size of the nonuniform output weights.
    print(f"Effective sample size: {effective_sample_size:.8f}")

# -----------------------------
# Main example
# -----------------------------

# Run the complete initialization, APM reduction, diagnostics, and plotting example.
def main():
    # Create one reproducible NumPy random-number generator.
    rng = np.random.default_rng(random_seed)
    # Initialize the original 1 keV blackbody photon population.
    (
        positions,
        energy_keV,
        energy_J,
        frequency_Hz,
        wavelength_m,
        weights,
        momentum_keVc,
    ) = initialize_photons_in_cell(N_particles, kT_keV, rng)
    # Print the equivalent blackbody temperature.
    print(f"Temperature T = {T_K:.8e} K")
    # Print the blackbody radiation energy density.
    print(f"Radiation energy density u_rad = {u_rad:.8e} J/m^3")
    # Print the target radiation energy represented in the cell.
    print(f"Total radiation energy in cell (target): {E_total_cell:.8e} J")
    # Print the initial uniform photon macro-weight.
    print(f"Weight per macroparticle (original): {weights[0]:.8e}")
    # Verify the initialized weighted radiation energy.
    print(f"Check: sum(w_i * E_i) = {np.sum(weights * energy_J):.8e} J")
    # Plot the original photon positions.
    plot_positions(
        positions,
        energy_keV,
        weights,
        title_suffix="(original, Fig. 1)",
    )
    # Plot the original weighted energy distribution.
    plot_weighted_energy_histogram(
        energy_keV,
        weights,
        kT_keV,
        title_suffix="(original, Fig. 2)",
    )
    # Plot the original position-momentum correlations.
    plot_position_momentum_correlations(
        positions,
        momentum_keVc,
        weights,
        title_suffix="(original, Fig. 3)",
    )
    # Start the timer immediately before the APM routine.
    start_time = time.perf_counter_ns()
    # Apply the accelerated batched exact-moment APM.
    (
        positions_new,
        energy_keV_new,
        energy_J_new,
        momentum_keVc_new,
        weights_new,
        total_energy_old,
        total_energy_new,
    ) = batched_exact_moment_photon_apm(
        positions,
        energy_keV,
        energy_J,
        momentum_keVc,
        weights,
        N_target,
        rng,
    )
    # Stop the timer immediately after the APM routine.
    elapsed_time_ms = (time.perf_counter_ns() - start_time) * 1.0e-6
    # Plot the reduced photon positions with marker size indicating weight.
    plot_positions(
        positions_new,
        energy_keV_new,
        weights_new,
        title_suffix="(batched exact moments, Fig. 4)",
    )
    # Plot the reduced weighted energy distribution.
    plot_weighted_energy_histogram(
        energy_keV_new,
        weights_new,
        kT_keV,
        title_suffix="(batched exact moments, Fig. 5)",
    )
    # Plot the reduced position-momentum correlations.
    plot_position_momentum_correlations(
        positions_new,
        momentum_keVc_new,
        weights_new,
        title_suffix="(batched exact moments, Fig. 6)",
    )
    # Print all before-and-after and exact-moment diagnostics.
    report_error_metrics(
        positions,
        positions_new,
        momentum_keVc,
        momentum_keVc_new,
        energy_J,
        energy_J_new,
        energy_keV,
        energy_keV_new,
        weights,
        weights_new,
        total_energy_old,
        total_energy_new,
        elapsed_time_ms,
    )
    # Display all six figures after the terminal diagnostics are complete.
    plt.show()


# Execute the example only when this file is run as a script.
if __name__ == "__main__":
    # Call the main driver function.
    main()
