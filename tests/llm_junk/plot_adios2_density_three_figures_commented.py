#!/usr/bin/env python3
"""Create three density figures from ion and electron ADIOS2 BP data.

This script produces exactly three PNG files:

1. ``density_2x2_zeros_included_whitebelowvmin.png``
   A 2x2 comparison containing the native ion density, native electron
   density, smoothed ion density, and smoothed electron density.

2. ``density_smoothed_combined_half_figure_square.png``
   A square composite containing only the smoothed ion data for Z < 0
   and only the smoothed electron data for Z >= 0.

3. ``ion_density_x_slice_avg6cells_about_z0.png``
   A one-dimensional native ion-density lineout through X, averaged over
   the six native Z cells closest to Z = 0.

The script accepts either conventional ADIOS2 BP directories or the split
component-file prefixes used by the uploaded data. For example, ``--ion i``
refers to files named ``i_md.idx``, ``i_md.0``, ``i_data.0``, and ``i_mmd.0``.

Required packages:

    python -m pip install adios2 numpy matplotlib

Example command:

    python plot_adios2_density_commented.py \
        --ion /path/to/i \
        --electron /path/to/e \
        --output-dir /path/to/plots
"""

# Enable modern type annotations while retaining compatibility with older Python versions.
from __future__ import annotations

# Import argparse so command-line input paths and plotting options can be specified clearly.
import argparse

# Import contextlib so a temporary BP directory can be managed with a context manager.
import contextlib

# Import os so symbolic links or hard links can be created for split BP component files.
import os

# Import Path for reliable platform-independent filesystem path handling.
from pathlib import Path

# Import shutil as a final fallback for copying files when links cannot be created.
import shutil

# Import tempfile so split BP components can be assembled in a temporary directory.
import tempfile

# Import Iterator for the return type of the BP-directory context manager.
from typing import Iterator

# Import the ADIOS2 Python reader used to load the density arrays and attributes.
import adios2

# Import Matplotlib's plotting interface for all figure construction and saving.
import matplotlib.pyplot as plt

# Import logarithmic normalization for the density color scales.
from matplotlib.colors import LogNorm

# Import NumPy for array manipulation, smoothing, masking, and coordinate construction.
import numpy as np


# List the files that can be present inside one ADIOS2 BP directory.
BP_COMPONENTS = ("md.idx", "md.0", "data.0", "mmd.0", "profiling.json")

# Set the default lower limit of both logarithmic density color scales.
DEFAULT_VMIN = 1.0e24

# Set the default upper limit of both logarithmic density color scales.
DEFAULT_VMAX = 1.0e28

# Use a 5-by-5 boxcar by default, matching the previous plots in this workflow.
DEFAULT_SMOOTH_WINDOW = 5

# Apply the boxcar twice by default to remove isolated holes more completely.
DEFAULT_SMOOTH_PASSES = 2

# Convert meters to micrometers whenever physical axes are displayed.
METERS_TO_MICROMETERS = 1.0e6


@contextlib.contextmanager
def bp_dataset(source: str | Path) -> Iterator[Path]:
    """Yield a readable ADIOS2 BP directory for a directory or split-file prefix."""

    # Convert the supplied path to an absolute Path object.
    source_path = Path(source).expanduser().resolve()

    # Use an existing BP directory directly when one was supplied.
    if source_path.is_dir():
        # Yield the directory to the caller without creating temporary files.
        yield source_path

        # Stop after the caller finishes using the directory.
        return

    # Build the expected paths for the four component files required to read the BP data.
    required_paths = [
        source_path.parent / f"{source_path.name}_{component}"
        for component in BP_COMPONENTS[:4]
    ]

    # Identify any required component files that are missing.
    missing_paths = [path for path in required_paths if not path.exists()]

    # Report a precise error if the supplied prefix does not resolve to a complete dataset.
    if missing_paths:
        # Format the missing paths on separate indented lines for readability.
        missing_text = "\n  ".join(str(path) for path in missing_paths)

        # Raise an error rather than allowing ADIOS2 to fail with a less useful message.
        raise FileNotFoundError(
            f"Could not find a BP directory at {source_path}, and these "
            f"required split component files are missing:\n  {missing_text}"
        )

    # Create a temporary parent directory that is removed automatically afterward.
    with tempfile.TemporaryDirectory(prefix=f"{source_path.name}_bp_") as temporary_root:
        # Construct the BP directory name expected by ADIOS2 inside the temporary directory.
        bp_directory = Path(temporary_root) / f"{source_path.name}.bp"

        # Create the temporary BP directory itself.
        bp_directory.mkdir()

        # Process each possible BP component, including the optional profiling file.
        for component in BP_COMPONENTS:
            # Locate the uploaded split component file beside the supplied prefix.
            original_path = source_path.parent / f"{source_path.name}_{component}"

            # Skip optional files that are not present.
            if not original_path.exists():
                # Continue to the next component without treating the optional file as an error.
                continue

            # Construct the normal component filename inside the temporary BP directory.
            target_path = bp_directory / component

            # Prefer a symbolic link so the large data file is not duplicated.
            try:
                # Create a symbolic link from the temporary BP directory to the original file.
                os.symlink(original_path, target_path)

            # Fall back when symbolic links are unavailable or disallowed.
            except OSError:
                # Prefer a hard link as the next most efficient option.
                try:
                    # Create a hard link without copying the file contents.
                    os.link(original_path, target_path)

                # Fall back to an actual file copy only when neither link type is possible.
                except OSError:
                    # Copy the component while preserving its metadata.
                    shutil.copy2(original_path, target_path)

        # Yield the assembled BP directory while the temporary directory still exists.
        yield bp_directory


def read_attribute(reader: adios2.FileReader, name: str) -> np.ndarray:
    """Read one ADIOS2 attribute and return it as a NumPy array."""

    # Ask ADIOS2 for the named attribute and normalize its result to a NumPy array.
    return np.asarray(reader.read_attribute(name))


def read_density(
    source: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, int]:
    """Read the X-Z density plane, coordinate limits, simulation time, and step."""

    # Resolve a conventional BP directory regardless of how the input was supplied.
    with bp_dataset(source) as bp_directory:
        # Open the BP dataset with the high-level ADIOS2 file reader.
        with adios2.FileReader(str(bp_directory)) as reader:
            # Retrieve the variable metadata so the expected field can be validated.
            available_variables = reader.available_variables()

            # Stop with a clear message if the dataset does not contain Density.
            if "Density" not in available_variables:
                # Raise a keyed error because the requested field is absent.
                raise KeyError(f"No 'Density' variable exists in {bp_directory}")

            # Read the full three-dimensional density field as double precision.
            density_xyz = np.asarray(reader.read("Density"), dtype=np.float64)

            # Read the physical X-axis limits stored as an ADIOS2 attribute.
            x_range_m = read_attribute(reader, "x_range").astype(np.float64)

            # Read the physical Z-axis limits stored as an ADIOS2 attribute.
            z_range_m = read_attribute(reader, "z_range").astype(np.float64)

            # Read and reduce the simulation time to one Python floating-point value.
            simulation_time_s = float(np.asarray(reader.read("Time")).squeeze())

            # Read and reduce the simulation step to one Python integer value.
            simulation_step = int(np.asarray(reader.read("Step")).squeeze())

    # Confirm that the stored density has the expected X-Y-Z dimensionality.
    if density_xyz.ndim != 3:
        # Report the actual shape to make incompatible data easy to diagnose.
        raise ValueError(
            f"Expected Density(X,Y,Z) to be three-dimensional; found {density_xyz.shape}."
        )

    # Confirm that this dataset contains exactly one Y plane.
    if density_xyz.shape[1] != 1:
        # Stop rather than silently selecting one plane from a genuinely three-dimensional field.
        raise ValueError(
            "This script expects Density(X,Y,Z) with one Y cell; "
            f"found shape {density_xyz.shape}."
        )

    # Remove the singleton Y dimension to obtain a two-dimensional X-Z array.
    density_xz = density_xyz[:, 0, :]

    # Return all data needed for plotting and labeling.
    return density_xz, x_range_m, z_range_m, simulation_time_s, simulation_step


def box_sum_2d(array: np.ndarray, window: int) -> np.ndarray:
    """Return the centered two-dimensional box sum using a fast integral image."""

    # Compute the number of cells required on each side of an odd-width window.
    radius = window // 2

    # Pad the array with zeros so the output retains exactly the original shape.
    padded = np.pad(
        array,
        ((radius, radius), (radius, radius)),
        mode="constant",
        constant_values=0.0,
    )

    # Add one extra zero row and column to simplify integral-image indexing.
    integral = np.pad(
        padded,
        ((1, 0), (1, 0)),
        mode="constant",
        constant_values=0.0,
    )

    # Form the cumulative sum first along X and then along Z.
    integral = integral.cumsum(axis=0).cumsum(axis=1)

    # Use four integral-image corners to obtain every centered window sum at once.
    summed = (
        integral[window:, window:]
        - integral[:-window, window:]
        - integral[window:, :-window]
        + integral[:-window, :-window]
    )

    # Return an array with the same shape as the unpadded input.
    return summed


def boxcar_smooth_include_zeros(
    density_xz: np.ndarray,
    window: int,
    passes: int,
) -> np.ndarray:
    """Apply a true boxcar average in which zero-valued cells remain in the mean."""

    # Require a positive odd window so every output cell has a centered neighborhood.
    if window < 1 or window % 2 == 0:
        # Reject invalid window sizes with a direct explanation.
        raise ValueError("The smoothing window must be a positive odd integer.")

    # Require at least one pass so the operation is meaningful.
    if passes < 1:
        # Reject zero or negative pass counts explicitly.
        raise ValueError("The number of smoothing passes must be at least one.")

    # Copy the input so the original native density array is never modified.
    smoothed = np.array(density_xz, dtype=np.float64, copy=True)

    # Count every cell in the boxcar denominator, including zero-valued holes.
    window_area = float(window * window)

    # Apply the requested number of identical boxcar passes.
    for _pass_index in range(passes):
        # Replace negative numerical noise with zero before performing the physical average.
        nonnegative = np.where(smoothed > 0.0, smoothed, 0.0)

        # Compute the sum over every centered boxcar window efficiently.
        neighborhood_sum = box_sum_2d(nonnegative, window)

        # Divide by the full number of cells, so zero-valued holes dilute the average.
        smoothed = neighborhood_sum / window_area

    # Return the full-resolution smoothed density field.
    return smoothed


def coordinate_centers(range_m: np.ndarray, number_of_cells: int) -> np.ndarray:
    """Construct cell-center coordinates from a stored physical axis range."""

    # Linearly distribute coordinates across the stored endpoints.
    return np.linspace(
        float(range_m[0]),
        float(range_m[1]),
        number_of_cells,
        dtype=np.float64,
    )


def mask_below_display_minimum(
    density_xz: np.ndarray,
    vmin: float,
) -> np.ma.MaskedArray:
    """Mask all values below vmin so zero and low nonzero values display as white."""

    # Copy the input to prevent later plotting operations from changing the source array.
    display_density = np.array(density_xz, dtype=np.float64, copy=True)

    # Build a mask for zeros, negative values, and positive values below the color scale.
    below_minimum = display_density < vmin

    # Return a masked array whose masked pixels will use the colormap's bad color.
    return np.ma.masked_array(display_density, mask=below_minimum)


def make_density_colormap(name: str):
    """Return a private colormap copy with masked values rendered in white."""

    # Copy the requested Matplotlib colormap so global colormap state is not changed.
    colormap = plt.get_cmap(name).copy()

    # Assign pure white to masked values, including zeros and values below vmin.
    colormap.set_bad("white")

    # Return the configured colormap.
    return colormap


def configure_density_axes(ax: plt.Axes) -> None:
    """Apply the common Z-horizontal and flipped-X-vertical axis convention."""

    # Label the horizontal coordinate as Z in micrometers.
    ax.set_xlabel(r"$Z\;(\mu\mathrm{m})$")

    # Label the vertical coordinate as X in micrometers.
    ax.set_ylabel(r"$X\;(\mu\mathrm{m})$")

    # Display Z conventionally from -15 micrometers at left to +15 at right.
    ax.set_xlim(-15.0, 15.0)

    # Display X with +15 at the bottom and -15 at the top, as requested.
    ax.set_ylim(15.0, -15.0)

    # Use the same seven major Z ticks on every panel.
    ax.set_xticks([-15, -10, -5, 0, 5, 10, 15])

    # Use the same seven major X ticks on every panel.
    ax.set_yticks([15, 10, 5, 0, -5, -10, -15])

    # Preserve one physical unit in Z per physical unit in X.
    ax.set_aspect("equal", adjustable="box")


def draw_density(
    ax: plt.Axes,
    density_xz: np.ndarray,
    x_range_m: np.ndarray,
    z_range_m: np.ndarray,
    title: str,
    colormap_name: str,
    vmin: float,
    vmax: float,
):
    """Draw one density field and return the image object used by a colorbar."""

    # Mask all values below vmin so they render as white rather than the vmin color.
    plotted_density = mask_below_display_minimum(density_xz, vmin)

    # Create the requested colormap and configure masked values to be white.
    colormap = make_density_colormap(colormap_name)

    # Convert the stored Z limits from meters to micrometers.
    z_min_um = float(z_range_m[0]) * METERS_TO_MICROMETERS

    # Convert the stored upper Z limit from meters to micrometers.
    z_max_um = float(z_range_m[1]) * METERS_TO_MICROMETERS

    # Convert the stored X limits from meters to micrometers.
    x_min_um = float(x_range_m[0]) * METERS_TO_MICROMETERS

    # Convert the stored upper X limit from meters to micrometers.
    x_max_um = float(x_range_m[1]) * METERS_TO_MICROMETERS

    # Draw the X-Z array with Z horizontal and X vertical.
    image = ax.imshow(
        plotted_density,
        origin="lower",
        extent=[z_min_um, z_max_um, x_min_um, x_max_um],
        interpolation="nearest",
        aspect="equal",
        norm=LogNorm(vmin=vmin, vmax=vmax, clip=True),
        cmap=colormap,
    )

    # Place the panel-specific title above the axes.
    ax.set_title(title)

    # Apply the shared coordinate directions, ticks, labels, and equal aspect ratio.
    configure_density_axes(ax)

    # Return the image so its normalization and colormap can drive a colorbar.
    return image


def save_two_by_two_figure(
    ion_native: np.ndarray,
    electron_native: np.ndarray,
    ion_smoothed: np.ndarray,
    electron_smoothed: np.ndarray,
    ion_x_range_m: np.ndarray,
    ion_z_range_m: np.ndarray,
    electron_x_range_m: np.ndarray,
    electron_z_range_m: np.ndarray,
    simulation_time_s: float,
    simulation_step: int,
    smoothing_window: int,
    smoothing_passes: int,
    vmin: float,
    vmax: float,
    output_path: Path,
) -> None:
    """Save the native-versus-smoothed 2x2 ion/electron comparison figure."""

    # Create four square data panels arranged in two rows and two columns.
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(14.5, 11.5),
        constrained_layout=True,
    )

    # Draw the native ion density in the upper-left panel with the plasma colormap.
    ion_native_image = draw_density(
        axes[0, 0],
        ion_native,
        ion_x_range_m,
        ion_z_range_m,
        (
            f"Ion density — native {ion_native.shape[0]}×{ion_native.shape[1]}\n"
            f"step {simulation_step}, t = {simulation_time_s:.3e} s"
        ),
        "plasma",
        vmin,
        vmax,
    )

    # Draw the native electron density in the upper-right panel with viridis.
    electron_native_image = draw_density(
        axes[0, 1],
        electron_native,
        electron_x_range_m,
        electron_z_range_m,
        (
            f"Electron density — native {electron_native.shape[0]}×"
            f"{electron_native.shape[1]}\n"
            f"step {simulation_step}, t = {simulation_time_s:.3e} s"
        ),
        "viridis",
        vmin,
        vmax,
    )

    # Draw the full-resolution smoothed ion density in the lower-left panel.
    draw_density(
        axes[1, 0],
        ion_smoothed,
        ion_x_range_m,
        ion_z_range_m,
        (
            f"Ion density — smoothed full resolution\n"
            f"{smoothing_window}×{smoothing_window} boxcar × {smoothing_passes}; "
            "zeros included"
        ),
        "plasma",
        vmin,
        vmax,
    )

    # Draw the full-resolution smoothed electron density in the lower-right panel.
    draw_density(
        axes[1, 1],
        electron_smoothed,
        electron_x_range_m,
        electron_z_range_m,
        (
            f"Electron density — smoothed full resolution\n"
            f"{smoothing_window}×{smoothing_window} boxcar × {smoothing_passes}; "
            "zeros included"
        ),
        "viridis",
        vmin,
        vmax,
    )

    # Add one plasma colorbar shared by the two ion panels in the left column.
    ion_colorbar = figure.colorbar(
        ion_native_image,
        ax=[axes[0, 0], axes[1, 0]],
        pad=0.02,
        shrink=0.92,
    )

    # Label the ion colorbar with physical units.
    ion_colorbar.set_label(r"Ion density $(\mathrm{m}^{-3})$")

    # Add one viridis colorbar shared by the two electron panels in the right column.
    electron_colorbar = figure.colorbar(
        electron_native_image,
        ax=[axes[0, 1], axes[1, 1]],
        pad=0.02,
        shrink=0.92,
    )

    # Label the electron colorbar with physical units.
    electron_colorbar.set_label(r"Electron density $(\mathrm{m}^{-3})$")

    # Add a concise figure-level title describing the masking threshold.
    figure.suptitle(
        "Ion and electron density: native and smoothed "
        f"(values below {vmin:.0e} m$^{{-3}}$ shown white)",
        fontsize=16,
    )

    # Create the output directory if it does not already exist.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save the complete 2x2 figure at publication-friendly resolution.
    figure.savefig(output_path, dpi=220)

    # Close the figure to release memory before building the second figure.
    plt.close(figure)


def crop_z_half(
    density_xz: np.ndarray,
    z_range_m: np.ndarray,
    keep_negative: bool,
) -> tuple[np.ndarray, np.ndarray]:
    """Crop the density to Z < 0 or Z >= 0 without rescaling the retained data."""

    # Construct one physical Z coordinate for every density column.
    z_coordinates_m = coordinate_centers(z_range_m, density_xz.shape[1])

    # Select only negative-Z cells for the ion half when requested.
    if keep_negative:
        # Keep the original ion data only where Z is strictly less than zero.
        keep_mask = z_coordinates_m < 0.0

    # Select only nonnegative-Z cells for the electron half otherwise.
    else:
        # Keep the original electron data only where Z is zero or positive.
        keep_mask = z_coordinates_m >= 0.0

    # Ensure that the requested half actually contains at least one data column.
    if not np.any(keep_mask):
        # Raise an error if the coordinate range does not cross zero as expected.
        raise ValueError("The stored Z range does not contain the requested half-domain.")

    # Crop the density columns without interpolating, stretching, or copying the other half.
    cropped_density = density_xz[:, keep_mask]

    # Extract the physical Z coordinates corresponding to the retained columns.
    cropped_z_coordinates_m = z_coordinates_m[keep_mask]

    # Define the retained physical extent from its first and last actual coordinates.
    cropped_z_range_m = np.array(
        [cropped_z_coordinates_m[0], cropped_z_coordinates_m[-1]],
        dtype=np.float64,
    )

    # Return both the genuinely cropped data and its genuinely cropped coordinate extent.
    return cropped_density, cropped_z_range_m


def save_square_split_figure(
    ion_smoothed: np.ndarray,
    electron_smoothed: np.ndarray,
    ion_x_range_m: np.ndarray,
    ion_z_range_m: np.ndarray,
    electron_x_range_m: np.ndarray,
    electron_z_range_m: np.ndarray,
    vmin: float,
    vmax: float,
    output_path: Path,
) -> None:
    """Save a square figure with ion Z<0 on the left and electron Z>=0 on the right."""

    # Crop the ion array to its actual negative-Z columns only.
    ion_left, ion_left_z_range_m = crop_z_half(
        ion_smoothed,
        ion_z_range_m,
        keep_negative=True,
    )

    # Crop the electron array to its actual nonnegative-Z columns only.
    electron_right, electron_right_z_range_m = crop_z_half(
        electron_smoothed,
        electron_z_range_m,
        keep_negative=False,
    )

    # Mask ion values below vmin so low nonzero values and zeros appear white.
    plotted_ion = mask_below_display_minimum(ion_left, vmin)

    # Mask electron values below vmin using the same display rule.
    plotted_electron = mask_below_display_minimum(electron_right, vmin)

    # Prepare the plasma colormap for the negative-Z ion half.
    ion_colormap = make_density_colormap("plasma")

    # Prepare the viridis colormap for the nonnegative-Z electron half.
    electron_colormap = make_density_colormap("viridis")

    # Convert the ion X limits from meters to micrometers.
    ion_x_min_um = float(ion_x_range_m[0]) * METERS_TO_MICROMETERS

    # Convert the upper ion X limit from meters to micrometers.
    ion_x_max_um = float(ion_x_range_m[1]) * METERS_TO_MICROMETERS

    # Convert the negative-Z ion limits to micrometers.
    ion_z_min_um = float(ion_left_z_range_m[0]) * METERS_TO_MICROMETERS

    # Convert the ion half-domain's upper Z limit to micrometers.
    ion_z_max_um = float(ion_left_z_range_m[1]) * METERS_TO_MICROMETERS

    # Convert the electron X limits from meters to micrometers.
    electron_x_min_um = float(electron_x_range_m[0]) * METERS_TO_MICROMETERS

    # Convert the upper electron X limit from meters to micrometers.
    electron_x_max_um = float(electron_x_range_m[1]) * METERS_TO_MICROMETERS

    # Convert the nonnegative-Z electron limits to micrometers.
    electron_z_min_um = float(electron_right_z_range_m[0]) * METERS_TO_MICROMETERS

    # Convert the electron half-domain's upper Z limit to micrometers.
    electron_z_max_um = float(electron_right_z_range_m[1]) * METERS_TO_MICROMETERS

    # Create a square canvas so the complete output image has a 1:1 width-to-height ratio.
    figure = plt.figure(figsize=(12.0, 12.0), constrained_layout=True)

    # Create one central axes for both cropped halves.
    axes = figure.add_subplot(111)

    # Explicitly force the axes box itself to remain square.
    axes.set_box_aspect(1.0)

    # Draw only the actual negative-Z ion columns in the left half of the shared axes.
    ion_image = axes.imshow(
        plotted_ion,
        origin="lower",
        extent=[ion_z_min_um, ion_z_max_um, ion_x_min_um, ion_x_max_um],
        interpolation="nearest",
        aspect="equal",
        norm=LogNorm(vmin=vmin, vmax=vmax, clip=True),
        cmap=ion_colormap,
    )

    # Draw only the actual nonnegative-Z electron columns in the right half.
    electron_image = axes.imshow(
        plotted_electron,
        origin="lower",
        extent=[
            electron_z_min_um,
            electron_z_max_um,
            electron_x_min_um,
            electron_x_max_um,
        ],
        interpolation="nearest",
        aspect="equal",
        norm=LogNorm(vmin=vmin, vmax=vmax, clip=True),
        cmap=electron_colormap,
    )

    # Apply the shared coordinate limits, flipped X direction, labels, and equal scaling.
    configure_density_axes(axes)

    # Reassert a square axes box after applying the common axis settings.
    axes.set_box_aspect(1.0)

    # Draw a thin divider exactly at the boundary between the two retained datasets.
    axes.axvline(0.0, color="black", linewidth=0.9)

    # Describe which species is shown on each side of the composite.
    axes.set_title("Smoothed density: ion for Z < 0, electron for Z ≥ 0", fontsize=15)

    # Add the plasma colorbar on the left side of the square composite.
    ion_colorbar = figure.colorbar(
        ion_image,
        ax=axes,
        location="left",
        pad=0.03,
        shrink=0.80,
    )

    # Label the left colorbar as ion density.
    ion_colorbar.set_label(r"Ion density $(\mathrm{m}^{-3})$")

    # Add the viridis colorbar on the right side of the square composite.
    electron_colorbar = figure.colorbar(
        electron_image,
        ax=axes,
        location="right",
        pad=0.03,
        shrink=0.80,
    )

    # Label the right colorbar as electron density.
    electron_colorbar.set_label(r"Electron density $(\mathrm{m}^{-3})$")

    # Create the output directory if needed.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save the square composite figure.
    figure.savefig(output_path, dpi=220)

    # Close the figure to release its memory and file handles.
    plt.close(figure)


def averaged_x_slice_about_z_zero(
    ion_density_xz: np.ndarray,
    ion_x_range_m: np.ndarray,
    ion_z_range_m: np.ndarray,
    number_of_z_cells: int = 6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Average the native ion density over the Z cells closest to zero."""

    # Require at least one native Z cell in the averaging window.
    if number_of_z_cells < 1:
        # Stop immediately when the requested averaging width is invalid.
        raise ValueError("number_of_z_cells must be at least 1.")

    # Prevent a request for more Z cells than the native ion array contains.
    if number_of_z_cells > ion_density_xz.shape[1]:
        # Report the available number of native Z columns in the error message.
        raise ValueError(
            f"Requested {number_of_z_cells} Z cells, but only "
            f"{ion_density_xz.shape[1]} are available."
        )

    # Construct one native physical X coordinate for every ion-density row.
    x_coordinates_m = coordinate_centers(
        ion_x_range_m,
        ion_density_xz.shape[0],
    )

    # Construct one native physical Z coordinate for every ion-density column.
    z_coordinates_m = coordinate_centers(
        ion_z_range_m,
        ion_density_xz.shape[1],
    )

    # Sort all native Z-column indices by their absolute distance from Z = 0.
    indices_by_distance = np.argsort(np.abs(z_coordinates_m))

    # Retain exactly the six native Z columns nearest the central plane by default.
    selected_indices = indices_by_distance[:number_of_z_cells]

    # Sort the retained indices into increasing physical Z order for reporting.
    selected_indices = np.sort(selected_indices)

    # Extract only the selected native ion-density columns around Z = 0.
    selected_density_columns = ion_density_xz[:, selected_indices]

    # Average the six selected native columns along the Z direction.
    averaged_density_profile = selected_density_columns.mean(axis=1)

    # Convert the native X coordinates from meters to micrometers for plotting.
    x_coordinates_um = x_coordinates_m * METERS_TO_MICROMETERS

    # Convert the selected native Z-cell centers to micrometers for annotation.
    selected_z_coordinates_um = (
        z_coordinates_m[selected_indices] * METERS_TO_MICROMETERS
    )

    # Return the linear X coordinates, averaged density, and exact Z cells used.
    return x_coordinates_um, averaged_density_profile, selected_z_coordinates_um


def save_ion_x_slice_figure(
    ion_density_xz: np.ndarray,
    ion_x_range_m: np.ndarray,
    ion_z_range_m: np.ndarray,
    output_path: Path,
    number_of_z_cells: int = 6,
) -> None:
    """Save the requested logarithmic ion-density lineout versus linear X."""

    # Build the one-dimensional native ion-density profile around Z = 0.
    (
        x_coordinates_um,
        averaged_density_profile,
        selected_z_coordinates_um,
    ) = averaged_x_slice_about_z_zero(
        ion_density_xz,
        ion_x_range_m,
        ion_z_range_m,
        number_of_z_cells=number_of_z_cells,
    )

    # Replace nonpositive values with NaN because logarithmic axes require positivity.
    plotted_density_profile = np.where(
        averaged_density_profile > 0.0,
        averaged_density_profile,
        np.nan,
    )

    # Create one rectangular figure for the one-dimensional lineout.
    figure, axes = plt.subplots(
        figsize=(8.0, 6.0),
        constrained_layout=True,
    )

    # Draw the ion-density line in orange, as requested.
    axes.plot(
        x_coordinates_um,
        plotted_density_profile,
        color="orange",
        linewidth=2.0,
    )

    # Label the horizontal coordinate as X in micrometers.
    axes.set_xlabel(r"$X\;(\mu\mathrm{m})$")

    # Label the logarithmic vertical coordinate as ion density.
    axes.set_ylabel(r"Ion density $(\mathrm{m}^{-3})$")

    # Display ion density on a logarithmic vertical axis.
    axes.set_yscale("log")

    # Restrict the linear X axis to the requested interval from -7 to +13 micrometers.
    axes.set_xlim(-7.0, 13.0)

    # Restrict the logarithmic density axis to the requested displayed interval.
    axes.set_ylim(4.0e24, 2.0e27)

    # Add a light major grid to help read values from the logarithmic profile.
    axes.grid(
        True,
        which="major",
        alpha=0.25,
    )

    # Mark X = 0 with a thin neutral vertical reference line.
    axes.axvline(
        0.0,
        color="black",
        linewidth=0.8,
        alpha=0.5,
    )

    # Describe the native six-cell average used to form the lineout.
    axes.set_title(
        "Ion-density slice through X averaged over 6 native Z cells about Z = 0"
    )

    # Format each selected native Z-cell center for an explicit figure annotation.
    z_label = ", ".join(
        f"{value:.3f}"
        for value in selected_z_coordinates_um
    )

    # Record the exact six Z-cell centers used in the average inside the plot.
    axes.text(
        0.02,
        0.02,
        f"Averaged Z-cell centers (μm): {z_label}",
        transform=axes.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        bbox=dict(
            boxstyle="round",
            facecolor="white",
            alpha=0.8,
            edgecolor="0.8",
        ),
    )

    # Create the selected output directory when it does not already exist.
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Save the completed one-dimensional ion-density figure at high resolution.
    figure.savefig(
        output_path,
        dpi=220,
    )

    # Close the figure to release memory and prevent duplicate displays later.
    plt.close(figure)


def parse_arguments() -> argparse.Namespace:
    """Parse command-line paths and plotting controls."""

    # Create the command-line parser using the module documentation as help text.
    parser = argparse.ArgumentParser(description=__doc__)

    # Require the ion BP directory or split-component prefix.
    parser.add_argument(
        "--ion",
        required=True,
        help="Ion BP directory or split-component prefix, such as /path/to/i.",
    )

    # Require the electron BP directory or split-component prefix.
    parser.add_argument(
        "--electron",
        required=True,
        help="Electron BP directory or split-component prefix, such as /path/to/e.",
    )

    # Allow both output figures to be written to a chosen directory.
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for the three output PNG files.",
    )

    # Allow the lower logarithmic display threshold to be changed.
    parser.add_argument(
        "--vmin",
        type=float,
        default=DEFAULT_VMIN,
        help="Minimum displayed density; lower values are white.",
    )

    # Allow the upper logarithmic display threshold to be changed.
    parser.add_argument(
        "--vmax",
        type=float,
        default=DEFAULT_VMAX,
        help="Maximum displayed density.",
    )

    # Allow the odd boxcar width to be changed from its default value of five.
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=DEFAULT_SMOOTH_WINDOW,
        help="Odd boxcar width used at full native resolution.",
    )

    # Allow the number of repeated boxcar passes to be adjusted.
    parser.add_argument(
        "--smooth-passes",
        type=int,
        default=DEFAULT_SMOOTH_PASSES,
        help="Number of full-resolution boxcar smoothing passes.",
    )

    # Parse the actual arguments supplied by the user.
    arguments = parser.parse_args()

    # Validate the logarithmic normalization before any expensive data processing occurs.
    if arguments.vmin <= 0.0 or arguments.vmax <= arguments.vmin:
        # Reject invalid limits because LogNorm requires 0 < vmin < vmax.
        parser.error("Require 0 < --vmin < --vmax for logarithmic plotting.")

    # Return the validated command-line namespace.
    return arguments


def main() -> None:
    """Read both datasets, smooth them, and write the three requested figures."""

    # Parse all input paths and plotting options from the command line.
    arguments = parse_arguments()

    # Read the native ion density and its coordinate and time metadata.
    ion_density, ion_x_range_m, ion_z_range_m, ion_time_s, ion_step = read_density(
        arguments.ion
    )

    # Read the native electron density and its coordinate and time metadata.
    (
        electron_density,
        electron_x_range_m,
        electron_z_range_m,
        electron_time_s,
        electron_step,
    ) = read_density(arguments.electron)

    # Verify that the ion and electron outputs correspond to the same simulation step.
    if electron_step != ion_step:
        # Stop rather than comparing fields from different steps accidentally.
        raise ValueError(
            f"Ion step {ion_step} does not match electron step {electron_step}."
        )

    # Verify that the ion and electron outputs correspond to the same simulation time.
    if not np.isclose(electron_time_s, ion_time_s, rtol=1.0e-12, atol=0.0):
        # Stop rather than comparing fields from different times accidentally.
        raise ValueError(
            f"Ion time {ion_time_s:.16e} s does not match "
            f"electron time {electron_time_s:.16e} s."
        )

    # Smooth the ion field at its native resolution while including zeros in every average.
    ion_smoothed = boxcar_smooth_include_zeros(
        ion_density,
        window=arguments.smooth_window,
        passes=arguments.smooth_passes,
    )

    # Smooth the electron field with the identical native-resolution operation.
    electron_smoothed = boxcar_smooth_include_zeros(
        electron_density,
        window=arguments.smooth_window,
        passes=arguments.smooth_passes,
    )

    # Define the exact filename of the first requested 2x2 figure.
    two_by_two_path = (
        arguments.output_dir / "density_2x2_zeros_included_whitebelowvmin.png"
    )

    # Generate and save the first figure containing all four density panels.
    save_two_by_two_figure(
        ion_native=ion_density,
        electron_native=electron_density,
        ion_smoothed=ion_smoothed,
        electron_smoothed=electron_smoothed,
        ion_x_range_m=ion_x_range_m,
        ion_z_range_m=ion_z_range_m,
        electron_x_range_m=electron_x_range_m,
        electron_z_range_m=electron_z_range_m,
        simulation_time_s=ion_time_s,
        simulation_step=ion_step,
        smoothing_window=arguments.smooth_window,
        smoothing_passes=arguments.smooth_passes,
        vmin=arguments.vmin,
        vmax=arguments.vmax,
        output_path=two_by_two_path,
    )

    # Define the exact filename of the second requested square split figure.
    square_split_path = (
        arguments.output_dir / "density_smoothed_combined_half_figure_square.png"
    )

    # Generate and save the square figure using only the correct Z half of each species.
    save_square_split_figure(
        ion_smoothed=ion_smoothed,
        electron_smoothed=electron_smoothed,
        ion_x_range_m=ion_x_range_m,
        ion_z_range_m=ion_z_range_m,
        electron_x_range_m=electron_x_range_m,
        electron_z_range_m=electron_z_range_m,
        vmin=arguments.vmin,
        vmax=arguments.vmax,
        output_path=square_split_path,
    )

    # Report the native ion array size so the user can confirm no interpolation was used.
    print(f"Ion native shape: {ion_density.shape}")

    # Report the native electron array size for the same reason.
    print(f"Electron native shape: {electron_density.shape}")

    # Report the first output path clearly.
    print(f"Wrote 2x2 figure: {two_by_two_path}")

    # Define the exact filename of the third requested ion-density lineout.
    ion_slice_path = (
        arguments.output_dir / "ion_density_x_slice_avg6cells_about_z0.png"
    )

    # Generate the third figure from the native ion density without interpolation.
    save_ion_x_slice_figure(
        ion_density_xz=ion_density,
        ion_x_range_m=ion_x_range_m,
        ion_z_range_m=ion_z_range_m,
        output_path=ion_slice_path,
        number_of_z_cells=6,
    )

    # Report the second output path clearly.
    print(f"Wrote square split figure: {square_split_path}")

    # Report the third output path clearly.
    print(f"Wrote ion X-slice figure: {ion_slice_path}")


# Run main only when this file is executed directly from the command line.
if __name__ == "__main__":
    # Execute the complete read-smooth-plot workflow.
    main()
