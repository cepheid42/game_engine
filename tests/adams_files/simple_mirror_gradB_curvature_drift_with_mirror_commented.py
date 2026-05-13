import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

# ============================================================
# Grad-B and curvature drift test in a simple analytical mirror
# ============================================================
#
# Field:
#   uniform guide field + two finite end solenoids.
#
# The pusher uses Bx, By, Bz.  The field is generated from an
# analytical on-axis field B0(z), with the lowest-order
# axisymmetric divergence-free expansion:
#
#   Bx = Br x/r,  By = Br y/r
#   Br = -r/2 dB0/dz
#   Bz = B0 - r^2/4 d2B0/dz2
#
# Drift theory:
#   v_gradB = (m v_perp^2 / (2 q B^3)) B x grad(B)
#   v_curv  = (m v_parallel^2 / (q B)) b x kappa
#   kappa   = (b dot grad) b
#
# In an axisymmetric mirror, these drifts are mainly azimuthal.
# We compare the measured gyrocenter azimuthal drift speed with
# the predicted grad-B + curvature drift speed.
# ============================================================

# Elementary charge [C].
QE = 1.602176634e-19
# Deuteron mass [kg].
M_D = 3.343583719e-27
# Vacuum permeability [H/m].
MU0 = 4.0e-7 * np.pi
# Deuteron charge-to-mass ratio [C/kg].
QOM_D = QE / M_D

# Test-particle kinetic energy [keV].
E_KEV = 14.0

# Uniform guide field [T].
B_GUIDE = 0.30          # T
# Approximate target on-axis field near each mirror throat [T].
B_MIRROR_TARGET = 1.50  # T
# Axial center location of each end solenoid [m].
Z_SOL = 2.0             # m
# Effective axial solenoid length [m].
SOL_LENGTH = 0.40       # m
# Effective solenoid radius [m].
SOL_RADIUS = 0.45       # m

# Off-axis launch is essential: at r=0 the azimuthal drift is not visible
# as a guiding-center circle.
# Off-axis particle launch point; x0 is nonzero so azimuthal drift is visible.
X0, Y0, Z0 = 0.22, 0.0, 0.0

# Launch pitch angle for the drift benchmark [deg].
PITCH_ANGLE_DEG = 60.0
# Axial loss boundary [m].
Z_LOSS = 2.55
# Simple cylindrical wall radius [m].
R_WALL = 1.25

# Boris-pusher timestep [s].
DT = 1.0e-9
# Maximum simulated time [s].
TMAX = 2.5e-5
# Store trajectory and diagnostic data every N timesteps.
SAMPLE_EVERY = 2

# Output directory for plots.
OUTDIR = Path(__file__).resolve().parent / "simple_mirror_drift_outputs_quieter_commented"
# Create the output directory if it does not already exist.
OUTDIR.mkdir(parents=True, exist_ok=True)


# Convert nonrelativistic kinetic energy in keV into particle speed.
def speed_from_energy_kev(Ekev):
    # Use E = 1/2 m v^2, with keV converted to joules.
    return np.sqrt(2.0 * Ekev * 1.0e3 * QE / M_D)


# Dimensionless on-axis field shape for one finite solenoid.
def finite_solenoid_shape(z, zc):
    # Left axial edge of the finite solenoid.
    z1 = zc - 0.5 * SOL_LENGTH
    # Right axial edge of the finite solenoid.
    z2 = zc + 0.5 * SOL_LENGTH
    # Distance from the evaluation point to the left edge.
    s1 = z - z1
    # Distance from the evaluation point to the right edge.
    s2 = z - z2
    return s1 / np.sqrt(SOL_RADIUS**2 + s1**2) - s2 / np.sqrt(SOL_RADIUS**2 + s2**2)


# First axial derivative of the finite-solenoid shape.
def finite_solenoid_shape_d1(z, zc):
    # Left axial edge of the finite solenoid.
    z1 = zc - 0.5 * SOL_LENGTH
    # Right axial edge of the finite solenoid.
    z2 = zc + 0.5 * SOL_LENGTH
    # Distance from the evaluation point to the left edge.
    s1 = z - z1
    # Distance from the evaluation point to the right edge.
    s2 = z - z2
    return (
        SOL_RADIUS**2 / (SOL_RADIUS**2 + s1**2)**1.5
        - SOL_RADIUS**2 / (SOL_RADIUS**2 + s2**2)**1.5
    )


# Second axial derivative of the finite-solenoid shape.
def finite_solenoid_shape_d2(z, zc):
    # Left axial edge of the finite solenoid.
    z1 = zc - 0.5 * SOL_LENGTH
    # Right axial edge of the finite solenoid.
    z2 = zc + 0.5 * SOL_LENGTH
    # Distance from the evaluation point to the left edge.
    s1 = z - z1
    # Distance from the evaluation point to the right edge.
    s2 = z - z2
    return (
        -3.0 * SOL_RADIUS**2 * s1 / (SOL_RADIUS**2 + s1**2)**2.5
        +3.0 * SOL_RADIUS**2 * s2 / (SOL_RADIUS**2 + s2**2)**2.5
    )


# Choose solenoid ampere-turns that produce the requested mirror field.
def choose_ampere_turns():
    # Include both the nearby solenoid and the opposite solenoid at the mirror center.
    shape_at_mirror = finite_solenoid_shape(Z_SOL, +Z_SOL) + finite_solenoid_shape(Z_SOL, -Z_SOL)
    return (B_MIRROR_TARGET - B_GUIDE) * 2.0 * SOL_LENGTH / (MU0 * shape_at_mirror)


# Total ampere-turns for each identical end solenoid.
NI_SOL = choose_ampere_turns()


# Analytical on-axis magnetic field B0(z).
def B0_axis(z):
    # Convert scalar or array input to a NumPy array.
    z = np.asarray(z, dtype=float)
    # Convert total ampere-turns into effective ampere-turns per meter.
    nI = NI_SOL / SOL_LENGTH
    return B_GUIDE + 0.5 * MU0 * nI * (
        finite_solenoid_shape(z, +Z_SOL) + finite_solenoid_shape(z, -Z_SOL)
    )


# Analytical first derivative dB0/dz.
def dB0_dz(z):
    # Convert scalar or array input to a NumPy array.
    z = np.asarray(z, dtype=float)
    # Convert total ampere-turns into effective ampere-turns per meter.
    nI = NI_SOL / SOL_LENGTH
    return 0.5 * MU0 * nI * (
        finite_solenoid_shape_d1(z, +Z_SOL) + finite_solenoid_shape_d1(z, -Z_SOL)
    )


# Analytical second derivative d^2B0/dz^2.
def d2B0_dz2(z):
    # Convert scalar or array input to a NumPy array.
    z = np.asarray(z, dtype=float)
    # Convert total ampere-turns into effective ampere-turns per meter.
    nI = NI_SOL / SOL_LENGTH
    return 0.5 * MU0 * nI * (
        finite_solenoid_shape_d2(z, +Z_SOL) + finite_solenoid_shape_d2(z, -Z_SOL)
    )


# Compute the axisymmetric magnetic field in cylindrical components.
def B_cyl(r, z):
    # Convert radius input to a NumPy array.
    r = np.asarray(r, dtype=float)
    # Convert scalar or array input to a NumPy array.
    z = np.asarray(z, dtype=float)
    # Determine the common vectorized shape for r and z.
    shape = np.broadcast_shapes(r.shape, z.shape)
    # Broadcast radius onto the common evaluation grid.
    rr = np.broadcast_to(r, shape)
    # Broadcast axial position onto the common evaluation grid.
    zz = np.broadcast_to(z, shape)

    # Evaluate the on-axis field.
    B0 = B0_axis(zz)
    # Evaluate the first axial derivative.
    B1 = dB0_dz(zz)
    # Evaluate the second axial derivative.
    B2 = d2B0_dz2(zz)

    # Lowest-order divergence-free radial component.
    Br = -0.5 * rr * B1
    # Lowest-order off-axis correction to the axial component.
    Bz = B0 - 0.25 * rr**2 * B2
    return Br, Bz


# Convert the axisymmetric field to Cartesian Bx, By, Bz at particle positions.
def B_xyz_at_points(X):
    # Extract x positions.
    x = X[:, 0]
    # Extract y positions.
    y = X[:, 1]
    # Extract z positions.
    z = X[:, 2]
    # Cylindrical radius from Cartesian position.
    r = np.hypot(x, y)
    # Evaluate the cylindrical field.
    Br, Bz = B_cyl(r, z)
    # Allocate Cartesian magnetic-field array.
    B = np.zeros_like(X)
    # Avoid division by zero on the symmetry axis.
    mask = r > 0.0
    # Convert radial component into Bx.
    B[mask, 0] = Br[mask] * x[mask] / r[mask]
    # Convert radial component into By.
    B[mask, 1] = Br[mask] * y[mask] / r[mask]
    # Assign the axial component Bz.
    B[:, 2] = Bz
    return B


# Convenience wrapper for evaluating Bx, By, Bz at one point.
def B_xyz_single(x, y, z):
    return B_xyz_at_points(np.array([[x, y, z]], dtype=float))[0]


# Advance particle positions and velocities with the magnetic-only Boris pusher.
def boris_push(X, V, dt):
    # Gather Bx, By, Bz at the particle positions.
    B = B_xyz_at_points(X)
    # Boris half-rotation vector.
    t = QOM_D * B * (0.5 * dt)
    # Squared magnitude of the rotation vector.
    t2 = np.sum(t * t, axis=1)
    # Boris s vector for the second rotation cross product.
    s = 2.0 * t / (1.0 + t2)[:, None]
    # Intermediate rotated velocity.
    vp = V + np.cross(V, t)
    # Final rotated velocity.
    Vnew = V + np.cross(vp, s)
    # Update position using the new velocity.
    Xnew = X + Vnew * dt
    return Xnew, Vnew


# Flag particles that leave the axial domain or hit the radial wall.
def lost_mask(X):
    r = np.hypot(X[:, 0], X[:, 1])
    return (np.abs(X[:, 2]) >= Z_LOSS) | (r >= R_WALL)


# Compute magnetic moment for one particle.
def magnetic_moment_single(v, B):
    # Magnetic-field magnitude.
    # Magnetic-field magnitude.
    Bmag = np.linalg.norm(B)
    b = B / max(Bmag, 1e-300)
    # Parallel velocity at this point.
    vpar = np.dot(v, b)
    # Perpendicular speed squared.
    vperp2 = max(np.dot(v, v) - vpar**2, 0.0)
    return M_D * vperp2 / (2.0 * max(Bmag, 1e-300))


# Estimate the first-order guiding-center position from particle phase-space coordinates.
def guiding_center_position(x, v):
    # Evaluate magnetic field at the particle position.
    B = B_xyz_single(x[0], x[1], x[2])
    # Magnetic-field magnitude.
    Bmag = np.linalg.norm(B)
    # Magnetic-field unit vector.
    b = B / Bmag
    # Cyclotron frequency for a positive deuteron.
    Omega = QOM_D * Bmag
    # Positive-charge guiding-center estimate:
    # R_gc = r + v x b / Omega.
    return x + np.cross(v, b) / Omega


# Compute numerical grad(|B|) by centered differences.
def grad_B(x, h=1.0e-4):
    """Numerical grad |B| in Cartesian coordinates."""
    # Allocate Cartesian gradient vector.
    grad = np.zeros(3)
    # Differentiate one Cartesian direction at a time.
    for j in range(3):
        # Small finite-difference displacement.
        dx = np.zeros(3)
        dx[j] = h
        Bp = np.linalg.norm(B_xyz_single(*(x + dx)))
        Bm = np.linalg.norm(B_xyz_single(*(x - dx)))
        grad[j] = (Bp - Bm) / (2.0 * h)
    return grad


# Return the local magnetic-field unit vector.
def unit_b(x):
    # Evaluate magnetic field at the particle position.
    B = B_xyz_single(x[0], x[1], x[2])
    return B / np.linalg.norm(B)


# Compute magnetic-field curvature kappa = (b dot grad)b.
def curvature_kappa(x, h=1.0e-4):
    """Compute kappa = (b dot grad) b by finite difference along b."""
    b0 = unit_b(x)
    bp = unit_b(x + h * b0)
    bm = unit_b(x - h * b0)
    return (bp - bm) / (2.0 * h)


# Evaluate grad-B drift, curvature drift, and their sum at the gyrocenter.
def drift_theory(x_gc, v):
    """Return grad-B, curvature, and total drift velocities at the guiding center."""
    # Magnetic field at the gyrocenter.
    B = B_xyz_single(x_gc[0], x_gc[1], x_gc[2])
    # Magnetic-field magnitude.
    Bmag = np.linalg.norm(B)
    # Magnetic-field unit vector.
    b = B / Bmag

    # Parallel velocity at this point.
    vpar = np.dot(v, b)
    # Perpendicular speed squared.
    vperp2 = max(np.dot(v, v) - vpar**2, 0.0)

    # Magnetic-field strength gradient.
    gB = grad_B(x_gc)
    # Magnetic-field curvature vector.
    kappa = curvature_kappa(x_gc)

    v_gradB = M_D * vperp2 / (2.0 * QE * Bmag**3) * np.cross(B, gB)
    v_curv = M_D * vpar**2 / (QE * Bmag) * np.cross(b, kappa)
    return v_gradB, v_curv, v_gradB + v_curv


# Project a Cartesian vector onto the local azimuthal direction.
def cylindrical_vphi_from_vector(x, v):
    r = np.hypot(x[0], x[1])
    if r <= 0.0:
        return 0.0
    ephi = np.array([-x[1] / r, x[0] / r, 0.0])
    return np.dot(v, ephi)


# Run the full-particle orbit and store drift-theory diagnostics.
def run_orbit():
    # Particle speed from the kinetic energy.
    vmag = speed_from_energy_kev(E_KEV)
    # Pitch angle in radians.
    pitch = np.radians(PITCH_ANGLE_DEG)

    # One-particle position array.
    X = np.array([[X0, Y0, Z0]], dtype=float)
    # Initial velocity with perpendicular x and parallel z components.
    V = np.array([[vmag * np.sin(pitch), 0.0, vmag * np.cos(pitch)]], dtype=float)

    # Raw particle trajectory storage.
    xs, ys, zs, ts = [], [], [], []
    # Guiding-center trajectory storage.
    xgcs, ygcs, zgcs = [], [], []
    # Gyrocenter azimuthal angle history.
    phi_gc = []
    vphi_theory_total = []
    vphi_theory_gradB = []
    vphi_theory_curv = []
    vphi_measured_inst = []
    mu_rel = []

    mu0 = magnetic_moment_single(V[0], B_xyz_single(X[0,0], X[0,1], X[0,2]))

    # First pass: store guiding center and theory at sample points.
    for i in range(int(np.ceil(TMAX / DT)) + 1):
        if i % SAMPLE_EVERY == 0:
            x = X[0].copy()
            v = V[0].copy()
            xgc = guiding_center_position(x, v)

            vgB, vcurv, vtot = drift_theory(xgc, v)

            xs.append(x[0]); ys.append(x[1]); zs.append(x[2]); ts.append(i * DT)
            xgcs.append(xgc[0]); ygcs.append(xgc[1]); zgcs.append(xgc[2])
            phi_gc.append(np.arctan2(xgc[1], xgc[0]))
            vphi_theory_gradB.append(cylindrical_vphi_from_vector(xgc, vgB))
            vphi_theory_curv.append(cylindrical_vphi_from_vector(xgc, vcurv))
            vphi_theory_total.append(cylindrical_vphi_from_vector(xgc, vtot))
            mu_rel.append((magnetic_moment_single(v, B_xyz_single(*x)) - mu0) / mu0)

        if lost_mask(X)[0]:
            break

        X, V = boris_push(X, V, DT)

    ts = np.array(ts)
    xgcs = np.array(xgcs)
    ygcs = np.array(ygcs)
    zgcs = np.array(zgcs)
    phi_gc = np.unwrap(np.array(phi_gc))
    r_gc = np.hypot(xgcs, ygcs)

    # Measured azimuthal drift from finite-difference of gyrocenter angle.
    dphi_dt = np.gradient(phi_gc, ts)
    vphi_measured = r_gc * dphi_dt

    return {
        "t": ts,
        "x": np.array(xs),
        "y": np.array(ys),
        "z": np.array(zs),
        "xgc": xgcs,
        "ygc": ygcs,
        "zgc": zgcs,
        "r_gc": r_gc,
        "phi_gc": phi_gc,
        "vphi_measured": vphi_measured,
        "vphi_gradB": np.array(vphi_theory_gradB),
        "vphi_curv": np.array(vphi_theory_curv),
        "vphi_total": np.array(vphi_theory_total),
        "mu_rel": np.array(mu_rel),
    }


# Smooth noisy time series with a simple moving average.
def moving_average(y, n=21):
    if len(y) < n:
        return y.copy()
    kernel = np.ones(n) / n
    return np.convolve(y, kernel, mode="same")



# Predict the mirror point from the adiabatic mirror condition.
def predicted_mirror_point_for_pitch(alpha_deg, z_axis, B_axis, B_launch):
    """Adiabatic mirror point from B(z_m) = B_launch / sin^2(alpha)."""
    target = B_launch / np.sin(np.radians(alpha_deg))**2
    if target > np.max(B_axis):
        return np.nan, target
    idx = np.where(B_axis >= target)[0]
    if len(idx) == 0:
        return np.nan, target
    i = idx[0]
    if i == 0:
        return z_axis[0], target
    z0, z1 = z_axis[i - 1], z_axis[i]
    b0, b1 = B_axis[i - 1], B_axis[i]
    z_m = z0 + (target - b0) * (z1 - z0) / (b1 - b0)
    return float(z_m), float(target)


# Compare the numerical maximum |z| with the adiabatic mirror point.
def compute_mirror_point_diagnostic(data):
    """Compare numerical bounce point to adiabatic mirror theory for this orbit."""
    z_axis = np.linspace(0.0, Z_LOSS, 2001)
    B_axis = np.abs(B0_axis(z_axis))
    B_launch = np.linalg.norm(B_xyz_single(X0, Y0, Z0))
    z_theory, B_target = predicted_mirror_point_for_pitch(PITCH_ANGLE_DEG, z_axis, B_axis, B_launch)

    # Use the gyrocenter z trajectory for the numerical mirror point.
    # This avoids some gyro-scale contamination in the raw particle z.
    z_num_gc = float(np.max(np.abs(data["zgc"])))
    z_num_particle = float(np.max(np.abs(data["z"])))

    return {
        "z_axis": z_axis,
        "B_axis": B_axis,
        "B_launch": B_launch,
        "B_target": B_target,
        "z_theory": z_theory,
        "z_num_gc": z_num_gc,
        "z_num_particle": z_num_particle,
        "dz_gc": z_num_gc - z_theory,
        "dz_particle": z_num_particle - z_theory,
    }


# Generate field, orbit, drift, mirror-point, and conservation plots.
def plot_results(data):
    # Field topology and orbit.
    # X grid for the XZ field plot.
    x = np.linspace(-1.2, 1.2, 241)
    # Z grid for the XZ field plot.
    z = np.linspace(-2.7, 2.7, 381)
    # 2-D plotting grid.
    Xg, Zg = np.meshgrid(x, z, indexing="xy")
    # Cylindrical radius in the XZ plane.
    Rg = np.abs(Xg)
    # Cylindrical magnetic field on the plotting grid.
    Brg, Bzg = B_cyl(Rg, Zg)
    # Cartesian Bx on the plotting grid.
    Bxg = np.zeros_like(Xg)
    # Avoid division by zero on axis.
    mask = Rg > 0
    # Convert radial field to Bx in the XZ plane.
    Bxg[mask] = Brg[mask] * Xg[mask] / Rg[mask]
    # Magnetic-field magnitude for the color map.
    Bmag = np.sqrt(Bxg**2 + Bzg**2)
    # Clip the displayed magnitude range to improve color contrast.
    Bplot = np.clip(Bmag, 0.05, 0.8)

    fig, ax = plt.subplots(figsize=(7.2, 7.0))
    # Log-spaced magnetic-field contour levels.
    levels = np.geomspace(0.05, 0.8, 55)
    # Filled contour plot of |B|.
    cf = ax.contourf(Xg, Zg, Bplot, levels=levels, norm=LogNorm(vmin=0.05, vmax=0.8))
    # Magnetic-field magnitude colorbar.
    cb = fig.colorbar(cf, ax=ax)
    # Add more labeled field values to the colorbar.
    cb.set_ticks([0.05, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80])
    # Use explicit labels so the log colorbar is easy to read.
    cb.set_ticklabels(["0.05", "0.06", "0.08", "0.10", "0.15", "0.20", "0.30", "0.40", "0.60", "0.80"])
    # Label the colorbar quantity and units.
    cb.set_label("|B| [T]")
    # Overlay magnetic-field streamlines.
    ax.streamplot(x, z, Bxg, Bzg, density=1.0, linewidth=0.55)
    # Overlay the raw particle orbit.
    ax.plot(data["x"], data["z"], linewidth=0.8, label="particle orbit")
    # Overlay the estimated gyrocenter orbit.
    ax.plot(data["xgc"], data["zgc"], linewidth=2.0, label="gyrocenter estimate")
    ax.axhline(+Z_SOL, linestyle="--", linewidth=1)
    ax.axhline(-Z_SOL, linestyle="--", linewidth=1)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    ax.set_title("Grad-B/curvature drift test: orbit and gyrocenter in XZ")
    ax.set_aspect("equal")
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-2.7, 2.7)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "drift_orbit_xz.png", dpi=190)
    plt.close(fig)

    # XY projection: shows azimuthal precession.
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    ax.plot(data["x"], data["y"], linewidth=0.8, label="particle orbit")
    ax.plot(data["xgc"], data["ygc"], linewidth=2.0, label="gyrocenter")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Azimuthal drift in XY projection")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "drift_orbit_xy.png", dpi=190)
    plt.close(fig)

    # v_phi comparison.
    t_us = data["t"] * 1.0e6
    vmeas = moving_average(data["vphi_measured"], 31)
    vtot = moving_average(data["vphi_total"], 31)
    vgB = moving_average(data["vphi_gradB"], 31)
    vcurv = moving_average(data["vphi_curv"], 31)

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.plot(t_us, vmeas, label="measured gyrocenter v_phi")
    ax.plot(t_us, vtot, "--", label="theory: grad-B + curvature")
    ax.plot(t_us, vgB, ":", label="theory: grad-B only")
    ax.plot(t_us, vcurv, "-.", label="theory: curvature only")
    ax.set_xlabel("time [microseconds]")
    ax.set_ylabel("azimuthal drift speed v_phi [m/s]")
    ax.set_title("Measured gyrocenter drift vs drift theory")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "drift_vphi_comparison.png", dpi=190)
    plt.close(fig)

    # Integrated azimuthal angle comparison.
    phi_meas = data["phi_gc"] - data["phi_gc"][0]
    # Integrate theory dphi/dt = vphi/r.
    dphi_dt_theory = data["vphi_total"] / np.maximum(data["r_gc"], 1e-12)
    dt = np.gradient(data["t"])
    phi_theory = np.cumsum(dphi_dt_theory * dt)
    phi_theory -= phi_theory[0]

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.plot(t_us, phi_meas, label="measured gyrocenter angle")
    ax.plot(t_us, phi_theory, "--", label="integrated drift theory")
    ax.set_xlabel("time [microseconds]")
    ax.set_ylabel("azimuthal angle change Δphi [rad]")
    ax.set_title("Azimuthal precession: measured vs drift theory")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "drift_phi_comparison.png", dpi=190)
    plt.close(fig)

    # Magnetic moment.
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.plot(t_us, data["mu_rel"])
    ax.axhline(0.0, linewidth=1.0)
    ax.set_xlabel("time [microseconds]")
    ax.set_ylabel("(mu - mu0)/mu0")
    ax.set_title("Magnetic moment conservation in drift test")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTDIR / "drift_mu_conservation.png", dpi=190)
    plt.close(fig)

    # On-axis field.
    z = np.linspace(-2.7, 2.7, 1001)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(z, np.abs(B0_axis(z)))
    ax.axhline(B_GUIDE, linestyle="--", label="guide field")
    ax.axhline(B_MIRROR_TARGET, linestyle=":", label="mirror target")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("|B(r=0,z)| [T]")
    ax.set_title("Simple mirror field for drift benchmark")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "drift_axis_field.png", dpi=190)
    plt.close(fig)

    # Mirror-point comparison to adiabatic theory.
    mirror = compute_mirror_point_diagnostic(data)

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(mirror["z_axis"], mirror["B_axis"], label=r"$|B(r=0,z)|$")
    ax.axhline(mirror["B_target"], linestyle="--", label=r"$B_{\rm launch}/\sin^2\alpha$")
    ax.axvline(mirror["z_theory"], linestyle="--", label="theory mirror point")
    ax.axvline(mirror["z_num_gc"], linestyle=":", label="gyrocenter numerical mirror point")
    ax.axvline(mirror["z_num_particle"], linestyle="-.", label="particle numerical max |z|")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("|B| [T]")
    ax.set_title("Mirror-point comparison to adiabatic theory")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "drift_mirror_point_comparison.png", dpi=190)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(data["t"] * 1.0e6, np.abs(data["zgc"]), label=r"gyrocenter $|z|$")
    ax.plot(data["t"] * 1.0e6, np.abs(data["z"]), alpha=0.5, label=r"particle $|z|$")
    ax.axhline(mirror["z_theory"], linestyle="--", label="theory mirror point")
    ax.set_xlabel("time [microseconds]")
    ax.set_ylabel(r"$|z|$ [m]")
    ax.set_title("Numerical bounce position vs adiabatic mirror point")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "drift_mirror_point_time_history.png", dpi=190)
    plt.close(fig)

    # Error summary.
    sl = slice(20, -20) if len(vmeas) > 80 else slice(None)
    rms_err = np.sqrt(np.mean((vmeas[sl] - vtot[sl])**2))
    rms_theory = np.sqrt(np.mean(vtot[sl]**2))
    mean_meas = np.mean(vmeas[sl])
    mean_theory = np.mean(vtot[sl])
    phi_err = phi_meas[-1] - phi_theory[-1]

    return rms_err, rms_theory, mean_meas, mean_theory, phi_meas[-1], phi_theory[-1], phi_err, mirror


# Run the benchmark when executed as a script.
if __name__ == "__main__":
    data = run_orbit()
    rms_err, rms_theory, mean_meas, mean_theory, phi_meas, phi_theory, phi_err, mirror = plot_results(data)

    B_launch = np.linalg.norm(B_xyz_single(X0, Y0, Z0))
    v = speed_from_energy_kev(E_KEV)
    rho = v * np.sin(np.radians(PITCH_ANGLE_DEG)) / (QOM_D * B_launch)

    print("Grad-B and curvature drift benchmark")
    print(f"  E = {E_KEV:.1f} keV deuterons")
    print(f"  pitch angle = {PITCH_ANGLE_DEG:.1f} deg")
    print(f"  launch point = ({X0:.3f}, {Y0:.3f}, {Z0:.3f}) m")
    print(f"  B_launch = {B_launch:.4f} T")
    print(f"  launch Larmor radius = {rho:.4f} m")
    print(f"  NI per end solenoid = {NI_SOL:.3e} A-turns")
    print()
    print("Drift comparison")
    print(f"  mean measured gyrocenter v_phi = {mean_meas:.3e} m/s")
    print(f"  mean theory v_phi              = {mean_theory:.3e} m/s")
    print(f"  RMS(v_meas - v_theory)         = {rms_err:.3e} m/s")
    print(f"  RMS(v_theory)                  = {rms_theory:.3e} m/s")
    print(f"  relative RMS error             = {rms_err / max(rms_theory, 1e-300):.3f}")
    print()
    print("Integrated azimuthal precession")
    print(f"  measured Delta phi = {phi_meas:.4f} rad")
    print(f"  theory Delta phi   = {phi_theory:.4f} rad")
    print(f"  difference         = {phi_err:+.4f} rad")
    print()
    print("Mirror-point comparison")
    print(f"  B_launch                  = {mirror['B_launch']:.5f} T")
    print(f"  B_launch/sin^2(alpha)     = {mirror['B_target']:.5f} T")
    print(f"  theory mirror point       = {mirror['z_theory']:.5f} m")
    print(f"  numerical gyrocenter max |z| = {mirror['z_num_gc']:.5f} m")
    print(f"  numerical particle max |z|   = {mirror['z_num_particle']:.5f} m")
    print(f"  gyrocenter error          = {mirror['dz_gc']:+.3e} m")
    print(f"  particle error            = {mirror['dz_particle']:+.3e} m")
    print()
    print(f"Saved plots in: {OUTDIR}")
