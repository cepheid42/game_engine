import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

# ============================================================
# Simple analytical axisymmetric mirror field
# ============================================================
# Field model:
#   1) uniform guide field B_GUIDE e_z
#   2) two identical finite end solenoids at z = +/- Z_SOL
#
# The off-axis field is computed from the on-axis field using the
# lowest-order axisymmetric, divergence-free expansion:
#
#   B_r(r,z) = - r/2 dB0/dz
#   B_z(r,z) = B0(z) - r^2/4 d^2B0/dz^2
#
# This is intended as a clean orbit-code benchmark, not an engineering
# reconstruction of GDMT.
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

# Deliberately low guide field so the 14 keV deuteron Larmor radius is visible.
# Uniform guide field; chosen low enough to make the Larmor radius visible.
B_GUIDE = 0.08          # T
# Approximate target magnetic field near each end mirror.
B_MIRROR_TARGET = 0.40  # T, mirror ratio about 5
# Axial locations of the two end solenoids.
Z_SOL = 2.0             # m, end-solenoid centers
# Effective axial length of each finite solenoid.
SOL_LENGTH = 0.40       # m
# Effective radius of each finite solenoid.
SOL_RADIUS = 0.45       # m

# Initial particle position.
X0, Y0, Z0 = 0.10, 0.0, 0.0
# Axial loss boundary.
Z_LOSS = 2.55
# Simple cylindrical radial wall boundary.
R_WALL = 1.25

# Boris-pusher timestep.
DT = 1.0e-8
# Maximum integration time.
TMAX = 1.2e-5
# Store diagnostics every N timesteps.
SAMPLE_EVERY = 2

# Representative launch pitch angles [deg].
PITCH_ANGLES = [15.0, 25.0, 30.0, 45.0, 60.0, 80.0]

# Output directory for figures.
OUTDIR = Path(__file__).resolve().parent / "simple_analytic_mirror_outputs_commented"
# Create the output directory if needed.
OUTDIR.mkdir(parents=True, exist_ok=True)


# Convert a nonrelativistic kinetic energy in keV into speed.
def speed_from_energy_kev(Ekev):
    # Use E = 1/2 m v^2 after converting keV to joules.
    return np.sqrt(2.0 * Ekev * 1.0e3 * QE / M_D)


# Dimensionless on-axis field shape for a finite solenoid.
def finite_solenoid_shape(z, zc):
    """Dimensionless on-axis finite-solenoid shape per mu0*nI/2."""
    # Left edge of the solenoid.
    z1 = zc - 0.5 * SOL_LENGTH
    # Right edge of the solenoid.
    z2 = zc + 0.5 * SOL_LENGTH
    # Distance from evaluation point to left edge.
    s1 = z - z1
    # Distance from evaluation point to right edge.
    s2 = z - z2
    return s1 / np.sqrt(SOL_RADIUS**2 + s1**2) - s2 / np.sqrt(SOL_RADIUS**2 + s2**2)


# First derivative of the finite-solenoid shape with respect to z.
def finite_solenoid_shape_d1(z, zc):
    # Left edge of the solenoid.
    z1 = zc - 0.5 * SOL_LENGTH
    # Right edge of the solenoid.
    z2 = zc + 0.5 * SOL_LENGTH
    # Distance from evaluation point to left edge.
    s1 = z - z1
    # Distance from evaluation point to right edge.
    s2 = z - z2
    return (
        SOL_RADIUS**2 / (SOL_RADIUS**2 + s1**2)**1.5
        - SOL_RADIUS**2 / (SOL_RADIUS**2 + s2**2)**1.5
    )


# Second derivative of the finite-solenoid shape with respect to z.
def finite_solenoid_shape_d2(z, zc):
    # Left edge of the solenoid.
    z1 = zc - 0.5 * SOL_LENGTH
    # Right edge of the solenoid.
    z2 = zc + 0.5 * SOL_LENGTH
    # Distance from evaluation point to left edge.
    s1 = z - z1
    # Distance from evaluation point to right edge.
    s2 = z - z2
    return (
        -3.0 * SOL_RADIUS**2 * s1 / (SOL_RADIUS**2 + s1**2)**2.5
        +3.0 * SOL_RADIUS**2 * s2 / (SOL_RADIUS**2 + s2**2)**2.5
    )


# Choose the solenoid ampere-turns so the mirror field is near the target value.
def choose_ampere_turns():
    """Choose NI so B(0,Z_SOL) approximately equals B_MIRROR_TARGET."""
    # Add self and opposite-solenoid contributions at the right mirror center.
    shape_at_mirror = finite_solenoid_shape(Z_SOL, +Z_SOL) + finite_solenoid_shape(Z_SOL, -Z_SOL)
    # B = B_GUIDE + mu0*(NI/L)/2*shape
    return (B_MIRROR_TARGET - B_GUIDE) * 2.0 * SOL_LENGTH / (MU0 * shape_at_mirror)


# Ampere-turns used for each identical end solenoid.
NI_SOL = choose_ampere_turns()


# On-axis field B0(z), including the guide field and both finite solenoids.
def B0_axis(z):
    # Support scalar or array z inputs.
    z = np.asarray(z, dtype=float)
    # Effective current density in ampere-turns per meter.
    nI = NI_SOL / SOL_LENGTH
    return B_GUIDE + 0.5 * MU0 * nI * (
        finite_solenoid_shape(z, +Z_SOL) + finite_solenoid_shape(z, -Z_SOL)
    )


# First derivative dB0/dz.
def dB0_dz(z):
    # Support scalar or array z inputs.
    z = np.asarray(z, dtype=float)
    # Effective current density in ampere-turns per meter.
    nI = NI_SOL / SOL_LENGTH
    return 0.5 * MU0 * nI * (
        finite_solenoid_shape_d1(z, +Z_SOL) + finite_solenoid_shape_d1(z, -Z_SOL)
    )


# Second derivative d^2B0/dz^2.
def d2B0_dz2(z):
    # Support scalar or array z inputs.
    z = np.asarray(z, dtype=float)
    # Effective current density in ampere-turns per meter.
    nI = NI_SOL / SOL_LENGTH
    return 0.5 * MU0 * nI * (
        finite_solenoid_shape_d2(z, +Z_SOL) + finite_solenoid_shape_d2(z, -Z_SOL)
    )


# Build the axisymmetric magnetic field in cylindrical components.
def B_cyl(r, z):
    r = np.asarray(r, dtype=float)
    # Support scalar or array z inputs.
    z = np.asarray(z, dtype=float)
    # Find the common shape for vectorized r,z evaluation.
    shape = np.broadcast_shapes(r.shape, z.shape)
    # Broadcast radius to that common shape.
    rr = np.broadcast_to(r, shape)
    # Broadcast z to that common shape.
    zz = np.broadcast_to(z, shape)

    # Evaluate the on-axis field at all z values.
    B0 = B0_axis(zz)
    # Evaluate dB0/dz.
    B1 = dB0_dz(zz)
    # Evaluate d^2B0/dz^2.
    B2 = d2B0_dz2(zz)

    # Lowest-order radial field required by div(B)=0.
    Br = -0.5 * rr * B1
    # Lowest-order off-axis correction to Bz.
    Bz = B0 - 0.25 * rr**2 * B2
    return Br, Bz


# Convert the axisymmetric field to Cartesian Bx, By, Bz for particle pushing.
def B_xyz(X):
    # Particle x coordinates.
    x = X[:, 0]
    # Particle y coordinates.
    y = X[:, 1]
    # Particle z coordinates.
    z = X[:, 2]
    # Cylindrical radius.
    r = np.hypot(x, y)
    # Evaluate cylindrical field components.
    Br, Bz = B_cyl(r, z)
    # Allocate Cartesian field array.
    B = np.zeros_like(X)
    # Avoid division by zero on the axis.
    mask = r > 0.0
    # Convert radial field to Bx.
    B[mask, 0] = Br[mask] * x[mask] / r[mask]
    # Convert radial field to By.
    B[mask, 1] = Br[mask] * y[mask] / r[mask]
    # Copy the axial field into Bz.
    B[:, 2] = Bz
    return B


# Advance particles by one timestep using the magnetic-only Boris pusher.
def boris_push(X, V, dt):
    # Gather Bx, By, Bz at particle positions.
    B = B_xyz(X)
    # Boris half-rotation vector.
    t = QOM_D * B * (0.5 * dt)
    # Magnitude squared of t.
    t2 = np.sum(t * t, axis=1)
    # Boris s vector.
    s = 2.0 * t / (1.0 + t2)[:, None]
    # First cross product of the rotation.
    vp = V + np.cross(V, t)
    # Complete the velocity rotation.
    Vnew = V + np.cross(vp, s)
    # Move particles with the updated velocity.
    Xnew = X + Vnew * dt
    return Xnew, Vnew


# Compute magnetic moment mu = m v_perp^2/(2B).
def magnetic_moment(V, B):
    # Magnetic-field magnitude.
    Bmag = np.linalg.norm(B, axis=1)
    # Unit vector along B with zero-division protection.
    b = B / np.maximum(Bmag, 1.0e-300)[:, None]
    # Parallel velocity.
    vpar = np.sum(V * b, axis=1)
    # Total speed squared.
    v2 = np.sum(V * V, axis=1)
    # Perpendicular speed squared.
    vperp2 = np.maximum(v2 - vpar**2, 0.0)
    return M_D * vperp2 / (2.0 * np.maximum(Bmag, 1.0e-300))


# Return True for particles outside the axial domain or radial wall.
def lost_mask(X):
    r = np.hypot(X[:, 0], X[:, 1])
    return (np.abs(X[:, 2]) >= Z_LOSS) | (r >= R_WALL)


# Compute on-axis field, effective mirror ratio, and loss-cone angle.
def axis_quantities():
    z = np.linspace(0.0, Z_LOSS, 2001)
    Baxis = np.abs(B0_axis(z))
    B0 = np.linalg.norm(B_xyz(np.array([[X0, Y0, Z0]], dtype=float))[0])
    Rm = np.max(Baxis) / B0
    alpha_lc = np.degrees(np.arcsin(1.0 / np.sqrt(Rm)))
    return z, Baxis, B0, Rm, alpha_lc


# Predict mirror point using B_mirror = B0/sin^2(alpha).
def predicted_turning_point(alpha_deg, z_axis, B_axis, B0):
    target = B0 / np.sin(np.radians(alpha_deg))**2
    if target > np.max(B_axis):
        return np.nan
    idx = np.where(B_axis >= target)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    if i == 0:
        return z_axis[0]
    z0, z1 = z_axis[i - 1], z_axis[i]
    b0, b1 = B_axis[i - 1], B_axis[i]
    return z0 + (target - b0) * (z1 - z0) / (b1 - b0)


# Integrate representative particle orbits at several pitch angles.
def run_representative():
    vmag = speed_from_energy_kev(E_KEV)
    z_axis, B_axis, B0, Rm, alpha_lc = axis_quantities()
    results = []

    for pitch in PITCH_ANGLES:
        alpha = np.radians(pitch)
        X = np.array([[X0, Y0, Z0]], dtype=float)
        V = np.array([[vmag * np.sin(alpha), 0.0, vmag * np.cos(alpha)]], dtype=float)

        mu0 = magnetic_moment(V, B_xyz(X))[0]

        xs, ys, zs, ts, mus = [], [], [], [], []
        zturn = 0.0
        lost = False

        for i in range(int(np.ceil(TMAX / DT)) + 1):
            zturn = max(zturn, abs(X[0, 2]))
            if i % SAMPLE_EVERY == 0:
                xs.append(X[0, 0])
                ys.append(X[0, 1])
                zs.append(X[0, 2])
                ts.append(i * DT)
                mus.append(magnetic_moment(V, B_xyz(X))[0])

            if lost_mask(X)[0]:
                lost = True
                break

            X, V = boris_push(X, V, DT)

        results.append({
            "pitch": pitch,
            "x": np.array(xs),
            "y": np.array(ys),
            "z": np.array(zs),
            "t": np.array(ts),
            "mu": np.array(mus),
            "mu0": mu0,
            "lost": lost,
            "status": "lost" if lost else "trapped",
            "zturn_num": np.nan if lost else zturn,
            "zturn_pred": predicted_turning_point(pitch, z_axis, B_axis, B0),
        })

    return results, z_axis, B_axis, B0, Rm, alpha_lc


# Scan launch pitch angle to estimate the numerical loss cone.
def pitch_scan():
    vmag = speed_from_energy_kev(E_KEV)
    scan = np.arange(10.0, 45.01, 0.5)
    alpha = np.radians(scan)

    X = np.zeros((len(scan), 3))
    X[:, 0] = X0
    X[:, 1] = Y0
    X[:, 2] = Z0

    V = np.zeros_like(X)
    V[:, 0] = vmag * np.sin(alpha)
    V[:, 2] = vmag * np.cos(alpha)

    lost = np.zeros(len(scan), dtype=bool)
    active = np.ones(len(scan), dtype=bool)

    for _ in range(int(np.ceil(TMAX / DT))):
        if not np.any(active):
            break
        Xa, Va = X[active], V[active]
        Xa, Va = boris_push(Xa, Va, DT)
        X[active], V[active] = Xa, Va
        newly_lost = lost_mask(X)
        lost |= newly_lost
        active &= ~newly_lost

    trapped = ~lost
    lost_angles = scan[~trapped]
    trapped_angles = scan[trapped]
    lower = lost_angles.max() if lost_angles.size else np.nan
    upper = trapped_angles.min() if trapped_angles.size else np.nan
    mid = 0.5 * (lower + upper) if np.isfinite(lower) and np.isfinite(upper) else np.nan
    return scan, trapped, lower, upper, mid


# Run the benchmark and generate all figures.
def plot_all():
    results, z_axis, B_axis, B0, Rm, alpha_lc = run_representative()

    x = np.linspace(-1.3, 1.3, 241)
    z = np.linspace(-2.7, 2.7, 381)
    Xg, Zg = np.meshgrid(x, z, indexing="xy")
    Rg = np.abs(Xg)
    Brg, Bzg = B_cyl(Rg, Zg)
    Bxg = np.zeros_like(Xg)
    mask = Rg > 0.0
    Bxg[mask] = Brg[mask] * Xg[mask] / Rg[mask]
    Bmag = np.sqrt(Bxg**2 + Bzg**2)
    Bplot = np.clip(Bmag, 0.03, 0.6)

    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    levels = np.geomspace(0.03, 0.6, 50)
    cf = ax.contourf(Xg, Zg, Bplot, levels=levels, norm=LogNorm(vmin=0.03, vmax=0.6))
    cb = fig.colorbar(cf, ax=ax)
    # Label several magnetic-field values on the log colorbar.
    cb.set_ticks([0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60])
    # Use explicit decimal labels so the colorbar is easier to read.
    cb.set_ticklabels(["0.03", "0.04", "0.05", "0.06", "0.08", "0.10", "0.15", "0.20", "0.30", "0.40", "0.60"])
    # Label the plotted magnetic-field magnitude.
    cb.set_label("|B| [T]")
    ax.streamplot(x, z, Bxg, Bzg, density=1.0, linewidth=0.55)
    for res in results:
        ax.plot(res["x"], res["z"], linewidth=1.2, label=f"{res['pitch']:.0f}° {res['status']}")
    ax.axhline(+Z_SOL, linestyle="--", linewidth=1)
    ax.axhline(-Z_SOL, linestyle="--", linewidth=1)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    ax.set_title("Simple analytical mirror: XZ trajectories with visible Larmor radius")
    ax.set_aspect("equal")
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-2.7, 2.7)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUTDIR / "simple_mirror_trajectories_xz.png", dpi=190)
    plt.close(fig)

    fig = plt.figure(figsize=(8.2, 7.0))
    ax = fig.add_subplot(111, projection="3d")
    for res in results:
        ax.plot(res["x"], res["y"], res["z"], label=f"{res['pitch']:.0f}° {res['status']}")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_title("Simple analytical mirror: 3D orbits")
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_zlim(-2.7, 2.7)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "simple_mirror_trajectories_3d.png", dpi=190)
    plt.close(fig)

    zfull = np.linspace(-2.7, 2.7, 1001)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(zfull, np.abs(B0_axis(zfull)))
    ax.axhline(B_GUIDE, linestyle="--", label="guide field")
    ax.axhline(B_MIRROR_TARGET, linestyle=":", label="mirror target")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("|B(r=0,z)| [T]")
    ax.set_title("Simple analytical mirror: on-axis field")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "simple_mirror_axis_field.png", dpi=190)
    plt.close(fig)

    trapped = [r for r in results if not r["lost"] and np.isfinite(r["zturn_pred"])]
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    if trapped:
        ax.plot([r["pitch"] for r in trapped], [r["zturn_pred"] for r in trapped], marker="o", label="adiabatic prediction")
        ax.plot([r["pitch"] for r in trapped], [r["zturn_num"] for r in trapped], marker="s", label="Boris orbit")
    ax.set_xlabel("pitch angle [deg]")
    ax.set_ylabel(r"turning point $|z_{turn}|$ [m]")
    ax.set_title("Simple analytical mirror: turning point comparison")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "simple_mirror_turning_points.png", dpi=190)
    plt.close(fig)

    scan, trapped_scan, lower, upper, mid = pitch_scan()
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(scan, trapped_scan.astype(int), marker="o", markersize=3, linewidth=1)
    ax.axvline(alpha_lc, linestyle="--", label=f"adiabatic alpha_LC={alpha_lc:.1f} deg")
    if np.isfinite(mid):
        ax.axvline(mid, linestyle=":", label=f"numerical alpha_LC~{mid:.1f} deg")
    ax.set_xlabel("pitch angle [deg]")
    ax.set_ylabel("classification (1=trapped, 0=lost)")
    ax.set_ylim(-0.1, 1.1)
    ax.set_title("Simple analytical mirror: loss-cone scan")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "simple_mirror_pitch_scan.png", dpi=190)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    for res in results:
        mu_rel = (res["mu"] - res["mu0"]) / res["mu0"]
        ax.plot(res["t"] * 1e6, mu_rel, label=f"{res['pitch']:.0f}° {res['status']}")
    ax.axhline(0, linewidth=1)
    ax.set_xlabel("time [microseconds]")
    ax.set_ylabel(r"$(\mu-\mu_0)/\mu_0$")
    ax.set_title("Simple analytical mirror: magnetic moment conservation")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "simple_mirror_mu_conservation.png", dpi=190)
    plt.close(fig)

    return results, z_axis, B_axis, B0, Rm, alpha_lc


if __name__ == "__main__":
    results, z_axis, B_axis, B0, Rm, alpha_lc = plot_all()
    print("Simple analytical mirror parameters")
    print(f"  E = {E_KEV:.1f} keV deuterons")
    print(f"  B_guide = {B_GUIDE:.3f} T")
    print(f"  B_mirror_target = {B_MIRROR_TARGET:.3f} T")
    print(f"  NI per end solenoid = {NI_SOL:.3e} A-turns")
    print(f"  effective mirror ratio = {Rm:.3f}")
    print(f"  adiabatic loss cone = {alpha_lc:.2f} deg")
    v = speed_from_energy_kev(E_KEV)
    rho60 = v * np.sin(np.radians(60)) / (QOM_D * B0)
    print(f"  approximate launch rho_L at 60 deg = {rho60:.3f} m")
    print()
    print("Representative orbit summary")
    for r in results:
        mu_rel = (r["mu"] - r["mu0"]) / r["mu0"]
        print(
            f"  {r['pitch']:5.1f} deg {r['status']:7s} "
            f"z_pred={r['zturn_pred'] if np.isfinite(r['zturn_pred']) else np.nan:6.3f} m "
            f"z_num={r['zturn_num'] if np.isfinite(r['zturn_num']) else np.nan:6.3f} m "
            f"max|dmu/mu|={np.max(np.abs(mu_rel)):.3e}"
        )
    print(f"\nSaved plots in: {OUTDIR}")
