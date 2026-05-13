
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# E x B and polarization-drift benchmarks with a
# leapfrog-consistent Boris push
# ============================================================
# This script contains three focused tests:
#
#   1) Uniform static E x B drift:
#        E = E0 xhat, B = B0 zhat.
#
#   2) Strong-gradient nonuniform static E x B-like motion:
#        E_x(x) = E0 + G x, B = B0 zhat.
#
#   3) Magnetized polarization drift:
#        E_x(t) = Eamp sin(omega t), B = B0 zhat.
#
# The second case is analytically solvable because the electric field is
# linear in x.  The particle executes a modified cyclotron orbit while
# drifting in y, and the local drift at the oscillation center is much larger
# than in Case 1 so the gradient effect is visually obvious.
#
# The third case is compared with the exact finite-frequency driven solution.
# The first-order polarization drift and the accompanying time-dependent
# E(t) x B drift are also plotted as low-frequency approximations.
#
# ============================================================

# -----------------------------
# Constants
# -----------------------------

# Elementary charge [C].
QE = 1.602176634e-19

# Deuteron mass [kg].
M_D = 3.343583719e-27

# Deuteron charge-to-mass ratio [C/kg].
QOM_D = QE / M_D

# Deuteron kinetic energy used to set the perpendicular speed scale [keV].
E_KEV = 14.0

# Uniform guide magnetic field [T].
B0 = 0.50

# Uniform-case electric-field strength [V/m].
E0_UNIFORM = 2.0e4

# Base electric-field strength for the nonuniform linear-E test [V/m].
E0_LINEAR = 2.0e4

# Cyclotron frequency [rad/s].
OMEGA = QOM_D * B0

# Cyclotron period [s].
T_GYRO = 2.0 * np.pi / OMEGA

# Timestep small enough that both orbit phase and drift agree cleanly with theory.
DT = 2.0e-10

# Output directory.
OUTDIR = Path(__file__).resolve().parent / "static_exb_three_case_outputs"
OUTDIR.mkdir(parents=True, exist_ok=True)


# -----------------------------
# Shared helpers
# -----------------------------

def speed_from_energy_kev(Ekev):
    # Convert kinetic energy from keV to speed using E = 1/2 m v^2.
    return np.sqrt(2.0 * Ekev * 1.0e3 * QE / M_D)


def boris_velocity_step(v_half_minus, x_integer, dt, electric_field_func):
    # Advance v_{n-1/2} to v_{n+1/2} using E(x_n) and uniform B.
    # This is the leapfrog-compatible Boris force step.

    # Electric field at the integer-time position.
    E = electric_field_func(x_integer)

    # Uniform magnetic field vector.
    B = np.array([0.0, 0.0, B0])

    # First half electric acceleration.
    v_minus = v_half_minus + QOM_D * E * (0.5 * dt)

    # Boris t vector.
    tvec = QOM_D * B * (0.5 * dt)

    # Boris s vector.
    svec = 2.0 * tvec / (1.0 + np.dot(tvec, tvec))

    # First magnetic rotation.
    v_prime = v_minus + np.cross(v_minus, tvec)

    # Second magnetic rotation.
    v_plus = v_minus + np.cross(v_prime, svec)

    # Second half electric acceleration.
    v_half_plus = v_plus + QOM_D * E * (0.5 * dt)

    # Return half-step velocity.
    return v_half_plus


def reconstruct_centered_velocity(v_half_minus, x_integer, electric_field_func):
    # Reconstruct an approximate centered velocity v_n from v_{n-1/2}.
    # We do this by advancing only a half timestep using the same force law.
    return boris_velocity_step(v_half_minus, x_integer, 0.5 * DT, electric_field_func)


def run_leapfrog_orbit(X0, V0_centered, tmax, electric_field_func, exact_velocity_func_for_init):
    # Number of integer timesteps.
    nsteps = int(round(tmax / DT))

    # Integer-time sample array.
    t = np.arange(nsteps + 1) * DT

    # Allocate position history.
    Xhist = np.zeros((nsteps + 1, 3))

    # Allocate centered velocity history used only for diagnostics.
    Vcenter_hist = np.zeros((nsteps + 1, 3))

    # Initialize v_{-1/2} from the exact solution where possible.
    # This removes the start-up cyclotron transient from the benchmark.
    v_half = exact_velocity_func_for_init(-0.5 * DT).copy()

    # Initial position at t=0.
    x = X0.copy()

    # Main leapfrog loop.
    for n in range(nsteps + 1):

        # Store x_n.
        Xhist[n] = x

        # Store reconstructed centered velocity v_n.
        Vcenter_hist[n] = reconstruct_centered_velocity(v_half, x, electric_field_func)

        # Stop after final sample.
        if n == nsteps:
            break

        # Advance velocity from v_{n-1/2} to v_{n+1/2}.
        v_half_new = boris_velocity_step(v_half, x, DT, electric_field_func)

        # Advance position using v_{n+1/2}.
        x = x + v_half_new * DT

        # Shift half-step velocity.
        v_half = v_half_new

    # Return integer-time diagnostics.
    return t, Xhist, Vcenter_hist


def guiding_center_uniform_B(X, Vcenter, vE_local):
    # Unit vector along B.
    b = np.array([0.0, 0.0, 1.0])

    # Velocity relative to the chosen drift velocity.
    u = Vcenter - vE_local

    # Uniform-field positive-charge guiding-center estimate.
    return X + np.cross(u, b) / OMEGA


def rms_error(numerical, exact):
    # RMS absolute error.
    return np.sqrt(np.mean((numerical - exact)**2))


def rms_relative_error(numerical, exact):
    # RMS absolute error.
    err = rms_error(numerical, exact)

    # RMS scale of exact quantity.
    scale = np.sqrt(np.mean(exact**2))

    # Relative error.
    rel = err / max(scale, 1.0e-300)

    # Return absolute, scale, and relative.
    return err, scale, rel


# ============================================================
# Case 1: uniform static E x B drift
# ============================================================

def uniform_electric_field(x):
    # Constant electric field E = E0 xhat.
    return np.array([E0_UNIFORM, 0.0, 0.0])


def uniform_exact_solution(t, X0, V0_centered):
    # Exact E x B drift velocity.
    vE = np.array([0.0, -E0_UNIFORM / B0, 0.0])

    # Initial cyclotron velocity in the drift frame.
    u0 = V0_centered - vE

    # Initial cyclotron components.
    u0x = u0[0]
    u0y = u0[1]

    # Cyclotron phase.
    th = OMEGA * np.asarray(t)

    # Allocate exact arrays.
    Xex = np.zeros((len(np.atleast_1d(t)), 3))
    Vex = np.zeros_like(Xex)

    # Exact cyclotron velocity.
    ux = u0x * np.cos(th) + u0y * np.sin(th)
    uy = -u0x * np.sin(th) + u0y * np.cos(th)

    # Exact velocity.
    Vex[:, 0] = ux + vE[0]
    Vex[:, 1] = uy + vE[1]
    Vex[:, 2] = V0_centered[2]

    # Exact position.
    Xex[:, 0] = X0[0] + (u0x / OMEGA) * np.sin(th) + (u0y / OMEGA) * (1.0 - np.cos(th)) + vE[0] * t
    Xex[:, 1] = X0[1] + (u0x / OMEGA) * (np.cos(th) - 1.0) + (u0y / OMEGA) * np.sin(th) + vE[1] * t
    Xex[:, 2] = X0[2] + V0_centered[2] * t

    # Return exact arrays.
    return Xex, Vex, vE


def run_uniform_case():
    # 14 keV perpendicular speed.
    vgyro = speed_from_energy_kev(E_KEV)

    # Exact E x B drift.
    vE = np.array([0.0, -E0_UNIFORM / B0, 0.0])

    # Initial position.
    X0 = np.array([0.0, 0.0, 0.0])

    # Centered initial velocity: gyromotion plus exact E x B drift.
    V0_centered = np.array([vgyro, vE[1], 0.0])

    # Final time is an integer number of gyroperiods.
    tmax = 20 * T_GYRO

    # Exact velocity helper for half-step initialization.
    def exact_velocity_at_time(tt):
        _, Vtmp, _ = uniform_exact_solution(np.array([tt]), X0, V0_centered)
        return Vtmp[0]

    # Run leapfrog Boris.
    t, X, Vc = run_leapfrog_orbit(X0, V0_centered, tmax, uniform_electric_field, exact_velocity_at_time)

    # Exact solution.
    Xex, Vex, vE = uniform_exact_solution(t, X0, V0_centered)

    # Numerical and exact guiding centers.
    Rgc = guiding_center_uniform_B(X, Vc, vE)
    Rgc_ex = guiding_center_uniform_B(Xex, Vex, vE)

    # Plot trajectory.
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.plot(X[:, 0], X[:, 1], linewidth=1.0, label="Boris particle")
    ax.plot(Xex[:, 0], Xex[:, 1], "--", linewidth=1.3, label="exact particle")
    ax.plot(Rgc[:, 0], Rgc[:, 1], linewidth=2.0, label="Boris guiding center")
    ax.plot(Rgc_ex[:, 0], Rgc_ex[:, 1], "--", linewidth=1.5, label="exact guiding center")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Case 1: uniform static $E\\times B$ drift")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "case1_uniform_exb_trajectory.png", dpi=190)
    plt.close(fig)

    # Plot Xgc and Ygc.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(t * 1e6, Rgc[:, 0], label="Boris $X_{gc}$")
    ax.plot(t * 1e6, Rgc_ex[:, 0], "--", label="exact $X_{gc}$")
    ax.plot(t * 1e6, Rgc[:, 1], label="Boris $Y_{gc}$")
    ax.plot(t * 1e6, Rgc_ex[:, 1], "--", label="exact $Y_{gc}$")
    ax.set_xlabel("time [µs]")
    ax.set_ylabel("guiding-center coordinate [m]")
    ax.set_title("Case 1: guiding-center comparison")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(OUTDIR / "case1_uniform_exb_guiding_center.png", dpi=190)
    plt.close(fig)

    # Summarize errors.
    abs_xgc = rms_error(Rgc[:, 0], Rgc_ex[:, 0])
    abs_ygc = rms_error(Rgc[:, 1], Rgc_ex[:, 1])
    final_xgc = Rgc[-1, 0] - Rgc_ex[-1, 0]
    final_ygc = Rgc[-1, 1] - Rgc_ex[-1, 1]

    # Return diagnostics.
    return {
        "vE_y": vE[1],
        "rho": vgyro / OMEGA,
        "abs_xgc_rms": abs_xgc,
        "abs_ygc_rms": abs_ygc,
        "final_xgc_error": final_xgc,
        "final_ygc_error": final_ygc,
    }


# ============================================================
# Case 2: nonuniform linear E_x(x) with static B
# ============================================================

def make_linear_case_parameters():
    # Select a visible but stable field gradient.
    # The dimensionless parameter gamma/Omega^2 controls the modified gyrofrequency.
    eps_grad = 0.10

    # Convert gamma/Omega^2 to G = dEx/dx.
    G = eps_grad * OMEGA**2 / QOM_D

    # Modified oscillator frequency.
    kappa = np.sqrt(OMEGA**2 - QOM_D * G)

    # Choose the oscillation center/equilibrium x location.
    xeq = 0.30

    # Use 14 keV speed as the transverse speed scale.
    vgyro = speed_from_energy_kev(E_KEV)

    # Oscillation amplitude for x(t)=xeq + A cos(kappa t).
    A = vgyro / kappa

    # Constant C = vy + Omega x.
    C = ((OMEGA**2 - QOM_D * G) * xeq - QOM_D * E0_LINEAR) / OMEGA

    # Mean y drift.
    vbar = C - OMEGA * xeq

    # Return all parameters.
    return eps_grad, G, kappa, xeq, vgyro, A, C, vbar


def linear_electric_field_factory(G):
    # Build an electric-field function for E_x = E0 + G x.
    def linear_electric_field(x):
        # Return Ex, Ey, Ez at one particle position.
        return np.array([E0_LINEAR + G * x[0], 0.0, 0.0])
    return linear_electric_field


def linear_exact_solution(t, X0, params):
    # Unpack parameters.
    eps_grad, G, kappa, xeq, vgyro, A, C, vbar = params

    # Convert t to array.
    t = np.asarray(t)

    # Allocate exact arrays.
    Xex = np.zeros((len(np.atleast_1d(t)), 3))
    Vex = np.zeros_like(Xex)

    # Exact x and vx.
    Xex[:, 0] = xeq + A * np.cos(kappa * t)
    Vex[:, 0] = -A * kappa * np.sin(kappa * t)

    # Exact y and vy.
    Xex[:, 1] = X0[1] + vbar * t - (OMEGA * A / kappa) * np.sin(kappa * t)
    Vex[:, 1] = vbar - OMEGA * A * np.cos(kappa * t)

    # z remains zero.
    Xex[:, 2] = 0.0
    Vex[:, 2] = 0.0

    # Return exact arrays.
    return Xex, Vex


def run_linear_nonuniform_case():
    # Get exact-solution parameters.
    params = make_linear_case_parameters()

    # Unpack parameters.
    eps_grad, G, kappa, xeq, vgyro, A, C, vbar = params

    # Build the electric field function.
    Efunc = linear_electric_field_factory(G)

    # Initial position at maximum x excursion.
    X0 = np.array([xeq + A, 0.0, 0.0])

    # Initial centered velocity from the exact solution.
    V0_centered = np.array([0.0, C - OMEGA * X0[0], 0.0])

    # Run for an integer number of modified gyroperiods.
    tmax = 20 * 2.0 * np.pi / kappa

    # Exact velocity helper for half-step initialization.
    def exact_velocity_at_time(tt):
        _, Vtmp = linear_exact_solution(np.array([tt]), X0, params)
        return Vtmp[0]

    # Run leapfrog Boris.
    t, X, Vc = run_leapfrog_orbit(X0, V0_centered, tmax, Efunc, exact_velocity_at_time)

    # Exact solution.
    Xex, Vex = linear_exact_solution(t, X0, params)

    # For the linear case, the exact orbit itself is the cleanest theory comparison.
    # A local ExB drift estimate at xeq is also useful as a guiding-center-scale reference.
    vE_xeq = -(E0_LINEAR + G * xeq) / B0

    # Mean numerical y drift over the integer number of modified periods.
    mean_vy_num = (X[-1, 1] - X[0, 1]) / (t[-1] - t[0])

    # Error diagnostics.
    _, _, rel_x = rms_relative_error(X[:, 0], Xex[:, 0])
    _, _, rel_y = rms_relative_error(X[:, 1], Xex[:, 1])
    _, _, rel_vx = rms_relative_error(Vc[:, 0], Vex[:, 0])
    _, _, rel_vy = rms_relative_error(Vc[:, 1], Vex[:, 1])

    # Plot particle trajectory.
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.plot(X[:, 0], X[:, 1], linewidth=1.0, label="Boris particle")
    ax.plot(Xex[:, 0], Xex[:, 1], "--", linewidth=1.3, label="exact linear-field orbit")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Case 2: nonuniform static E, linear $E_x(x)$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "case2_linear_nonuniform_exb_trajectory.png", dpi=190)
    plt.close(fig)

    # Plot x comparison.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(t * 1e6, X[:, 0], label="Boris x")
    ax.plot(t * 1e6, Xex[:, 0], "--", label="exact x")
    ax.set_xlabel("time [µs]")
    ax.set_ylabel("x [m]")
    ax.set_title("Case 2: x motion in linear nonuniform E")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "case2_linear_nonuniform_exb_x_comparison.png", dpi=190)
    plt.close(fig)

    # Plot y comparison.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(t * 1e6, X[:, 1], label="Boris y")
    ax.plot(t * 1e6, Xex[:, 1], "--", label="exact y")
    ax.set_xlabel("time [µs]")
    ax.set_ylabel("y [m]")
    ax.set_title("Case 2: y drift in linear nonuniform E")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "case2_linear_nonuniform_exb_y_comparison.png", dpi=190)
    plt.close(fig)

    # Plot vy comparison.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(t * 1e6, Vc[:, 1], label="Boris centered $v_y$")
    ax.plot(t * 1e6, Vex[:, 1], "--", label="exact $v_y$")
    ax.axhline(vE_xeq, linestyle=":", color="k", label="local $E\\times B$ at $x_{eq}$")
    ax.set_xlabel("time [µs]")
    ax.set_ylabel("$v_y$ [m/s]")
    ax.set_title("Case 2: velocity comparison")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "case2_linear_nonuniform_exb_vy_comparison.png", dpi=190)
    plt.close(fig)

    # Plot a simple reference showing how much stronger the nonuniform drift is
    # than the uniform Case 1 drift.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.axhline(-E0_UNIFORM / B0, linestyle="--", label="Case 1 uniform $E\\times B$")
    ax.axhline(vE_xeq, linestyle="-", label="Case 2 local $E(x_{eq})\\times B$")
    ax.plot(t * 1e6, Vex[:, 1], alpha=0.7, label="Case 2 exact $v_y(t)$")
    ax.set_xlabel("time [µs]")
    ax.set_ylabel("$v_y$ [m/s]")
    ax.set_title("Uniform vs nonuniform $E\\times B$ speed")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "case2_uniform_vs_nonuniform_drift_speed.png", dpi=190)
    plt.close(fig)

    # Return diagnostics.
    return {
        "eps_grad": eps_grad,
        "G": G,
        "kappa_over_omega": kappa / OMEGA,
        "xeq": xeq,
        "A": A,
        "vbar": vbar,
        "local_vE_xeq": vE_xeq,
        "mean_vy_num": mean_vy_num,
        "rel_x": rel_x,
        "rel_y": rel_y,
        "rel_vx": rel_vx,
        "rel_vy": rel_vy,
        "final_x_error": X[-1, 0] - Xex[-1, 0],
        "final_y_error": X[-1, 1] - Xex[-1, 1],
    }



# ============================================================
# Case 3: magnetized polarization drift with sinusoidal E(t)
# ============================================================

def make_polarization_case_parameters():
    # Electric-field amplitude for the time-dependent field [V/m].
    Eamp = 3.0e4

    # Frequency of the imposed electric field [Hz].
    freq = 2.5e5

    # Angular frequency [rad/s].
    omega = 2.0 * np.pi * freq

    # Small parameter omega/Omega.
    eps = omega / OMEGA

    # Number of electric-field periods to simulate.
    n_periods = 8

    # Final time is an integer number of E-field periods.
    tmax = n_periods / freq

    # Return all parameters.
    return Eamp, freq, omega, eps, n_periods, tmax


def polarization_electric_field_factory(Eamp, omega):
    # Build E(t) = Eamp sin(omega t) xhat.
    def polarization_electric_field(x, t):
        # Return Ex, Ey, Ez at one particle position and time.
        return np.array([Eamp * np.sin(omega * t), 0.0, 0.0])
    return polarization_electric_field


def boris_velocity_step_time_dependent(v_half_minus, x_integer, t_integer, dt, electric_field_time_func):
    # Advance v_{n-1/2} to v_{n+1/2} for time-dependent E(t) and uniform B.
    # This is the same Boris leapfrog logic, but now the two half electric kicks
    # use E at the beginning and end of the force step.

    # Electric field at the beginning of this force step.
    E_old = electric_field_time_func(x_integer, t_integer)

    # Uniform magnetic field vector.
    B = np.array([0.0, 0.0, B0])

    # First half electric acceleration.
    v_minus = v_half_minus + QOM_D * E_old * (0.5 * dt)

    # Boris t vector.
    tvec = QOM_D * B * (0.5 * dt)

    # Boris s vector.
    svec = 2.0 * tvec / (1.0 + np.dot(tvec, tvec))

    # First magnetic rotation.
    v_prime = v_minus + np.cross(v_minus, tvec)

    # Second magnetic rotation.
    v_plus = v_minus + np.cross(v_prime, svec)

    # Electric field at the end of this force step.
    E_new = electric_field_time_func(x_integer, t_integer + dt)

    # Second half electric acceleration.
    v_half_plus = v_plus + QOM_D * E_new * (0.5 * dt)

    # Return v_{n+1/2}.
    return v_half_plus


def reconstruct_centered_velocity_time_dependent(v_half_minus, x_integer, t_integer, electric_field_time_func):
    # Reconstruct v_n from v_{n-1/2} by advancing only half a timestep.
    return boris_velocity_step_time_dependent(
        v_half_minus,
        x_integer,
        t_integer - 0.5 * DT,
        0.5 * DT,
        electric_field_time_func,
    )


def run_time_dependent_leapfrog_orbit(X0, V0_half_minus, tmax, electric_field_time_func):
    # Number of integer timesteps.
    nsteps = int(round(tmax / DT))

    # Integer-time samples.
    t = np.arange(nsteps + 1) * DT

    # Allocate position history.
    Xhist = np.zeros((nsteps + 1, 3))

    # Allocate centered velocity history.
    Vcenter_hist = np.zeros((nsteps + 1, 3))

    # Initial integer-time position.
    x = X0.copy()

    # Initial half-step velocity v_{-1/2}.
    v_half = V0_half_minus.copy()

    # Main loop.
    for n in range(nsteps + 1):

        # Store position x_n.
        Xhist[n] = x

        # Reconstruct and store centered velocity v_n.
        Vcenter_hist[n] = reconstruct_centered_velocity_time_dependent(
            v_half,
            x,
            t[n],
            electric_field_time_func,
        )

        # Stop after final sample.
        if n == nsteps:
            break

        # Advance half-step velocity from v_{n-1/2} to v_{n+1/2}.
        v_half_new = boris_velocity_step_time_dependent(
            v_half,
            x,
            t[n],
            DT,
            electric_field_time_func,
        )

        # Advance position using v_{n+1/2}.
        x = x + v_half_new * DT

        # Shift half-step velocity.
        v_half = v_half_new

    # Return integer-time diagnostics.
    return t, Xhist, Vcenter_hist


def polarization_exact_velocity(t, Eamp, omega, eps, vpar):
    # Convert t to an array.
    t = np.asarray(t)

    # Finite-frequency denominator.
    denom = 1.0 - eps**2

    # Allocate velocity array.
    V = np.zeros((len(np.atleast_1d(t)), 3))

    # Exact driven vx. Its low-frequency limit is the polarization drift.
    V[:, 0] = (Eamp / B0) * (eps / denom) * np.cos(omega * t)

    # Exact driven vy. Its low-frequency limit is the E x B drift.
    V[:, 1] = -(Eamp / B0) * (1.0 / denom) * np.sin(omega * t)

    # Keep 14 keV motion along B, so it does not contaminate transverse drift.
    V[:, 2] = vpar

    # Return exact velocity.
    return V


def polarization_exact_position(t, X0, Eamp, omega, eps, vpar):
    # Convert t to an array.
    t = np.asarray(t)

    # Finite-frequency denominator.
    denom = 1.0 - eps**2

    # Allocate position array.
    X = np.zeros((len(np.atleast_1d(t)), 3))

    # Exact x position with x(0)=X0[0].
    X[:, 0] = X0[0] + (Eamp / B0) * (eps / denom) * np.sin(omega * t) / omega

    # Exact y position with y(0)=X0[1].
    X[:, 1] = X0[1] + (Eamp / B0) * (1.0 / denom) * (np.cos(omega * t) - 1.0) / omega

    # Exact z position from parallel motion.
    X[:, 2] = X0[2] + vpar * t

    # Return exact position.
    return X


def polarization_first_order_drifts(t, Eamp, omega):
    # Convert t to an array.
    t = np.asarray(t)

    # Allocate first-order ExB velocity.
    Vexb = np.zeros((len(np.atleast_1d(t)), 3))

    # Allocate first-order polarization velocity.
    Vpol = np.zeros_like(Vexb)

    # First-order E x B drift.
    Vexb[:, 1] = -(Eamp / B0) * np.sin(omega * t)

    # First-order polarization drift.
    Vpol[:, 0] = (M_D / (QE * B0**2)) * Eamp * omega * np.cos(omega * t)

    # Return both components.
    return Vexb, Vpol


def run_polarization_case():
    # Get parameters.
    Eamp, freq, omega, eps, n_periods, tmax = make_polarization_case_parameters()

    # Use 14 keV speed parallel to B.
    vpar = speed_from_energy_kev(E_KEV)

    # Time-dependent electric field.
    Efunc = polarization_electric_field_factory(Eamp, omega)

    # Initial position.
    X0 = np.array([0.0, 0.0, 0.0])

    # Exact half-step initial velocity v_{-1/2}.
    V0_half_minus = polarization_exact_velocity(
        np.array([-0.5 * DT]),
        Eamp,
        omega,
        eps,
        vpar,
    )[0]

    # Run leapfrog Boris.
    t, X, Vc = run_time_dependent_leapfrog_orbit(X0, V0_half_minus, tmax, Efunc)

    # Exact solution.
    Xex = polarization_exact_position(t, X0, Eamp, omega, eps, vpar)
    Vex = polarization_exact_velocity(t, Eamp, omega, eps, vpar)

    # First-order drift approximations.
    Vexb_first, Vpol_first = polarization_first_order_drifts(t, Eamp, omega)

    # Error diagnostics.
    _, _, rel_x = rms_relative_error(X[:, 0], Xex[:, 0])
    _, _, rel_y = rms_relative_error(X[:, 1], Xex[:, 1])
    _, _, rel_vx = rms_relative_error(Vc[:, 0], Vex[:, 0])
    _, _, rel_vy = rms_relative_error(Vc[:, 1], Vex[:, 1])

    # First-order polarization error relative to exact finite-frequency vx.
    _, _, rel_vpol_first = rms_relative_error(Vpol_first[:, 0], Vex[:, 0])

    # First-order ExB error relative to exact finite-frequency vy.
    _, _, rel_vexb_first = rms_relative_error(Vexb_first[:, 1], Vex[:, 1])

    # Plot transverse trajectory.
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.plot(X[:, 0], X[:, 1], linewidth=1.0, label="Boris transverse path")
    ax.plot(Xex[:, 0], Xex[:, 1], "--", linewidth=1.4, label="exact finite-frequency path")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Case 3: magnetized polarization drift trajectory")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "case3_polarization_trajectory.png", dpi=190)
    plt.close(fig)

    # Plot vx comparison.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(t * 1e6, Vc[:, 0], label="Boris centered $v_x$")
    ax.plot(t * 1e6, Vex[:, 0], "--", label="exact finite-frequency $v_x$")
    ax.plot(t * 1e6, Vpol_first[:, 0], ":", label="first-order polarization $v_x$")
    ax.set_xlabel("time [µs]")
    ax.set_ylabel("$v_x$ [m/s]")
    ax.set_title("Case 3: polarization drift velocity")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "case3_polarization_vx_comparison.png", dpi=190)
    plt.close(fig)

    # Plot vy comparison.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(t * 1e6, Vc[:, 1], label="Boris centered $v_y$")
    ax.plot(t * 1e6, Vex[:, 1], "--", label="exact finite-frequency $v_y$")
    ax.plot(t * 1e6, Vexb_first[:, 1], ":", label="first-order $E(t)\\times B$ $v_y$")
    ax.set_xlabel("time [µs]")
    ax.set_ylabel("$v_y$ [m/s]")
    ax.set_title("Case 3: accompanying $E(t)\\times B$ velocity")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "case3_polarization_vy_comparison.png", dpi=190)
    plt.close(fig)

    # Plot x position comparison.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(t * 1e6, X[:, 0], label="Boris x")
    ax.plot(t * 1e6, Xex[:, 0], "--", label="exact x")
    ax.set_xlabel("time [µs]")
    ax.set_ylabel("x [m]")
    ax.set_title("Case 3: polarization drift x-position")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "case3_polarization_x_comparison.png", dpi=190)
    plt.close(fig)

    # Return diagnostics.
    return {
        "Eamp": Eamp,
        "freq": freq,
        "omega_over_Omega": eps,
        "n_periods": n_periods,
        "rel_x": rel_x,
        "rel_y": rel_y,
        "rel_vx": rel_vx,
        "rel_vy": rel_vy,
        "rel_vpol_first": rel_vpol_first,
        "rel_vexb_first": rel_vexb_first,
        "final_x_error": X[-1, 0] - Xex[-1, 0],
        "final_y_error": X[-1, 1] - Xex[-1, 1],
        "final_z_error": X[-1, 2] - Xex[-1, 2],
    }


def main():
    # Run all three cases.
    r1 = run_uniform_case()
    r2 = run_linear_nonuniform_case()
    r3 = run_polarization_case()

    # Print shared setup.
    print("Static E x B two-case benchmark")
    print(f"  E_KEV = {E_KEV:.3f} keV")
    print(f"  14 keV speed = {speed_from_energy_kev(E_KEV):.6e} m/s")
    print(f"  B0 = {B0:.6e} T")
    print(f"  Omega = {OMEGA:.6e} rad/s")
    print(f"  Tgyro = {T_GYRO:.6e} s")
    print(f"  dt = {DT:.6e} s")
    print(f"  dt/Tgyro = {DT/T_GYRO:.6e}")
    print()

    # Print uniform case.
    print("Case 1: uniform static E x B")
    print(f"  E0 = {E0_UNIFORM:.6e} V/m")
    print(f"  theory vE_y = {r1['vE_y']:.6e} m/s")
    print(f"  Larmor radius = {r1['rho']:.6e} m")
    print(f"  RMS X_gc absolute error = {r1['abs_xgc_rms']:.3e} m")
    print(f"  RMS Y_gc absolute error = {r1['abs_ygc_rms']:.3e} m")
    print(f"  final X_gc error = {r1['final_xgc_error']:+.3e} m")
    print(f"  final Y_gc error = {r1['final_ygc_error']:+.3e} m")
    print()

    # Print nonuniform case.
    print("Case 2: linear nonuniform static E")
    print(f"  E0 = {E0_LINEAR:.6e} V/m")
    print(f"  G = dEx/dx = {r2['G']:.6e} V/m^2")
    print(f"  gamma/Omega^2 = {r2['eps_grad']:.6e}")
    print(f"  kappa/Omega = {r2['kappa_over_omega']:.6e}")
    print(f"  xeq = {r2['xeq']:.6e} m")
    print(f"  oscillation amplitude A = {r2['A']:.6e} m")
    print(f"  exact mean y drift vbar = {r2['vbar']:.6e} m/s")
    print(f"  local ExB drift at xeq = {r2['local_vE_xeq']:.6e} m/s")
    print(f"  numerical mean vy = {r2['mean_vy_num']:.6e} m/s")
    print(f"  relative x error = {r2['rel_x']:.3e}")
    print(f"  relative y error = {r2['rel_y']:.3e}")
    print(f"  relative vx error = {r2['rel_vx']:.3e}")
    print(f"  relative vy error = {r2['rel_vy']:.3e}")
    print(f"  final x error = {r2['final_x_error']:+.3e} m")
    print(f"  final y error = {r2['final_y_error']:+.3e} m")
    print()

    # Print polarization case.
    print("Case 3: magnetized polarization drift")
    print(f"  Eamp = {r3['Eamp']:.6e} V/m")
    print(f"  frequency = {r3['freq']:.6e} Hz")
    print(f"  omega/Omega = {r3['omega_over_Omega']:.6e}")
    print(f"  number of E-field periods = {r3['n_periods']}")
    print(f"  relative x error = {r3['rel_x']:.3e}")
    print(f"  relative y error = {r3['rel_y']:.3e}")
    print(f"  relative vx error = {r3['rel_vx']:.3e}")
    print(f"  relative vy error = {r3['rel_vy']:.3e}")
    print(f"  first-order polarization-vx relative error = {r3['rel_vpol_first']:.3e}")
    print(f"  first-order ExB-vy relative error = {r3['rel_vexb_first']:.3e}")
    print(f"  final x error = {r3['final_x_error']:+.3e} m")
    print(f"  final y error = {r3['final_y_error']:+.3e} m")
    print(f"  final z error = {r3['final_z_error']:+.3e} m")
    print()

    # Print output folder.
    print(f"Saved plots in: {OUTDIR}")


if __name__ == "__main__":
    main()
