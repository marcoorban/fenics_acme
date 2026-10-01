"""Post-process the DFG 2D-2 cylinder benchmark (Re=100) OpenFOAM run.

Reads the per-timestep forceCoeffs and pressure-probe series written by
`foamPostProcess` (see ../plotting.md for how those were generated), plots
Cd(t)/Cl(t), and derives the pressure drop and Strouhal number from the
saturated (post-transient) portion of the run.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

CASE_DIR = Path(__file__).parent / "openFoam"
PLOTS_DIR = Path(__file__).parent / "plots"

D = 0.1  # cylinder diameter [m]
Ubar = 1.0  # mean inlet velocity [m/s] (Re=100: Um=1.5, Ubar=(2/3)*Um)

# The lift-driven instability grows from t=0 and only reaches its limit
# cycle after several shedding periods (checked by eye: Cl/Cd amplitude
# is still climbing through t~3.5s and roughly flat by t~4s). Stats and
# the Strouhal number should only use data past this point.
T_SATURATED = 4.0  # for Cd/Cl extrema and pressure drop [s]
T_STROUHAL = 3.0  # wider window for period averaging (frequency is the
# same during growth as at saturation, so more zero-crossings = better
# averaging here, unlike the amplitude stats above) [s]


def load_force_coeffs():
    data = np.loadtxt(CASE_DIR / "postProcessing/forceCoeffsDict/0/forceCoeffs.dat")
    return data[:, 0], data[:, 2], data[:, 3]  # t, Cd, Cl


def load_pressure_probes():
    data = np.loadtxt(CASE_DIR / "postProcessing/pressureDropDict/0/p")
    return data[:, 0], data[:, 1], data[:, 2]  # t, p_front, p_rear


def plot_coefficients(t, Cd, Cl):
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    axes[0].plot(t, Cd)
    axes[0].set_ylabel("$C_d$")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(t, Cl, color="tab:orange")
    axes[1].set_ylabel("$C_l$")
    axes[1].set_xlabel("Time [s]")
    axes[1].grid(True, alpha=0.3)
    fig.suptitle("DFG 2D-2 cylinder benchmark (Re=100): force coefficients")
    fig.tight_layout()
    PLOTS_DIR.mkdir(exist_ok=True)
    out = PLOTS_DIR / "re100_force_coefficients.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def strouhal_from_zero_crossings(t, Cl, t_min):
    mask = t >= t_min
    tt, cc = t[mask], Cl[mask]
    crossings = []
    for i in range(len(cc) - 1):
        if cc[i] < 0 <= cc[i + 1]:
            frac = -cc[i] / (cc[i + 1] - cc[i])  # linear interp for sub-sample crossing time
            crossings.append(tt[i] + frac * (tt[i + 1] - tt[i]))
    periods = np.diff(crossings)
    period = periods.mean()
    St = (1.0 / period) * D / Ubar
    return St, period, periods, crossings


def _refine_peak_time(t, y, i):
    """Parabolic-vertex refinement of the local max of y at index i, using
    its two neighbors -- sub-sample-accurate peak time without needing a
    finer-resolution signal."""
    t0, t1, t2 = t[i - 1], t[i], t[i + 1]
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    denom = (t0 - t1) * (t0 - t2) * (t1 - t2)
    A = (t2 * (y1 - y0) + t1 * (y0 - y2) + t0 * (y2 - y1)) / denom
    B = (t2**2 * (y0 - y1) + t1**2 * (y2 - y0) + t0**2 * (y1 - y2)) / denom
    return t1 if A == 0 else -B / (2 * A)


def delta_p_at_half_period(t, Cl, tp, dp, t_min, period):
    """DFG benchmark definition: delta_p is reported at t0 + T/2, where t0
    is a time at which Cl attains a local maximum -- not a plain time-average
    over some window. Averages the result over every Cl-maximum found past
    t_min for which t0+T/2 still falls inside the recorded pressure series."""
    mask = t >= t_min
    idx = np.where(mask)[0]
    tt, cc = t[mask], Cl[mask]
    peak_idx_local = [i for i in range(1, len(cc) - 1) if cc[i] > cc[i - 1] and cc[i] > cc[i + 1]]
    peak_times = [_refine_peak_time(t, Cl, idx[i]) for i in peak_idx_local]

    samples = []
    for t0 in peak_times:
        t_half = t0 + period / 2
        if t_half <= tp.max():
            samples.append(np.interp(t_half, tp, dp))
    return np.mean(samples), samples, peak_times


def main():
    t, Cd, Cl = load_force_coeffs()
    plot_coefficients(t, Cd, Cl)

    mask = t >= T_SATURATED
    print(f"\nSaturated-cycle stats (t >= {T_SATURATED}s, {mask.sum()} samples):")
    print(f"  Cd: mean={Cd[mask].mean():.4f}  min={Cd[mask].min():.4f}  max={Cd[mask].max():.4f}")
    print(f"  Cl: min={Cl[mask].min():.4f}  max={Cl[mask].max():.4f}")

    St, period, periods, crossings = strouhal_from_zero_crossings(t, Cl, T_STROUHAL)
    print(f"\nStrouhal number (t >= {T_STROUHAL}s, {len(crossings)} zero-crossings, "
          f"{len(periods)} periods):")
    print(f"  periods: {np.array2string(periods, precision=4)}")
    print(f"  mean period T={period:.4f}s -> St = T^-1 * D / Ubar = {St:.4f}")

    tp, p_front, p_rear = load_pressure_probes()
    dp_series = p_front - p_rear

    dp_half, dp_samples, peak_times = delta_p_at_half_period(
        t, Cl, tp, dp_series, T_SATURATED, period
    )
    print(f"\nPressure drop at t0+T/2 (DFG benchmark definition, t0 = Cl local max, "
          f"t0 >= {T_SATURATED}s):")
    print(f"  Cl-maxima used as t0: {np.array2string(np.array(peak_times), precision=4)}")
    print(f"  delta_p samples: {np.array2string(np.array(dp_samples), precision=6)}")
    print(f"  delta_p = {dp_half:.6f}")


if __name__ == "__main__":
    main()
