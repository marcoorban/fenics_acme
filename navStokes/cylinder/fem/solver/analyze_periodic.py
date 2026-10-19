"""Post-process a coefficients.dat (t,C_D,C_L,delta_p) time series for a
periodic benchmark case (DFG 2D-2 convention, Schafer & Turek 1996):
  - C_D_max, C_L_max: peak values over the recorded series.
  - t*: the time at which C_L is maximum.
  - delta_p: evaluated at t* + period/2 (half a shedding period after peak
    lift), not at t* itself -- the benchmark's own definition for this
    case. period = D / (St * Ubar) from the given Strouhal number.
  - "at t* + period/2" means the recorded sample closest to that time,
    since coefficients.dat is only sampled every write_every steps, not
    every dt.

Usage: python analyze_periodic.py <coefficients.dat> <St> <D> <Um>
"""

import csv
import sys


def main():
    coeffs_path = sys.argv[1]
    St, D, Um = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])

    with open(coeffs_path) as f:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]

    cd_max = max(r["C_D"] for r in rows)
    cl_row = max(rows, key=lambda r: r["C_L"])
    cl_max = cl_row["C_L"]
    t_star = cl_row["t"]

    Ubar = (2 / 3) * Um
    period = D / (St * Ubar)
    t_eval = t_star + period / 2

    closest = min(rows, key=lambda r: abs(r["t"] - t_eval))

    print(f"{cd_max},{cl_max},{t_star},{period},{t_eval},{closest['t']},{closest['delta_p']}")


if __name__ == "__main__":
    main()
