# Plotting/post-processing a finished run

How to go from a completed OpenFOAM run (all time directories
reconstructed) to force-coefficient/pressure plots plus the derived
Strouhal number and pressure drop. Worked example: DFG 2D-2 (Re=100).
Full numbers live in `results.csv` and `log.md` — this file is the
method, not the results.

## 1. Get Cd/Cl/Cm and p at every saved time step

`foamRun.fish` only calls `foamPostProcess ... -latestTime` — fine for
a single final-time snapshot, useless for a time history. Re-run the
same function dicts *without* `-latestTime`: with no `-time`/`-latestTime`
flag, `foamPostProcess` walks every time directory found on disk and
appends a row per step to the function object's `.dat` file.

```
foamPostProcess -solver incompressibleFluid -func forceCoeffsDict
foamPostProcess -solver incompressibleFluid -func pressureDropDict
```

Gotchas (see `log.md` for the fuller writeup):
- Needs the full field reconstruction done first (`reconstructPar`,
  not `-latestTime`), one time directory per write.
- Takes noticeably longer than the solve itself for many saved steps —
  it's re-reading every time directory, not just the last one.
- Output lands at `postProcessing/<funcName>/<startTime>/<file>`:
  `forceCoeffs.dat` (columns `Time Cm Cd Cl Cl(f) Cl(r)`), probes'
  raw `p` file (columns `Time p_probe0 p_probe1 ...`, header comments
  giving each probe's coordinates).

## 2. Find the saturated (post-transient) window first

For an unsteady/periodic case, don't compute anything periodic
(Strouhal, amplitude, mean Δp) from `t=0` — the wake instability grows
from a near-symmetric start and takes several shedding cycles to reach
its limit cycle. Check by eye first:

```
for each 0.5s window:
    print min/max of Cl (and Cd) in that window
# amplitude climbing -> still transient; flat -> saturated
```

Use the first window where consecutive per-cycle amplitudes stop
changing as the cutoff for anything amplitude-sensitive (Cd/Cl
extrema, Δp). Frequency itself is unaffected by the amplitude still
growing, so a wider (earlier-starting) window is fine — and better,
since it gives more cycles to average — when only the period/Strouhal
number is needed.

## 3. Pressure drop

Just the difference of the two probe columns, evaluated per row:

```
dp[i] = p_front[i] - p_rear[i]
```

For a steady case this converges to one number. For a periodic case
report mean/min/max over the saturated window instead of a single
value — it oscillates with the shedding cycle same as Cd/Cl do.

## 4. Strouhal number, from positive-going zero-crossings of Cl(t)

FFT works but a short, coarsely-sampled series (here: 5s @ 0.02s
writeInterval) gives poor frequency resolution. Zero-crossing timing
with linear interpolation between the bracketing samples is more
precise for this case:

```
for consecutive samples (t_i, Cl_i), (t_{i+1}, Cl_{i+1}):
    if Cl_i < 0 <= Cl_{i+1}:
        frac = -Cl_i / (Cl_{i+1} - Cl_i)
        crossing_time = t_i + frac * (t_{i+1} - t_i)
period = mean(diff(crossing_times))
St = (1/period) * D / Ubar
```

Average over several periods (all crossings in the saturated-or-wider
window), don't just take one cycle.

## 5. Plots

Two complementary tools, matching the project's `#plotting` convention
(`log.md`) of Okabe-Ito colors / `pngcairo` / italic axis labels:

- **`plots.gp`** (gnuplot) — quick, declarative, matches every other
  case's validation figures. Reads the `.dat`/probe files directly
  (awk one-liner for the pressure-drop column subtraction, since
  gnuplot can't do arbitrary column math inline). Good for the
  zoomed-in, periodic-tail figures with a reference line for the
  published benchmark value.
- **`plot_results.py`** (numpy, no scipy needed) — for anything gnuplot
  can't do inline: the zero-crossing interpolation above, windowed
  amplitude stats, and a full `0`–`endTime` overview plot. Run with
  `uv run python3 cylinder/plot_results.py` (numpy/matplotlib are
  project deps in `pyproject.toml`).

Pick gnuplot when it's just "plot this column vs that column";
reach for the python script as soon as you need to compute something
(a crossing time, a windowed min/max, an FFT) rather than just draw it.
