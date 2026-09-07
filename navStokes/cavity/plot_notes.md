# plots.gp — styling notes (2026-09-07)

Backup of the original script: `plots.gp.bak`. Rendered `u.png`/`v.png`
(since renamed to `plots/u-re100.png`/`plots/v-re100.png` in the project
reorg) with `gnuplot plots.gp` to confirm the new version actually runs
and looks right before writing this up — no data or column-selection
changes, styling only. (Confirmed while doing this: `using 1:2` on the
Ghia tables already picks the Re=100 column, matching this case's
`constant/physicalProperties` — that mapping wasn't touched.)

Note: this file and `plots.gp.bak` briefly ended up stranded at the old
`navStokes/openFoam/cavity/` path after the reorg (they were untracked
when the directory move happened, so didn't move with it); moved back
here alongside `plots.gp` afterward. `plots.gp.bak` itself was found in
the OS Trash at that point — the original `plots.gp` was never committed
to git, so that trashed copy was the only surviving copy of the
pre-styling version.

## Changes made

- **Reference vs. present-result convention**: Ghia's data (17 discrete
  tabulated points) is now drawn as **markers only** (`with points`,
  open circle `pt 6`); the OpenFOAM profile (200-point `sampleDict`
  line) is drawn as a **line only** (`with lines`). Previously both were
  `with line`, which is misleading for Ghia's data since it implies a
  continuous measurement rather than 17 discrete grid values connected
  by straight segments.
- **Color palette**: black markers + blue (`#0072B2`) line — from the
  Okabe-Ito colorblind-safe palette, standard in a lot of physics/CFD
  journals. Previously gnuplot's default palette (unset, so default
  red/green cycling) which isn't colorblind-safe and reads as less
  "serious."
- **Font/terminal**: `pngcairo enhanced font "Helvetica,14"` (was: no
  font set, default gnuplot bitmap font, which looks noticeably
  low-resolution/unprofessional at presentation size). `enhanced` turns
  on markup so axis labels can be italicized.
- **Axis labels italicized**: `{/Italic y}`, `{/Italic u}`, etc. —
  physical-quantity symbols in italics is the standard convention
  (ISO 80000 / most journal style guides) that plain-text labels don't
  follow.
- **Border/tics**: `set tics in`, minor tics added (`set mxtics 2;
  set mytics 2`), border linewidth bumped to 1.2. Reads as a tighter,
  more deliberate figure than gnuplot's default outward-pointing tics
  with no minor gridlines.
- **Light gridlines**: very pale gray (`#e0e0e0`/`#eeeeee`) on both
  major and minor tics — aids reading values off the curve without
  the plot looking cluttered. This is a judgment call — see Options
  below.
- **Legend box**: `box opaque` so the legend has a visible border and
  doesn't let the curve show through behind it; explicit `top left` /
  `top right` placement chosen by eye to avoid overlapping the curves
  for this specific data shape.
- **`Re = 100` label** added directly on each figure (`set label ... at
  graph 0.05, 0.90`) — common in validation-figure convention so the
  case parameter is visible without relying on a caption.
- **Citation fix**: legend read `"Ghia (1981)"` in the u-plot (wrong
  year — the file itself and its header comment both say 1982) and
  plain `"Ghia"` in the v-plot (inconsistent with the other one). Both
  now read `"Ghia et al. (1982)"`, consistently.
- **`set xrange [0:1]`** added explicitly on both plots (was
  auto-scaled) — pins the physical domain extent exactly rather than
  leaving gnuplot's default small auto-margin outside \[0,1\].

## Options considered but *not* applied — flagging for your review

- **Vector output (PDF/EPS) instead of PNG.** Journals generally want
  vector graphics for line plots (infinite resolution, smaller file
  size, crisp text when embedded in LaTeX) rather than raster PNG.
  I left the terminal as `pngcairo` since that's what the original
  script used and swapping output format is a bigger change than
  "styling." If you want this, the swap is just:
  ```gnuplot
  set terminal epslatex size 5,3.75 color colortext standalone
  # or: set terminal pdfcairo enhanced font "Helvetica,14" size 5,3.75
  set output "u.tex"   # / "u.pdf"
  ```
  `epslatex` is the more "publication-native" option if the final
  destination is a LaTeX document (it typesets the labels in your
  document's own font via `\input{u.tex}`), `pdfcairo` is simpler if
  you just want a good-looking standalone PDF.
- **Grid on/off.** Added light gridlines (see above) — some journals
  (and some advisors) prefer *no* gridlines at all for a cleaner look,
  relying only on tics. Easy to remove: delete the `set grid ...` line.
- **Legend position.** Chose `top left` (u-plot) / `top right` (v-plot)
  by eye against the current data. If the curve shape changes much
  (different Re, different mesh), these may need to be moved — no
  automatic "avoid the data" placement in gnuplot.
- **Aspect ratio / figure size.** Kept the original `size 1200,900`
  (bumped from `1200,800` to give the added `Re=100` label headroom).
  Journals often want a specific aspect ratio (e.g. 4:3 or golden
  ratio) tied to their column width — not set here since I don't know
  the target publication's requirements.
