# Analysis macros, 2026-09-13

Every number recorded in `RecoTracker/CLAUDE.md` under the 2026-09-13 entries,
and every plot in `~/mic-dev/validation-slides`, was produced by one of these.
They were written in a session scratchpad, which is temporary — moved here so
the findings have a provenance that outlives it.

## How to run

Analysis needs only the dictionary library, not the mkFit build:

    . /home/matevz/root7.env ; unset DISPLAY
    cd /foo/matevz/mic-dev/current/src/standalone
    ROOT_INCLUDE_PATH=/foo/matevz/mic-dev/current/src \
    LD_LIBRARY_PATH=. root.exe -l -b -q \
      '../RecoTracker/MkFitCore/standalone/test/an-2026-09-13/<macro>.C("<input>.root")'

**Generation** needs the mkFit build and the `val_*` functions in
`RdfTrace/ValProp.cc`, driven from `--shell`. See CLAUDE.md for the per-study
invocation; the driver entry points are in `test/val-prop.C`.

## Macro -> finding

| macro | input | what it produced |
|---|---|---|
| `prof`, `tot`, `tots`, `tots2` | `val-bkfit*.root` | backward-fit chi2 profile, the seed-covariance scale scan, material on/off |
| `an-covx`, `an-norm`, `an-cv`, `an-chain` | `val-covxport.root` | S11 step 1: ensemble transport, rank, the curvilinear projection, the full chain |
| `an-plane`, `an-inc` | `val-covxport.root` | the decomposition-free plane-distance check; whether the incidence was applied |
| `an-herm`, `an-wx` | `val-search*.root` | Hermite vs the exact double solve; the wrong-crossing cases |
| `an-search`, `an-lay`, `an-res`, `an-cut` | `val-search-mat1.root` | inward-search funnel, per-layer q pulls, module-frame residuals, chi2 vs the production cut |
| `an-mat`, `an-fpix` | `val-search-mat{0,1}.root` | material on/off in the search; the TFPX/TEPX split |
| `an-seq`, `an-doom2`, `an-seed` | `val-search-mat1.root` | per-layer-search and per-seed funnels (NOTE: survivorship-biased, see CLAUDE.md) |
| `an-miss`, `an-pix` | `val-miss-pt5c.root` | why the true hit was not considered; the pixel-region verdict split |
| `an-beta` | `val-miss-pt5c.root` | best-hit efficiency vs \|eta\| — the 97.5 % -> 62.4 % collapse |
| `an-cut2` | `val-miss-pt5c.root` | which half of the pre-selection cut rejects — dq, 98 % at \|eta\| > 2 |
| `an-dqterm` | `val-miss-pt5c.root` | the dq cut decomposed: the track term's 9x collapse |
| `mkplots`, `runplots`, `mk350`, `split250`, `split340` | various | the deck's `data/*.root` histograms |

## Traps these macros embody

- **`search_id` is PER EVENT** — key any per-search map on `(event, search_id)`.
- **Group by the sim track's own hits** (`n_sim_hits_in_layer`), not by whether a
  hit was scanned, or the inclusive layer plan is counted as a failure.
- **Per-layer-search pickup is survivorship-biased**: a dead candidate makes no
  searches. Use the track-level macros for efficiency.
- **Never spot-check a density table by SAMPLING it -- integrate.** Material
  lives in thin dense spikes separated by large near-empty gaps, so a handful of
  sample points lands in the gaps essentially always. Sampling `matdiffuse_rz`'s
  `diffuse_x0` at three radii gave 3.33e-5/cm each, which multiplied out to
  ~1 % of the total and killed a correct plan for a day; integrating the same
  table gives 0.0551 at z = 0 (r 3-110) and 0.11-0.23 through the IT forward
  region, i.e. 15x more. The file's own units line ("density per cm, multiply by
  path length") invites the multiplication -- it is still wrong.
- **`core()` divides the IQR by 1.349, which is the GAUSSIAN factor.** For the
  ALONG-STRIP coordinate the hit error is **uniform**, not Gaussian: on [-a, a]
  the IQR is `a` while sigma is `a/sqrt(3)`, so IQR/sigma = **1.732** and using
  1.349 over-estimates sigma by **1.28x**. Wherever the uniform hit term
  dominates the residual, the measured "core" is 28 % too big and the covariance
  looks worse than it is -- a systematic in exactly the direction that
  manufactures a fake deficit.

  **Where it is safe, and for the right reason.** PHI IS ALSO UNIFORM in strips
  -- a single-strip cluster is uniform over the pitch, and charge interpolation
  only helps for multi-strip clusters, so "phi is Gaussian" is NOT the argument.
  What saves phi in the outer barrel is that the term is SMALL: ~26 um against
  ~345 um of track residual is 7 % in sigma and **0.5 % in variance**, so the
  convolution is overwhelmingly the Gaussian track term and the hit shape cannot
  matter whatever it is. In the PIXELS both legs hold: clusters are 3x2 and
  larger, so charge interpolation really does make the error Gaussian-ish --
  which is why `hl_fac = 3`, a 3-sigma convention, is the right treatment there
  and `sqrt(3)` the right one for a uniform strip.

  **MEASURED, `val_cluster_sizes()`** (rows = local x = precise/phi, cols =
  local y = coarse/q), one event:

  | region | <rows> | rows==1 | <cols> | cols==1 |
  |---|---|---|---|---|
  | PixB 0-3 | **3.27** | 18.9 % | **2.81** | 29.6 % |
  | TBPS-P 4/6/8 | 2.47 | 36.6 % | 1.00 | **100 %** |
  | TBPS-S 5/7/9 | 2.82 | 29.1 % | 1.00 | **100 %** |
  | TOB 2S 10-15 | 3.34 | 23.7 % | 1.00 | **100 %** |
  | fwd pix 16-27 | 2.18 | **54.1 %** | 1.40 | **72.5 %** |
  | TEC 28-37 | 2.42 | 36.5 % | 1.00 | 100 % |

  - **Pixel barrel is 3.27 x 2.81**, multi-cell in BOTH directions, so charge
    interpolation applies and `hl_fac = 3` is right. Numbers there are safe.
  - **Strips: cols are 100 % single-cell everywhere** -- pure uniform along the
    strip, exactly what `sqrt(3)` assumes. But ROWS average 2.5-3.3 with only
    23-37 % single, so PHI in strips is mostly charge-interpolated and
    Gaussian-ish after all.
  - **The FORWARD DISKS are the worst place for this, not the best**: 54 %
    single in phi and **72.5 % single in the coarse direction**, because tracks
    cross a disc near-normally and share little charge. So disk hit errors are
    uniform-dominated (IQR bias up to 1.28x) AND the endcap `half_length` uses
    `sqrt(exx + eyy)` -- the TRACE, not a marginal -- returning 4.24 sigma where
    3 was intended for a pixel with sigma_x ~ sigma_y. Two systematics of order
    1.3x and 1.4x on top of the ~1.2x effect being looked for.

  **Consequence: q is only cleanly measurable in the PIXEL BARREL.** Outside it,
  use phi.
- **`hit_q_half_len` carries a different factor for pixels and strips**:
  `hl_fac = 3` for pixels (a Gaussian 3-sigma convention) and `sqrt(3)` for
  strips (the uniform sigma), `HitStructures.cc`. Dividing by the wrong one
  mis-states sigma_hit by 1.7x; it cost one wrong table on 2026-09-17.
- **Do not rank-detect a covariance on raw eigenvalues** — it mixes cm^2, rad^2
  and GeV^-2. Scale to the correlation matrix first (`an-covx`).
