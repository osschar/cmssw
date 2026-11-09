# MkFinderV2p2 design notes

Design decisions behind `MkFinderV2p2` and the classes it uses (`MkBins`,
`LayerOfHits::HitInfo`, `V2p2Score.h`, `V2p2Config.h`). The code comments say
what a piece of code does and point to a section here for why. For what the
code does as a whole, read the overview at the top of `MkFinderV2p2.cc` first.
Each measurement quoted here names its sample, its size and what was counted.

Phase-2 layer names: the OT PS barrel is OT layers 1-3 (mkFit layers 4-9, one
mkFit layer per sensor of a module), the OT 2S barrel is OT layers 4-6 (mkFit
10-15), and the OT endcap has PS modules at small radius and 2S at large.

Contents:

1. [Samples](#samples)
2. [Configuration](#configuration)
3. [Layer processing](#layer-processing)
4. [Candidate pickup and stopping cuts](#candidate-pickup-and-stopping-cuts)
5. [Layer crossings and the Hermite cubic](#layer-crossings-and-the-hermite-cubic)
6. [Reachability](#reachability)
7. [Within-sensitive-region verdict](#within-sensitive-region-verdict)
8. [Search window](#search-window)
9. [Hit extents](#hit-extents)
10. [Line pre-cut](#line-pre-cut)
11. [Reduction and hit ordering](#reduction-and-hit-ordering)
12. [In-layer combinatorial search](#in-layer-combinatorial-search)
13. [Kalman hit acceptance](#kalman-hit-acceptance)
14. [End-of-layer selection](#end-of-layer-selection)
15. [Layer-step score](#layer-step-score)
16. [V2 compatibility](#v2-compatibility)

## Samples

All measurements use ttbar events with pileup 200 in the phase-2 tracker
geometry. The forward search runs from the CMSSW initialStep seeds, the inward
search from T5s from LST and has to add their pixel hits. The samples used
repeatedly are named here and by tag afterwards.

- **F30**: forward search, 30 events. Efficiency against the MTV-like
  denominator (29759 selected sim tracks). Comparisons are paired: same events,
  same seeds, one process.
- **I50**: inward search, 50 events, T5s from LST. Counted: pixel hits
  recovered, out of 211435, matched by exact (layer, index); and fully recovered
  tracks, those with all of them, out of 43358. Neither count uses truth
  matching.
- **T100**: timing of the forward search, 100 events, single thread, build
  without ROOT and without tracing, minimum of three repetitions. The quoted
  time brackets `(builder.*FindTracks)(IT_FwdSearch)` only.
- **H50**: the events of I50, used for detector properties such as hit
  multiplicities per layer crossing.

## Configuration

v2p2 has two kinds of parameters.

| | where | set by |
|---|---|---|
| global knobs: search window, pre-cut, in-layer search, stopping policy, score form | `src/V2p2Config.h`, namespace `Config::V2p2` | the code: constants, fixed at compile time |
| per-iteration tunings: the two score sets, the hit-acceptance cut, the hit caps | `interface/V2p2Params.h`, `IterationConfig::m_v2p2` | the python parameter set `v2p2` of `MkFitIterationConfigESProducer`; not the JSON |

The defaults of both are in the code; the python configuration only overrides
them. The global knobs are constants so that there is no mutable global state
and switched-off code compiles away. When an iteration needs its own value of a
knob, the knob moves to the per-iteration tunings, as the hit-acceptance
parameters did.

Outside CMSSW, in the standalone build, the global knobs are variables, and
both kinds can be set from the `mkFit` command line for scans. Diagnostics (the
`Diag` knobs and the policy counters) exist only there.

## Layer processing

The overview at the top of `MkFinderV2p2.cc` follows one seed through one
layer. This section explains how the code is organised around that, and why.

**Full vector lanes.** The track math runs on Matriplex batches of NN lanes
(8 with AVX, 16 with AVX-512), one track per lane. One seed rarely brings NN
candidates into a layer, so candidates from many seeds enter the layer
together and go through propagation, search window and hit pre-selection NN at
a time. Pre-selection uses a second batch shape: candidates bring different
numbers of hits, so it runs one lane per (candidate, hit) pair, filled across
candidates.

**Kalman updates across batches.** In the in-layer search each accepted hit can
be followed by more hits in the same layer, so the number of Kalman updates a
candidate needs is not known in advance (see "In-layer combinatorial search").
The updates are therefore queued, and a Kalman batch runs only when NN updates
are waiting, or at the end of the layer. The first hit of a path in the layer
and the later ones are queued separately. The first can use the crossing with
the module plane found in pre-selection; a later one has to find the crossing
from the state after the previous hit. Lanes are independent, so the result is
the same as running each step of each batch on its own, up to rounding.

Measured on T100 (phi3): the filled lanes went from 6.19 to 7.70 of 8 for first
hits and from 4.38 to 6.99 of 8 for later hits. The forward search became
5.3 % faster with `-mavx`, 4.6 % with AVX2 and 3.8 % with AVX-512.

**Finishing a seed early.** Because batches mix seeds, the work of one seed
ends at different times. For each seed the finder counts the candidates still
waiting for pre-selection, and the Kalman updates and paths still pending. When
both counts are zero, the seed's paths compete and its new candidates are
chosen right away instead of at the end of the layer, and the memory of its
paths is reused. This changes neither the result nor the time measurably. It
lowers the number of path nodes held in memory per layer from 50 to 12 on
average, at 184 bytes a node.

## Candidate pickup and stopping cuts

**Code:** `awaken_candidates()`, `begin_next_Ccrep_in_layer()`,
`stop_cuts_at_pickup()`, `fake_hit_index()`, `end_layer()`.

**Pickup.** A dormant CombCandidate wakes up when the plan reaches its pickup
layer, matched against either sub-layer of a paired entry. A seed whose last hit
sits in the layer being searched needs no special case. Its state is at that
hit, so the hit sits at dalpha = 0 and everything further along at dalpha > 0,
and the forward-only expansion never returns to it. `propagate_to_r()` picks the
crossing with the smaller |alpha|, so the entry crossing of such a candidate is
behind it, at negative dalpha, and the window also covers the part of the shell
behind the candidate. That costs scanned hits and gives no wrong answers.

**Stopping cuts.** Both cuts read only the candidate's own state, so they are
applied at pull-in, before any work is spent on the layer:

- `pT < minPtCut` from `IterationParams`;
- the looper stop, forward search only: pT below `Policy::looper_max_pt`, r above
  `Policy::looper_min_r`, and the transverse angle between position and momentum
  past `Policy::looper_max_angle` (pi/2 - 0.2, 78.5 deg). Past that angle the
  track crosses modules at grazing incidence, clusters get wide, and further
  hits degrade the measurement.

The angle test is on |posPhi - momPhi| with both angles wrapped to (-pi, pi]. An
angle A past the limit therefore appears as dphi > A, or as dphi < 2 pi - A when
the pair straddles the branch cut. The upper bound is derived from the lower one
in code. V1 writes it as the literal 4.512, which is pi + A and lets the
78.5-101.5 deg band through whenever the pair straddles. The test as written is
equivalent to cos(momPhi - posPhi) < sin(0.2).

A stopped candidate records why, in `TrackCand::StopReason_e`, because loopers
are of interest downstream (the phase-2 timing layer and HGCal). The trace
carries it for the final candidate of the forward search,
`TrCandMeta::fwd_stop_reason`; the export to CMSSW is not done yet.

**Holes.** A candidate that takes no hit in a layer records a hole HoT, chosen by
`fake_hit_index()` in V1's order: the hole limits (`maxHolesPerCand`,
`maxConsecHoles`) decide between a miss and a stop, and then the WSR overrides
with an edge or gap HoT, which do not count against the limits.

**Retiring.** `end_layer()` marks a CombCandidate finished once every TrackCand
under it has stopped.

## Layer crossings and the Hermite cubic

**Code:** `MkBins::prop_to_limits_in_order()`, `prop_to_layer_edges()`,
`preselect_hit_batch()`, `mini_propagators::Hermite3D` and `Hermite3DOnPlane`.

The candidate is propagated to the layer's two bounding surfaces in the order
the track meets them: `m_sp1` is the entry crossing and `m_sp2` the exit. These
are parameter-only propagations (`mini_propagators`), closed form in both r and
z. `m_isp` is left at `m_sp1`.

A cubic Hermite curve through the two crossings, with their momenta as
tangents, is the trajectory model inside the layer. For each scanned hit,
`Hermite3DOnPlane` solves the cubic onto the hit's own module plane with one
Newton step, and the dq/dphi cut and the state handed to the Kalman update are
both taken from that point (`h3_state`). The full helix-plane intersection is
transcendental for a tilted module, and the cubic is what makes a per-hit plane
solve affordable. The residual distance to the plane after the step is recorded
in the trace as `d_plane_h3`.

`MKFIT_TRACE_PROP_COMPARE` keeps a straight-line step onto the same plane as a
cross-check.

## Reachability

**Code:** `mini_propagators::propagate_to_r()`, `kReachMargin`, `kReachTolerance`.

In the transverse plane the track is a circle of radius R_c whose centre is at
distance d_c from the beam line, so it reaches only radii in
[|d_c - R_c|, d_c + R_c]. At both ends of that band the track grazes the
cylinder (dr/ds = 0) and the crossing solve degenerates. The target radius is
therefore kept `kReachMargin` = 0.1 cm inside the band. A target outside it is
clamped to the band edge, which returns the point of closest approach on the
trajectory, and sets the lane's `fail_flag` to 2. `fail_flag` = 1 means the solve
landed further than `kReachTolerance` from the target, which only round-off
should cause.

The margin is provisional. It was chosen on 10 events of the inward search
from T5 seeds: 0.1 cm catches 138 of 148 real failures and vetoes
none of 2087 healthy searches, 0.25 cm vetoes 7 and 1 cm vetoes 47. One length
for the whole detector is probably the wrong shape. Whether a track reaches a
layer is a statement about a track with an error, so the band should rather
scale with the track's expected dq, per layer. A retuning needs outward
searches, pixel seeds and the full pT range.

## Within-sensitive-region verdict

**Code:** `MkFinderV2p2::determine_wsr()`, `CandIncoming::m_wsr`.

The verdict is per candidate and per layer, and asks whether the segment the
track cuts through the layer's bounding shell lies in sensitive material. V1
tests one propagated point; v2p2 has both crossings and tests the segment.

- The two mini-propagator fail flags answer the radial half. Both failing means
  the track turns around before the layer: outside. Only the exit failing means
  it enters and turns around inside: at best an edge. Only the entry failing
  cannot happen, because the entry surface is the near one. In the endcap
  `propagate_to_z()` cannot fail, and the q test below carries the whole
  decision.
- `m_q_min` and `m_q_max` are the q extent of the segment. They are compared with
  the layer's q limits, widened by `Policy::wsr_n_sigma` sigma of the track's q
  error. Inside requires the whole widened segment to be inside. Erring towards
  edge is the safe side: an edge that was really inside loses one hole count,
  while an inside that was really an edge charges a hole against a candidate
  that crossed nothing.
- An endcap disc with an r hole gives a gap verdict when the segment touches the
  hole. No phase-2 layer sets `set_r_hole_range()` today, so this branch does not
  fire.

A candidate outside the layer is skipped: no bins are walked and no HoT is
added. The layer plans are deliberately inclusive (a transition plan is the
union over tracks), and without this every one of those layers would record a
hole. Measured on 20 events of the forward search: the verdict declines 64.5 %
of OT-barrel layer searches and 23.8 % of forward ones, and the sim track has a hit in the
declined layer in 0.05 % and 0.30 % of those.

## Search window

**Code:** `MkBins::determine_bin_windows()`, `MkBins::find_bin_ranges()`,
`MkFinderV2p2::preselect_hit_batch()`, `surface_referenced_dq()`,
`Config::V2p2::Window`.

The pre-selection window has two parts: a binnor **fetch**, a range of (phi, q)
bins, and a per-hit **cut**. The invariant is that the fetch covers everything
the cut accepts. A hit the cut would accept but the fetch never pulled is lost
silently.

### The cut

```
ddq   < dq_trk_fac   * dq_trk     + dq_hit_fac   * hit_q_half_length
ddphi < dphi_trk_fac * dphi_track + dphi_hit_fac * hit_phi_half_extent
```

`dq_track` and `dphi_track` are 3 sigma of the track's position error, from the
position block of the covariance at `m_sp2` (see "Track covariance at the
layer" below). `dq_trk` is `dq_track` referenced to the hit's module surface.
Each residual is taken at the hit's own module plane.

### Track covariance at the layer

`MkBins::transport_position_cov()` transports the position block of the
covariance from the previous hit to `m_sp2`, at the fixed path length the
mini-propagator reached it with. Only the position block is needed. At fixed
path length s the position is x = x0 + f(ipt, phi, theta; s), so

```
C_pos(s) = J C0 J^T,   J = [ P_in | dx/d(ipt, phi, theta) ]   (3 x 6)
```

with the derivatives in closed form from quantities the mini-propagator already
holds. With the helix in the turning angle a, k = 1/inv_k and p in GeV,

```
x = x0 + k (px sin a - py (1 - cos a)),  y = y0 + k (py sin a + px (1 - cos a)),
z = z0 + s cos(theta),                   a = s sin(theta) ipt inv_k
```

so, with (dx, dy) the displacement and p_end the transverse momentum at `m_sp2`,

```
d(x,y)/dphi   = (-dy, dx)
d(x,y)/dipt   = (-(dx, dy) + a k p_end) / ipt = k (px f1 - py f2, py f1 + px f2) / ipt
d(x,y)/dtheta = a k p_end cot(theta),     dz/dtheta = -a k / ipt
f1 = a cos a - sin a,   f2 = a sin a - (1 - cos a)
```

with p at the previous hit in the second form of the ipt column. As in
`errPropFromPathL_impl()` the result is curvilinear at both ends: `P_in`
projects the starting position onto the plane normal to the momentum there,
and the result is projected onto the plane normal to the momentum at `m_sp2`.
The dq surface reference below depends on that. The first form of the ipt
column is a difference of two O(s) terms that nearly cancel for a stiff track:
it lost 25 % in sigma_y at pT 137 GeV. The second form is used, with f1 and f2
from their series below |a| = 0.25, where they are O(a^3) and O(a^2).

No material enters. Material added during a step changes only the angular
terms, so none of it reaches the position block within that step.

This replaced a full `propagateHelixToPlaneMPlex` call, of which only the same
six elements were read. Against that call with uniform B, over 200508 lanes, the
sigmas agree to 7e-5 and the correlations to 1.6e-4 at the worst lane. Against
it with the parametric field at the starting point they differ by 0.4 % at p99.
The transport uses the same uniform B as `m_sp1`, `m_sp2` and the Hermite cubic.
On 20 events the quality-val track counts were identical in the default and two
varied configurations. The block costs 0.23 s against 0.68 s over those 20
events in a trace build, about 5 % of the build time.

**Two factors per coordinate.** Which term binds is a property of the layer.
The median `hit_q_half_length` (H50) is 0.0075 cm in the pixel barrel, 0.042 cm
on the P sensors and 0.80 cm on the S sensors of the OT PS barrel, and 2.51 cm
in the OT 2S barrel,
while the track term does not vary like that. Strips are containment-dominated
and pixels are covariance-dominated. One factor over both terms could not sit
in the right place for both, and a scan of it could not say which term it was
moving.

- `dq_hit_fac` multiplies the hit's own half-extent, so its floor is exactly
  1.0: below it the window no longer reaches the whole strip. The default 1.2
  keeps 20 % margin. Against 1.8, F30 gives -3 found tracks and -15 fakes, and
  I50 gives +55 recovered hits and -2 fully recovered tracks.
- `dq_trk_fac` is 1.5, in units of `dq_track`.

**dq surface reference.** The covariance is transported to a fixed path length,
not to a plane. The resulting q error describes where
the track is after travelling a distance s, not where it crosses the module. The
two differ by the ds degree of freedom. Sliding each point along the momentum
until it meets the surface is a linear map on the position block,

```
dx_s = (I - p^ n^T / (n^.p^)) dx
```

so the q variance becomes v^T C v with v = e_q - ((e_q.p^)/(n^.p^)) n^.
`surface_referenced_dq()` applies it with n^ the hit's own module normal. For a
radial track and a cylinder normal this amplifies sigma_q by 1/sin^2(theta).
With the module normal, the correction is almost nil for the OT PS barrel
modules, which face the interaction point, and tilted and flat layers use the same formula with no
branch. Measured on 20 events of the inward search, pixel barrel, as the ratio of the
measured residual core to the quoted sigma_q:

| \|eta\| | 0-0.8 | 0.8-1.6 | 1.6-2.0 | > 2.0 |
|---|---|---|---|---|
| without the reference | 1.64 | 4.66 | 17.0 | 46.0 |
| with the reference | 1.38 | 1.43 | 1.71 | 2.14 |

`MkBins::surface_reference_dq()` is the same correction with the layer's
cylinder or disc normal. It is off by default (`Diag::mkbins_surface_q`) and is
kept because it is the way a corrected `dq_track` reaches the trace.

**dphi jacobian.** `sigma_phi` comes from the xy covariance through
J = (-y, x)/r^2, so |J| = 1/r. One window covers the whole layer, so J is taken
at whichever crossing has the smaller radius, where `sigma_phi` is largest. The
ratio rout/rin is 1.33 at pixel layer 0.

**Per-hit phi extent and `dphi_trk_fac` = 2.** The hit term is the hit's own phi
extent from its covariance ([Hit extents](#hit-extents)), about 4e-5 rad in
the OT 2S barrel. The cut is therefore carried by the track term, and the track term needs
a factor 2:

| `dphi_trk_fac` | F30 found | F30 fakes | I50 fully recovered tracks |
|---|---|---|---|
| 1.0 | -41 (4.2 sigma) | +87 | -582 |
| 1.5 | -2 | +10 | -29 |
| 2.0 | +1 | +2 | identical |

The reference row is the flat tolerance of 0.0246 rad that the per-hit extent
replaced. It is still available with `Window::phi_per_hit` off, as
`Window::dphi_flat_rad`.

Above |eta| 0.8 the factor also compensates a phi covariance that is too small.
The ratio of the measured phi residual core to the quoted `sigma_phi`, 10 events
of the forward search, is 1.25 at |eta| < 0.8, 2.20 at 0.8-1.6, 3.40 at 1.6-2.0
and 2.66 above. The cause is the material description in the transition and
forward regions. When that is fixed, `dphi_trk_fac` should come down.

### The fetch

`find_bin_ranges()` derives the fetch from the cut, so a cut wider than the
fetch cannot be expressed. Two details:

- The cut is per hit, and the fetch runs before any hit is known, so the fetch
  uses the layer's largest hit extent (`max_hit_q_half_length()`,
  `max_hit_phi_half_extent()`).
- The phi range comes from the axis helper (`LayerOfHits::phiRangeBins()`),
  which returns a half-open range that covers the upper edge, with the mask
  applied after the increment. A hand-rolled pair of `phiBinChecked()` calls
  misses the bin holding the upper edge on every range. The phi half-width is
  clamped below pi first: a range on a circle is an arc, and a half-width of pi
  or more wraps to a small arc instead of the full circle.

Beyond the cut, the fetch adds whole bins on the bin index, so the margin
carries no float-to-bin rounding:

- `Window::phi_extra_bins` = 0. Against 1, F30 gives -1 found track and I50
  gives +1 recovered hit. One spare bin cost 20 % of build time before the line
  pre-cut (T100: 60.79 s against 48.38 s), since every fetched hit costs a
  plane solve per candidate.
- `Window::q_extra_bins` = 1. The per-hit q cut is surface-referenced and the
  fetch uses the raw `dq_track`, so at high |eta| the cut reaches past the
  fetch. With 0, I50 loses 1373 of 43358 fully recovered tracks, all in
  pixel-barrel hits of tracks that also cross the discs. F30 is unaffected. A
  fetch widened by the surface reference would remove the need for the spare
  bin. With the line pre-cut in, a spare bin costs little, since most fetched
  hits are rejected before the plane solve.

## Hit extents

**Code:** `LayerOfHits::HitInfo`, `hit_phi_half_extent_of()` and
`hit_r_half_extent_of()` in `HitStructures.cc`.

`HitInfo` caches, per hit, the extents the window and the pre-cut use. All are
`hl_fac` times a sigma from the hit's own covariance, with `hl_fac` = 3 for
pixels (a Gaussian 3 sigma) and sqrt(3) for strips, where a crossing uniform
along a segment of half-length L has sigma = L/sqrt(3).

- `q_half_length`: `hl_fac * sqrt(ezz)` in the barrel. For a strip along unit
  vector u the covariance is (L^2/3) u u^T, so this is L |u_z|, the strip's
  projection onto z. Module tilt is therefore already included. In the OT 2S barrel, where
  modules are flat, the H50 median is 2.5125 cm against a nominal 2.5. In the endcap
  it is `hl_fac * sqrt(exx + eyy)`, which is sigma_r^2 + r^2 sigma_phi^2. That is
  right for radial strips and a factor sqrt(2) too large for endcap pixels,
  where sigma_x and sigma_y are similar.
- `phi_half_extent`: from sigma_phi^2 = (y^2 exx - 2xy exy + x^2 eyy) / r^4, the
  across-strip error seen from the origin. It needs no barrel/endcap branch. For
  a 2S strip it is expected to be `hl_fac * pitch / sqrt(12)` / r, i.e.
  45 um / r for a 90 um pitch.
- `qbar_half_extent`: `hl_fac * sigma_r` in the barrel, 0 in the endcap. It is
  used by the line pre-cut for tilted strips.

## Line pre-cut

**Code:** `MkFinderV2p2::select_hits()`, the per-candidate `pc_*` block and the
per-hit test before the prefetch. Switches and slack factors in
`Config::V2p2::PreCut`, per-hit input `LayerOfHits::hit_qbar_half_extent()`.

**What it does.** Every hit fetched from the binnor used to go straight into the
pre-selection batch, where it costs a Hermite solve onto its module plane and a
lookup of the module geometry. Only then is the real dq/dphi cut applied. The
pre-cut rejects most of those hits earlier, using only the cached `HitInfo`
array and four numbers per candidate.

Between its two layer crossings `m_sp1` and `m_sp2` the track is taken as a
straight line in (qbar, q) and in (qbar, phi). Here qbar is r in the barrel and
z in the endcap. The line is evaluated at the hit's own qbar, and the hit is
dropped when its q or phi is further from the line than a tolerance that is
looser than the real cut.

**Tolerances**, with g the slope of the line:

- q: `dq_slack * dq_trk_fac * dq_track * (1 + g^2) + dq_hit_fac * hit_q_half_length + qbar_fac * |g| * hit_qbar_half_extent`
- phi: `dphi_slack * dphi_trk_fac * dphi_track + hit term of the real cut + qbar_fac * |g_phi| * hit_qbar_half_extent`

The tolerances are written in terms of the real cut's own factors, so the
pre-cut stays looser when those factors change.

**Why each term.**

- `(1 + g^2)` equals 1/sin^2(theta) for a track from the origin. It references
  the track's q error to the layer surface. The real cut does the same per hit
  (`surface_referenced_dq()`). Without it the pre-cut would drop hits that the
  real cut accepts at high |eta|.
- `hit_qbar_half_extent` is `hl_fac * sigma_r` of the hit, in the barrel. The
  qbar of a hit is its centroid radius. For a strip on a tilted module the true
  crossing lies somewhere along the strip, and the strip direction has a radial
  component, so the centroid radius is uncertain by up to L sin(tilt). Through
  the slope of the line that becomes a shift in q and in phi. In the endcap
  qbar is z, and the same term is not needed, so it is skipped there on a
  per-layer branch.

**Offline estimate that justified it.** `val_qprecut_*` in
`standalone/RdfTrace/ValProp.cc` replays the pre-cut over the trace. It counts,
for each scanned candidate-hit pair, whether the pre-cut would reject it and
whether the real cut accepts it. A false reject is a pair the real cut accepts
and the pre-cut drops. Production window, trace build:

| | forward: pairs | rejected | false rejects | inward: pairs | rejected | false rejects |
|---|---|---|---|---|---|---|
| q, slack 2 (as built) | 11.9 M | 79.3 % | 82 | 1.83 M | 85.8 % | 0 |
| q slack 2, and phi slack 1.5 (as built) | 11.9 M | 89.8 % | 405 | 1.83 M | 89.1 % | 15 |
| q, slack 2, without the surface factor | 11.9 M | 80.7 % | 17610 | 1.83 M | 89.0 % | 10381 |
| q, slack 1 | 11.9 M | 82.2 % | 2278 | 1.83 M | 89.7 % | 554 |
| q, slack 1, without the qbar term | 11.9 M | 84.3 % | 87659 | 1.83 M | 90.1 % | 3364 |

Forward is 5 events, inward 10 events. The real cut accepts
8.2 % and 7.6 % of the pairs. All
the q false rejects without the qbar term are in the tilted OT PS barrel layers.
The phi false rejects sit mostly in the OT 2S barrel.

**Physics, with the pre-cut in the finder.** Paired against the same
configuration with the pre-cut off, production window: F30 found tracks, fakes,
duplicates and the momentum-resolution width are identical to the unit, and so
are the I50 counts and truth-matched hit efficiencies, with q alone and with q
and phi.

**Time.** T100: 47.95 s without the pre-cut, 35.36 s with q, 33.81 s with q and
phi. V2 (`selectHitIndicesV2`) takes 62.95 s on the same events.

In one forward event the pre-cut rejects 2.87 M of 3.18 M scanned pairs
(90.4 %). The 306 k pairs left go through the plane solve, and 79 % of them
pass the real cut.

## Reduction and hit ordering

**Code:** `CandIncoming::m_pqueue`, `preselect_hit_batch()`,
`order_hits_along_path()`.

**Reduction.** Hits that pass the cut go into a bounded priority queue keyed on
ddphi, one queue per sub-layer, capped at `V2p2Params::max_presel_hits`. One queue
per sub-layer gives each sensor of a pair its own budget, so a busy sensor
cannot crowd out the other. The cap sits upstream of the in-layer search and
bounds what it can see. The queues and `m_layer_hits` use inline storage
(`InlineVec`) sized by `V2p2Params::max_presel_hits_limit`, because a
`CandIncoming` lives for one layer and a `std::vector` there allocated and freed
on every layer. The heaps run `std::push_heap` / `std::pop_heap` exactly as
`std::priority_queue` does. Overlap availability varies strongly across the
detector, so the cap should become per layer.

**Ordering.** The queues drain into one list per candidate, `m_layer_hits`,
sorted by `dir * dalpha`. dalpha is monotone in path length for all hits of one
candidate, since they are reached from the same state with the same curvature.
It is a signed turning angle, so its sign goes with the charge; `dir` is the
sign of its change from `m_sp1` (the near bounding surface in the search
direction) to `m_sp2` (the far one), which makes the key path order for both
charges and either search direction. Until 2026-10-01 `dir` was the search
direction alone, so for one charge every layer's hits were walked in reverse
path order. ddphi decides which hits survive; path order decides the order the
search walks them. The queues are copied out and sorted once per candidate. Two
hits on one module share the crossing with its plane and so have equal keys;
they go larger ddphi first, then by hit index, so the order does not depend on
the heap layout. A Kalman update moves the trajectory by about the hit
resolution, so it can only swap hits that are already that close.

The sub-layers of an OT PS barrel pair are radially interleaved (the tilted modules make
the partner sensor sit at larger or smaller r module by module), so path order
cannot be assumed from the layer index. Per-hit precision is not an ordering
key: `q_half_length` already distinguishes a P macro-pixel from an S strip,
and it enters the score.

The trace records `sub_rank` (rank by ddphi within the sub-layer) and
`full_rank` (across both), for measuring pre-selection against truth.

## In-layer combinatorial search

**Code:** `expand_in_layer()`, `harvest_nodes()`, `CandNode`,
`MkFinderV2p2::m_node_arena`.

**What it replaced.** The first version of v2p2 ran one Kalman update per
pre-selected hit, all from the candidate's incoming state, and kept the lowest
chi2. That takes at most one hit per layer. The best-hit path was removed once
the search replaced it. On H50 a track crossing an OT layer pair, barrel or
endcap, leaves both sensors' hits in 91-94 % of crossings, and 26-28 % of
crossings carry three or more hits, because an overlapping module is itself a
stack.

**What it does.** A path through the layer is an increasing sequence of
positions in `m_layer_hits`. Per candidate, the search is

```
frontier = { kalman(state, h) : h in m_layer_hits }               # depth 0
for depth = 1 .. max_hits_per_layer - 1, while frontier is not empty:
  next = { kalman(n.state, h) : n in frontier,
           h after n's hit in m_layer_hits, h not on a module of n's path }
  frontier = next
```

where `kalman()` keeps an outcome only if its acceptance chi2 passes
`hit_chi2_cut`. Every node of every depth enters the end-of-layer selection.
The code computes the same nodes in a different order, to keep the Matriplex
lanes full (see "Running Kalman batches" in "Layer processing").

Depth 0 starts from the candidate's state, and the Hermite has already solved
each crossing, so the propagation is handed the path length. Deeper nodes start
from an updated state, for which no crossing is known, and propagate-to-plane
solves it. The one-point Hermite (`Hermite3D::calculate_coeffs(sp, inv_k,
dalpha)`) would be cheaper there and is not yet used.

Forward-only traversal of an ordered list reaches every subset of the hits
exactly once. Taking both hits of an overlap therefore needs no special case,
and there is nothing to de-duplicate. Path order also puts the Kalman updates
in the order the track meets the hits.

A second hit from the same module as an earlier hit of the path is not taken.
That is a split cluster or another track's hit, and a split cluster is one
measurement: taking it twice shrinks the covariance without adding information.

`V2p2Params::max_hits_per_layer` (4) caps the hits per path. Forward, on 10
events, found tracks stop changing at depth 3. Inward on 20 events, fully
recovered tracks go from 82.4 % at depth 2 to 90.8 % at depth 4.

**Storage.** Nodes live in one `std::vector<CandNode>` per finder, reached by
index, with a single parent index each. A node exists to be walked back from,
so a surviving leaf can register its hits into the CombCandidate; nothing needs
to walk forward. `end_layer()` clears the vector and keeps its capacity. Slots
of a selected CombCandidate go to a free list, `m_node_free`, reused last in,
first out. A node with live children must not be reused, so slots are only ever
freed a whole CombCandidate at a time, after it has materialised.

The arena is per finder and not per CombCandidate because one Kalman batch
draws lanes from several candidates and CombCandidates. Breadth-first by depth
means a parent is always in the arena before its children are created. Nothing
depends on parents having lower indices: each depth's harvest is handed on as a
list of indices.

A `CandNode` carries a full `TrackState` and lives for one layer. A
`HoTNode` in `CombCandidate::m_hots` is 12 bytes and lives for the event. That
difference is why the two are separate.

**Effect.** F30, cap 3, against the best-hit path before its removal: found
tracks 18455 to 18846, fakes 3217 to 2765. Truth-matched hits per found track go
from 10.78 / 11.55 / 11.32 to 12.70 / 14.32 / 14.29 in barrel / transition /
endcap. Inward on 20
events, at depth 2: pixel hits recovered 75.7 % to 92.2 %,
fully recovered tracks 52.8 % to 82.4 %.

`maxCandsPerSeed` 3 is the production value with the search on. From 3 to 6
F30 gains 39 found tracks. Before the line pre-cut, T100 build time with the
search on was 60.79 s at cap 3 and 76.45 s at cap 4.

## Kalman hit acceptance

**Code:** `plane::chi2_scaled()` in `KalmanOpArgs::do_kalman_stuff()`,
`KalmanOpArgs::tsChi2Acc`, `harvest_nodes()`,
`V2p2Params::chi2_trk_fac` and `hit_chi2_cut`.

**The cut.** A hit is accepted if

```
chi2_acc = r^T (chi2_trk_fac^2 * C_trk + C_hit)^-1 r  <  hit_chi2_cut
```

with r the residual and C_trk, C_hit the track and hit covariances, all in the
module's local 2-D frame. Defaults: `chi2_trk_fac` 3, `hit_chi2_cut` 9.21, the
99 % point of chi2 with 2 degrees of freedom. The Kalman update and the
layer-step score use the plain chi2 with the unscaled covariance; only the
accept/reject decision uses `chi2_acc`.

**Why two factors.** The cut used to be the plain chi2 below a constant 30.
That one number stood for two things: a probability, and compensation for a
track covariance that is known to be too small, by 1.2-2x in sigma on the
inward search and more above |eta| 1.6. The hit covariance comes from the
cluster and matches the sensor pitch, so it is taken as it is. The factor on
the track covariance carries the deficit, and the cut is left as a probability.
The q window is split the same way, see "Search window".

**Measured.** F30, paired against the plain chi2 < 30:

| cut | chi2_trk_fac | d found | d fakes |
|---|---|---|---|
| 20 | 1 | -129 | +322 |
| 13.8 | 1 | -382 | +781 |
| 9.21 | 1 | -953 | +1612 |
| 30 | 1.5 | +110 | -225 |
| 13.8 | 2 | +95 | -187 |
| 9.21 | 2 | +19 | -44 |
| 20 | 2 | +115 | -233 |
| 13.8 | 3 | +118 | -242 |
| **9.21** | **3** | **+115** (8.3 sigma) | **-231** (12.3 sigma) |

Lowering the plain chi2 cut loses tracks and adds fakes at every value, so the
30 was compensating the covariance, not merely generous. With the track factor
the gain saturates at 2-3. It sits at |eta| 1.5-2.5, where the phi pull is
furthest from 1. The common-subset d(pT)/pT width does not move (0.02078 to
0.02085).

**These counts reward track length.** They use quality-val's association, in
which a correctly matched track with fewer than `nMinFoundHits` (10) hits is not
found and counts as a fake, as does a track on an unseeded pileup particle. With
MTV's association (a fake is a track matched to no sim track, no minimum length),
20 events, 9.21 / 3 against chi2 < 30 is -32 found and -4 fakes: no measurable
change. The acceptance cut is kept as the first-order replacement of the fixed
30. What removed v2p2's extra fakes is the hit bonus, see "Layer-step score".

I50: pixel hits recovered 198340 to 200841 (93.8 % to 95.0 %), fully
recovered tracks 39167 to 39798 (90.3 % to 91.8 %). True pixel-barrel hits on
barrel-only tracks go from 96.6 % to 97.3 %, which is the one-hit-per-layer
ceiling.

**Time.** T100 on phi3 (Skylake-SP, `-mavx`), forward search, with F30 against
the plain chi2 < 30:

| cut | chi2_trk_fac | d found | d fakes | s | |
|---|---|---|---|---|---|
| 30 | 1 | | | 22.62 | |
| 30 | 1.0001 | 0 | 0 | 22.66 | second chi2 only |
| 13.8 | 1.5 | -5 | +19 | 22.60 | same physics |
| 9.21 | 2 | +19 | -44 | 22.70 | +0.3 % |
| 9.21 | 2.5 | +94 | -187 | 22.97 | +1.6 % |
| 9.21 | 3 | +115 | -231 | 23.14 | +2.3 % |

The second 2x2 chi2 costs 0.2 %. At the old physics the time is the old time:
every pre-selected hit goes through the Kalman update either way, and the cut
only decides which results become nodes. Time grows with the true hits kept.
What reaches the Kalman is set upstream, by the window and the line pre-cut.

**Tried and not taken.**
- f a function of |eta| (a ramp from f_lo below 0.8 to 3 above 1.6): central
  f_lo 1 / 1.5 / 2 costs -245 / -76 / -24 F30 found and -2240 / -538 / -164 I50
  fully recovered tracks. The central region needs the full factor at the chi2
  level, although the phi pull there is 1.25.
- A separate f for the second and later hits within a layer: f 1.5 is a null
  forward and -53 on I50, and it saves no time, since the cost is set by what
  the first hit in the layer accepts.
- The gain is flat in pT above 0.9 GeV, so it is not multiple scattering.
- The hit covariance short instead of the track's
  (h^2 C_hit, as f = 1/h at cut h^2 * 9.21): h 2 costs -776 F30 found and -6601
  I50 fully recovered tracks. The deficit is on the track side.
- f 4 and 5: within noise of 3 in both directions, so one value serves forward
  and inward.

The factor 3 is a compensation and should come down when the track covariance
is fixed (the forward material model), as `Window::dphi_trk_fac` should.

## End-of-layer selection

**Code:** `select_and_materialise()`, `offer_best_short()`, `SelEntry`.

Everything that could continue a CombCandidate competes in one list on one
score:

- every in-layer path of every `CandIncoming`;
- every `CandIncoming` as the hole it would record, including those that did
  find paths;
- every TrackCand that did not become a `CandIncoming` (stopped, skipped at
  pull-in), at its unchanged score.

The list is partially sorted and the top `maxCandsPerSeed` survive. One sort is
possible because the score is additive over layer steps
([Layer-step score](#layer-step-score)): paths from different candidates are on
one scale. It runs at the end of the layer because one CombCandidate's
candidates can be spread over several batches.

**The hole competes with the hits.** A candidate may prefer a well-fitting wrong
hit to a hole, and no score of a single hit can tell. Offering the hole as a
competitor defers the decision until the paths can be compared. A candidate
outside the layer (WSR) is not a hole and competes unchanged.

**Reserve a hole slot** (`InLayer::reserve_hole_slot`, off by default). Taking a
hit always shrinks the covariance. When every survivor took a hit, no branch
remains that allows an earlier hit to have been wrong. This option keeps the
best decliner in the last slot. It is a beam policy, not a score term.

**Best short** (`InLayer::best_short`, off by default, outward only). A stopped
candidate cannot be extended, so it leaves the beam and the best one is kept on
the CombCandidate, the way V1's `CandCloner` does; `mergeCandsAndBestShortOne`
re-inserts it at the end if it still wins. Inward the trailing end of the track
is its head, and a truncated candidate is not a shorter track.

Survivors are built as copies before the CombCandidate is touched, because
several survivors can descend from one TrackCand.

`TrackCand::score_` holds the layer-step score during finding. At the end of the
search it is overwritten by `track_score_func`, which is what output and
cross-seed duplicate removal use. The two scores are never compared with each
other.

## Layer-step score

**Code:** `V2p2Score.h`, `LayerStepFeatures`, `Config::V2p2::Score`.

This score is separate from `track_score_func`, which stays as it is for V1/V2.
v2p2's final pick of a candidate per seed uses a track scorer built on it, see
"The final pick" below. It differs from `track_score_func` in three ways.

1. **Additive over layer steps.** A candidate's score is the sum of one term per
   layer. That is what makes the end-of-layer selection a single flat sort.
   `track_score_func` is additive only for some parameter choices: with the
   default scorer the hit bonus is linear in the hit count, so the total is
   quadratic, and the two hole penalties differ, so reclassifying a hole from
   tail to inside changes the total.
2. **No tail holes.** `TrackCand` reclassifies holes as inside or tail when a
   later hit is found, because the global formula is evaluated once at the end.
   Scoring each step when it happens has no such question. What the tail
   penalty approximated depends on direction: outward, trailing holes are at
   large radius where a track may leave the detector; inward, they are at small
   radius, where the track must have come from. The score has separate outward
   and inward parameter sets (`V2p2Params::score_outward`, `score_inward`) for
   this; today they differ in the hit efficiency (outward 0.10, with 0.14 in the
   OT PS barrel and 0.20 in the OT 2S barrel; inward 0.90).
3. **A feature struct.** `LayerStepFeatures` carries what the scores read of a
   layer step: hits taken, the chi2 sum, ln det V of the residual covariance,
   local hit density, hole kind, the layer and direction. A new term adds a
   member without changing a function signature: a fact of the step or of the
   candidate is filled in `fill_step_geometry()`, a per-hit one accumulated in
   `harvest_nodes()`. The struct can be written to the trace for fitting a
   score offline.

A term keyed on the step should take its origin from the layer of the
candidate's last found hit (`getLastFoundHitLyr()`), not from the plan
neighbour: a candidate that missed a layer steps across two.

**Two forms.** `Score::mode` 0 is linear: a bonus per hit and per overlap hit,
minus a chi2 weight, minus a penalty per hole kind. Mode 1 is a log-likelihood
ratio against "this layer produced no hit":

```
take hit j :  ln(eps/(1-eps)) - ln(2 pi) - chi2_j/2 - ln(det V_j)/2 - ln rho
take none  :  0
```

with eps the per-layer hit efficiency (`hit_eff`), V the 2x2 residual
covariance and rho the local hit density. The hit-versus-hole break-even is then
derived rather than tuned. It assumes chi2 is trustworthy, which it is to the
extent the covariance is.

**Terms fixed at their mean** (`Score::use_rho`, `use_detv`). rho and eps enter the
likelihood the same way, so replacing ln rho by a constant equal to its measured
mean removes only its variation. On F30 in mode 1 with eps 0.99, removing the
variation of rho costs 114 found tracks and removing that of det V costs 67.
The default is mode 1 with `hit_eff` 0.90, see "The default" below.

**The hit bonus.** In mode 0 a hit beats a hole when
`hit_bonus - chi2_weight * chi2 > -miss_penalty`, i.e. below chi2 =
`(hit_bonus + miss_penalty) / chi2_weight`. With the former bonus of 30 that was
chi2 38. Almost every hit that passed the acceptance cut then beat a hole. A
low-pT candidate that ran into a neighbouring particle's outer-tracker hits took
them. In CMSSW MTV on ttbar PU200 D121, v2p2 had a fake rate of 4.67 % against
1.98 % for V2. The extra fakes were tracks of 6-16 hits with a correct pixel
seed, and their outer-tracker hits belonged mostly to another particle. The
bonus is 3, so a hole wins above chi2 11.

MTV, `TrackBuilding/initialStep`, 100 events, V2 = CMSSW_20_1_0_pre2 as
released, v2p2 = this branch without the rin fix:

| | efficiency | fake rate | mean hits |
|---|---|---|---|
| V2 | 66.92 % | 1.98 % | 14.17 |
| v2p2, bonus 30 | 66.57 % | 4.67 % | 15.57 |
| v2p2, bonus 10 | 66.89 % | 3.27 % | 15.32 |
| v2p2, bonus 5 | 66.94 % | 2.79 % | 15.05 |
| v2p2, bonus 3 | 66.94 % | 2.57 % | 14.72 |

The remaining excess is at 1.7 < |eta| < 2.7 and pT < 0.5 GeV. Standalone, 50
events, MTV association: below a bonus of 3 the found and fake counts stop
improving, while tracks become shorter than V2's, the d(pT)/pT width grows and
duplicates rise.

**The default.** The linear score cannot remove the forward fakes without
shortening the tracks. The log-likelihood can: a hit has to be better than the
local hit density predicts for a random one, so a wrong hit in a dense region
loses to a hole. With `hit_eff` 0.95 it added fakes at 0.9 < |eta| < 1.7; 0.90
does not. Same MTV setup; the mode 1 row also has the in-layer path order fix
and v2p2's own backward fit, the mode 0 row is before both:

| | efficiency | fake rate | duplicates | mean hits |
|---|---|---|---|---|
| V2 | 66.92 % | 1.98 % | 0.18 % | 14.17 |
| v2p2, mode 0, bonus 3 | 66.94 % | 2.55 % | 0.12 % | 14.79 |
| v2p2, mode 1, hit_eff 0.90 | 67.03 % | 2.08 % | 0.15 % | 14.71 |

What is left of the fake excess is in the barrel (1.15 % against 0.93 %) and at
0.9 < |eta| < 1.7 (2.87 % against 2.36 %).

**The working point since.** A joint tuning (condor scans and CMSSW MTV, 2026-10)
set the outward hit efficiency to 0.10, with 0.14 in the OT PS barrel and 0.20 in
the OT 2S barrel (`V2p2ScoreParams::hit_eff_grp`); inward it stays 0.90. It is
conditional on the seeds, the duplicate cleaner and the final fit, and is to be
redone when any of them changes.

**The final pick.** At the end of the forward search each seed keeps one
candidate. The scorer that picks it is named per iteration
(`IterationConfig::m_final_pick_track_scorer_name`, default the search scorer).
`v2p2:llh` takes the summed layer-step score of the search and subtracts
`Policy::final_pick_hole_penalty` (8) per inside hole and
`final_pick_tail_penalty` (3) per tail hole, so the pick sees the holes the
layer steps scored as 0. CMSSW MTV, `TrackBuilding/initialStep`, D121 ttbar
PU200, 100 events, CMSSW_20_1_0_pre2, on this code as of 2026-10-30:

| | efficiency | fake rate | duplicates | mean hits |
|---|---|---|---|---|
| v2p2, working point, final pick `v2p2:llh` | 67.13 % | 1.72 % | 0.13 % | 14.37 |

## V2 compatibility

**Code:** `MkFinder::selectHitIndicesV2()`.

V2 is the production hit selection in CMSSW. It keeps its own window code, the
local `Bins` struct, as in upstream CMSSW, and does not use `MkBins`, so the
v2p2 window work does not move it. Its text differs from upstream in two
places only, both following interface changes on this branch: the
`propagate_to_plane()` call takes the module position and normal instead of a
`ModuleInfo`, and three renamed helpers inside `RNT_DUMP_MkF_SelHitIdcs`.

V2 does run on this branch's mini-propagators. `propagate_to_r()` is closed
form, where upstream iterates `Config::Niter` times without a convergence test,
and it clamps the target radius when the layer cannot be reached. The two agree
where the iteration converges.

The upstream window keeps two properties that the v2p2 fetch does not: the upper
phi edge of the bin range has no "+1", so the bin holding it is not scanned,
and a margin of 2.75 half-bins covers it; and the dphi jacobian is evaluated at
the track's incoming state rather than at the layer.
