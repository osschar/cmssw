sincos4 -- what it was, why it was removed, and how to get it back
==================================================================

Removed 2026-09-12 from the propagation code.  This file is the record; the
source no longer carries it.

WHAT IT WAS
-----------
A 4th-order Taylor pair, predating vdt:

    sin(x) ~ x - x^3/6
    cos(x) ~ 1 - x^2/2 + x^4/24

selected at run time by Config::useTrigApprox (and PropagationFlags::
use_trig_approx), and called at x = dalpha/2 -- the half turn-angle of a
propagation step.

WHERE IT WAS USED
-----------------
Five call sites:

    PropagationMPlexPlane.cc        (the only one that had already moved to vdt)
    PropagationMPlex.cc             x3
    PropagationMPlexEndcap.cc

The four R/Z sites fell back to libm std::sin/std::cos when the flag was off;
they now use vdt::fast_sincosf like the plane path already did.

Also removed: the flag itself, and the scalar helper.

THE CODE, AS IT WAS
-------------------
The scalar helper lived in the Matriplex library and is STILL THERE (see "WHAT
REMAINS" below) -- src/Matriplex/MatriplexCommon.h:

    namespace internal {
      template <typename T>
      void sincos4(const T x, T &sin, T &cos) {
        // Had this writen with explicit division by factorial.
        // The *whole* fitting test ran like 2.5% slower on MIC, sigh.

        const T x2 = x * x;
        cos = T(1.0) - T(0.5) * x2 + T(0.0416666666666666667) * x2 * x2;
        sin = x - T(0.166666666666666667) * x * x2;
      }
    }  // namespace internal

plus the Matriplex member and free function in src/Matriplex/Matriplex.h.

The switch, Config.h:

    constexpr bool useTrigApprox = true;

PropagationConfig.h carried only a note:

    // Could add: bool use_trig_approx       -- now Config::useTrigApprox = true

The five call sites, exactly as removed:

PropagationMPlexPlane.cc -- the only one already on vdt in the else branch:

    if constexpr (Config::useTrigApprox) {
      mpt::sincos4(0.5f * alpha, sinah, cosah);
    } else {
      mpt::fast_sincos(0.5f * alpha, sinah, cosah);
    }

PropagationMPlex.cc, three sites -- note the else branches call libm, not vdt:

    if constexpr (Config::useTrigApprox) {
      sincos4(ialpha * 0.5f, sinah, cosah);
    } else {
      cosah = std::cos(ialpha * 0.5f);
      sinah = std::sin(ialpha * 0.5f);
    }

    if constexpr (Config::useTrigApprox) {
      for (int n = nmin; n < nmax; ++n) {
        sincos4(id[n - nmin] * ipt[n - nmin] * kinv[n - nmin] * 0.5f, sinah[n - nmin], cosah[n - nmin]);
      }
    } else {
    #if !defined(__INTEL_COMPILER)
    #pragma omp simd
    #endif
      for (int n = nmin; n < nmax; ++n) {
        cosah[n - nmin] = std::cos(id[n - nmin] * ipt[n - nmin] * kinv[n - nmin] * 0.5f);
        sinah[n - nmin] = std::sin(id[n - nmin] * ipt[n - nmin] * kinv[n - nmin] * 0.5f);
      }
    }

    if constexpr (Config::useTrigApprox) {
    #pragma omp simd
      for (int n = nmin; n < nmax; ++n) {
        sincos4(alpha[n - nmin], sina[n - nmin], cosa[n - nmin]);
      }
    } else {
      for (int n = nmin; n < nmax; ++n) {
        cosa[n - nmin] = std::cos(alpha[n - nmin]);
        sina[n - nmin] = std::sin(alpha[n - nmin]);
      }
    }

PropagationMPlexEndcap.cc:

    if constexpr (Config::useTrigApprox) {
    #if !defined(__INTEL_COMPILER)
    #pragma omp simd
    #endif
      for (int n = 0; n < NN; ++n) {
        sincos4(alpha[n] * 0.5f, sinahTmp[n], cosahTmp[n]);
      }
    } else {
      for (int n = 0; n < NN; ++n) {
        cosahTmp[n] = std::cos(alpha[n] * 0.5f);
      }
    #if !defined(__INTEL_COMPILER)
    #pragma omp simd
    #endif
      for (int n = 0; n < NN; ++n) {
        sinahTmp[n] = std::sin(alpha[n] * 0.5f);
      }
    }

Each was replaced by a single vdt::fast_sincosf() call on the same argument
(mpt::fast_sincos on the plane path), and `#include "vdt/sincos.h"` added at
file scope in the two .cc files that did not already have it.


WHY IT WAS REMOVED -- IT WAS THE PROPAGATOR'S ENTIRE TURN-ANGLE DEPENDENCE
-------------------------------------------------------------------------
It is exact to float precision for the SMALL steps between adjacent layers,
which is what it was written for.  Beyond that the leading error x^5/120 at
x = dalpha/2, times the 2*R_c lever arm, makes the position error grow as
dalpha^5.

Median in-plane landing error of ONE step, against an exact double-precision
reference:

    |dalpha| [rad]     0.02     0.10     0.20      0.40      0.80
    sincos4 [cm]      9.1e-7   2.0e-6   9.8e-6   2.8e-4    7.7e-3
    vdt     [cm]      9.1e-7   2.0e-6   3.2e-6   4.0e-6    5.6e-6

The closed form predicts it: 2*R_c*(dalpha/2)^5/120 gives 1.5e-5 / 3.7e-4 /
9.2e-3 against the measured 9.8e-6 / 2.8e-4 / 7.7e-3 -- within 20-40 % over
three decades.

With vdt the error is FLAT in dalpha.  So the approximation was not a small
correction on top of the propagator: it WAS the propagator's whole turn-angle
dependence.

WHAT IT COST TO REMOVE: NOTHING MEASURABLE
------------------------------------------
20 events, LST-T5 sample: found 535, nH >= 80% 449, both ways; 0.88 s vs 0.90 s.

Read that carefully -- it is track COUNTS, per event and in total, over three
runs.  The states themselves were not diffed, so "identical" here means the
reconstruction outcome, not bit-level equality.

And parsFromPathL_impl already made three vdt calls against sincos4's one, so
the saving was a fraction of a fraction.

WHAT THE REAL STEP SIZES ARE, WHICH IS WHY THIS MATTERS LESS THAN IT LOOKS
--------------------------------------------------------------------------
Turn angle between consecutive hits, 56763 steps, sim tracks, 20 events:

    pT >= 0.5    median 0.0075   p90 0.164   p99 0.807   >0.4: 2.9%   >0.8: 1.0%
    pT >= 0.7    median 0.0056   p90 0.110   p99 0.379   >0.4: 0.89%  >0.8: 0.19%
    pT >= 2.0    median 0.0021   p90 0.042   p99 0.106   >0.4: 0.006% >0.8: 0

The median real step is dalpha ~ 0.006 rad, where sincos4 and vdt agree to the
digit.  The divergence above only bites in the tail -- and note large steps are
commonest in the BARREL (5.7% above 0.4 rad at |eta| < 0.8, against 1.0% at
|eta| 2.2-3), i.e. low-pT tracks curving between widely spaced barrel layers,
not the transition region.

WHAT REMAINS IN THE TREE
------------------------
The Matriplex library still carries its own implementation, now with no callers:

    src/Matriplex/MatriplexCommon.h   internal::sincos4()
    src/Matriplex/Matriplex.h         Matriplex::sincos4(), and the free function

Left alone deliberately -- that is library code, not propagation code, and
removing it would diverge from upstream Matriplex for no gain.

IF YOU WANT IT BACK
-------------------
Restore from git (the removal commit touches Config.h, PropagationConfig.h and
the three PropagationMPlex*.cc files), and gate it on step size rather than on
a global flag: it is a valid optimisation for short steps and a bad
approximation for long ones.  Only worth doing if microseconds are being
chased; on the evidence above it buys nothing at the tracking level.
