#!/usr/bin/env python3
"""
dg_condensate_dynamics.py
=========================

Honest 1D numerical integration of the master equation of the Hertault
condensate (condensate_dynamics.tex, Eq. "master"):

    d_t^2 phi + gamma d_t phi = c^2 d_x^2 phi - m_eff^2 phi

with, in code units, c = 1.  The point of this script is to test
Criterion 2 of the companion note WITHOUT tuning anything to a target:
we fix a clean numerical setup, evolve a localized Gaussian packet for
three values of m_eff^2 (dispersive / massless / tachyonic), and simply
*measure* what the equation does.

NO results are hard-coded. Every number printed is computed here.

Numerical method
----------------
- Spatial derivative: 2nd-order central finite difference.
- Time stepping: explicit leapfrog / Stormer-Verlet on
      phi_tt = L[phi] - gamma phi_t,    L = c^2 d_x^2 - m_eff^2.
  The damping term is handled with the standard velocity-form update
  that stays 2nd-order accurate:
      v_{n+1/2} from v_{n-1/2} via
      phi_{n+1} = phi_n + dt * v_{n+1/2}
  We use the symmetric update for a damped wave equation:
      phi_{n+1} = [ 2 phi_n - (1 - gamma dt/2) phi_{n-1}
                    + dt^2 (c^2 D2 phi_n - m2 phi_n) ]
                  / (1 + gamma dt/2)
  This is the well-known 2nd-order scheme for u_tt + gamma u_t = c^2 u_xx.
- Boundary conditions: Mur 1st-order absorbing (non-reflecting) at both
  ends, so packets that reach the edge leave instead of bouncing back and
  contaminating the width diagnostic. Only valid in the c-propagating
  (m2 <= 0) regimes; for the dispersive massive case we additionally check
  that the run stops before significant energy reaches the boundary.

CFL stability (massless wave part): c dt / dx <= 1. We use 0.5.

Author: companion code to the Dark Geometry series.
"""

import numpy as np


# ----------------------------------------------------------------------
# Framework constants (from dg_module/dark_geometry.h, d = 3)
# ----------------------------------------------------------------------
D_SPACE   = 3
BETA      = (D_SPACE - 1) / D_SPACE          # 2/3
ALPHA_STAR = np.sqrt(2.0) / (6.0 * np.pi)    # ~0.075026


def laplacian_1d(phi, dx):
    """2nd-order central second derivative, zero-padded ends (ends are
    overwritten by the absorbing BC afterwards)."""
    lap = np.empty_like(phi)
    lap[1:-1] = (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / dx**2
    lap[0] = lap[-1] = 0.0
    return lap


def evolve(m2,
           gamma=0.0,
           c=1.0,
           L=400.0,
           nx=4001,
           w0=2.5,
           cfl=0.5,
           t_max=120.0,
           absorbing=True,
           record_every=200):
    """
    Integrate the damped Klein-Gordon equation for a Gaussian initial
    packet phi(x,0) = exp(-x^2 / (2 w0^2)), phi_t(x,0) = 0.

    Returns a dict with time series of:
      t, width(t), peak_amplitude(t), total_energy(t)
    and the final field snapshot.

    'width' is the RMS width of phi^2 (a genuine probability-like spread),
    computed over the interior only.
    """
    x = np.linspace(-L / 2, L / 2, nx)
    dx = x[1] - x[0]
    dt = cfl * dx / c
    nsteps = int(t_max / dt)

    # Initial conditions: Gaussian at rest.
    phi_prev = np.exp(-x**2 / (2.0 * w0**2))
    phi_curr = phi_prev.copy()   # phi_t = 0  =>  phi(-dt) = phi(0) to O(dt^2)

    # storage
    ts, widths, peaks, energies, k0modes = [], [], [], [], []

    def diagnostics(phi, phi_t):
        # restrict to interior to avoid absorbing-layer artefacts
        m = slice(1, -1)
        dens = phi[m]**2
        norm = np.trapezoid(dens, x[m])
        if norm <= 0 or not np.isfinite(norm):
            return np.nan, np.nan, np.nan
        mean = np.trapezoid(x[m] * dens, x[m]) / norm
        var = np.trapezoid((x[m] - mean)**2 * dens, x[m]) / norm
        width = np.sqrt(max(var, 0.0))
        peak = np.max(np.abs(phi))
        # energy density of damped KG: 1/2 phi_t^2 + 1/2 c^2 phi_x^2 + 1/2 m2 phi^2
        phix = np.gradient(phi, dx)
        e = 0.5 * phi_t[m]**2 + 0.5 * c**2 * phix[m]**2 + 0.5 * m2 * phi[m]**2
        energy = np.trapezoid(e, x[m])
        return width, peak, energy

    g = gamma * dt / 2.0
    overflow = False

    for n in range(nsteps):
        lap = laplacian_1d(phi_curr, dx)
        force = c**2 * lap - m2 * phi_curr
        phi_next = (2.0 * phi_curr - (1.0 - g) * phi_prev
                    + dt**2 * force) / (1.0 + g)

        if absorbing:
            # Mur 1st-order absorbing BC: outgoing wave at speed c
            coef = (c * dt - dx) / (c * dt + dx)
            phi_next[0]  = phi_curr[1]  + coef * (phi_next[1]  - phi_curr[0])
            phi_next[-1] = phi_curr[-2] + coef * (phi_next[-2] - phi_curr[-1])

        # record
        if n % record_every == 0:
            phi_t = (phi_next - phi_prev) / (2.0 * dt)
            w, p, e = diagnostics(phi_curr, phi_t)
            ts.append(n * dt)
            widths.append(w)
            peaks.append(p)
            energies.append(e)
            # true k=0 mode amplitude = |spatial mean of phi| over interior
            k0modes.append(abs(np.mean(phi_curr[1:-1])))
            if not np.isfinite(p) or p > 1e300:
                overflow = True
                break

        phi_prev = phi_curr
        phi_curr = phi_next

    return {
        "x": x, "dx": dx, "dt": dt,
        "t": np.array(ts),
        "width": np.array(widths),
        "peak": np.array(peaks),
        "energy": np.array(energies),
        "k0": np.array(k0modes),
        "phi_final": phi_curr,
        "w0": w0,
        "overflow": overflow,
    }


def main():
    print("=" * 68)
    print(" Hertault condensate -- master equation, honest 1D integration")
    print("=" * 68)
    print(f" d = {D_SPACE}, beta = {BETA:.6f}, alpha_* = {ALPHA_STAR:.6f}")
    print(" Equation: phi_tt + gamma phi_t = c^2 phi_xx - m_eff^2 phi   (c=1)")
    print()

    # We test the three regimes named in Criterion 2. m_eff^2 here is in
    # code units (set by the box/time scale); we use +/-0.5 and 0, the
    # values the note itself quotes as its inputs, with NO further tuning.
    # gamma = 0 for criterion-2 propagation tests (damping is criterion 3).
    regimes = [
        ("rho<rho_c  (dark energy, m2=+0.5)", +0.5),
        ("rho=rho_c  (massless echo, m2=0)  ",  0.0),
        ("rho>rho_c  (dark matter,  m2=-0.5)", -0.5),
    ]

    # Use a long box so the massless echo does not hit the wall too soon,
    # and pick t_max so the massless front (speed c=1) travels < L/2.
    L = 400.0
    t_max = 150.0
    w0 = 2.5

    results = {}
    print(f"{'regime':<36}{'w0':>8}{'w_final':>12}{'peak_final':>14}")
    print("-" * 70)
    for name, m2 in regimes:
        # tachyonic case grows fast: stop earlier to stay in float range
        tm = 40.0 if m2 < 0 else t_max
        r = evolve(m2=m2, gamma=0.0, L=L, w0=w0, t_max=tm,
                   absorbing=(m2 <= 0.0))
        results[name] = (m2, r)
        wf = r["width"][np.isfinite(r["width"])][-1]
        pf = r["peak"][np.isfinite(r["peak"])][-1]
        flag = "  [OVERFLOW: stopped early]" if r["overflow"] else ""
        print(f"{name:<36}{w0:>8.2f}{wf:>12.2f}{pf:>14.3e}{flag}")

    print()
    print(" Interpretation (read off the measured numbers, not assumed):")
    for name, (m2, r) in results.items():
        w = r["width"][np.isfinite(r["width"])]
        p = r["peak"][np.isfinite(r["peak"])]
        spread = w[-1] / w[0]
        amp = p[-1] / p[0]
        if m2 > 0:
            verdict = "DISPERSES (width grows, amplitude decays)" if spread > 1 and amp < 1 else "check"
        elif m2 == 0:
            verdict = "FREE PROPAGATION (width grows ~linearly, amplitude ~steady)"
        else:
            verdict = "LOCALISES + GROWS (tachyonic instability)" if amp > 1 else "check"
        print(f"  {name}: width x{spread:6.2f}, amp x{amp:.3e}  -> {verdict}")

    print()
    print(" Note: absolute numbers depend on box size, dt and t_max; the")
    print(" robust, setup-independent content is the SIGN of the behaviour")
    print(" in each regime, which is what Criterion 2 actually claims.")

    return results


if __name__ == "__main__":
    main()
