#!/usr/bin/env python3
"""
make_figure.py
==============

Generates figures/condensate_regimes.png : the three regimes of the
Hertault condensate master equation, evolved from an identical Gaussian
packet. Top row: field snapshots. Bottom row: RMS width and peak
amplitude vs time. All curves are produced by the solver in
dg_condensate_dynamics.py ; nothing is drawn by hand.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "code"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dg_condensate_dynamics import evolve

REGIMES = [
    (r"$\rho<\rho_c$  (dark energy, $m^2_{\rm eff}=+0.5$)",  +0.5, "#2c7fb8", 150.0),
    (r"$\rho=\rho_c$  (massless echo, $m^2_{\rm eff}=0$)",     0.0, "#31a354", 150.0),
    (r"$\rho>\rho_c$  (dark matter, $m^2_{\rm eff}=-0.5$)",   -0.5, "#de2d26",  40.0),
]

L, w0 = 400.0, 2.5

runs = []
for label, m2, color, tmax in REGIMES:
    r = evolve(m2=m2, gamma=0.0, L=L, w0=w0, t_max=tmax,
               absorbing=(m2 <= 0.0))
    runs.append((label, m2, color, r))

fig, axes = plt.subplots(2, 3, figsize=(13, 7))

# top row: final snapshots (and initial dashed)
x0 = np.linspace(-L/2, L/2, 4001)
phi0 = np.exp(-x0**2 / (2*w0**2))
for j, (label, m2, color, r) in enumerate(runs):
    ax = axes[0, j]
    x = r["x"]
    pf = r["phi_final"]
    # normalize tachyonic snapshot for display (it has grown enormously)
    disp = pf / np.max(np.abs(pf)) if np.max(np.abs(pf)) > 0 else pf
    ax.plot(x0, phi0, "k--", lw=1, alpha=0.5, label="initial")
    ax.plot(x, disp, color=color, lw=1.5, label="final (norm.)")
    ax.set_xlim(-L/2, L/2)
    ax.set_title(label, fontsize=9)
    ax.set_xlabel("x")
    if j == 0:
        ax.set_ylabel(r"$\phi$ (normalized)")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(alpha=0.2)

# bottom row: width and amplitude vs time
axw = axes[1, 0]
axp = axes[1, 1]
for label, m2, color, r in runs:
    t = r["t"]; w = r["width"]; p = r["peak"]
    g = np.isfinite(w)
    axw.plot(t[g], w[g], color=color, lw=1.5, label=label.split("(")[0])
    gp = np.isfinite(p) & (p > 0)
    axp.semilogy(t[gp], p[gp], color=color, lw=1.5)
axw.set_xlabel("t"); axw.set_ylabel("RMS width  w(t)")
axw.set_title("Packet width", fontsize=9); axw.grid(alpha=0.2)
axw.legend(fontsize=6, loc="upper left")
axp.set_xlabel("t"); axp.set_ylabel("peak amplitude  (log)")
axp.set_title("Peak amplitude", fontsize=9); axp.grid(alpha=0.2)

# third bottom panel: measured vs analytic tachyonic growth rate
axr = axes[1, 2]
from dg_condensate_dynamics import evolve as _ev
m2s = np.array([-0.25, -0.5, -1.0, -2.0])
meas = []
for m2 in m2s:
    rr = _ev(m2=m2, gamma=0.0, L=200.0, nx=2001, w0=60.0,
             t_max=20.0, cfl=0.5, absorbing=False, record_every=20)
    p = rr["k0"]; t = rr["t"]
    good = np.isfinite(p) & (p > 0)
    sl, _ = np.polyfit(t[good][2:-2], np.log(p[good][2:-2]), 1)
    meas.append(sl)
meas = np.array(meas)
analytic = np.sqrt(np.abs(m2s))
axr.plot(analytic, analytic, "k-", lw=1, label=r"$\lambda=\sqrt{|m^2|}$ (exact)")
axr.plot(analytic, meas, "o", color="#de2d26", ms=7, label="measured (k=0)")
axr.set_xlabel(r"$\sqrt{|m^2_{\rm eff}|}$"); axr.set_ylabel(r"growth rate $\lambda$")
axr.set_title("Tachyonic growth rate", fontsize=9); axr.grid(alpha=0.2)
axr.legend(fontsize=7, loc="upper left")

fig.suptitle("Hertault condensate: damped Klein-Gordon master equation, "
             r"$\partial_t^2\phi+\gamma\partial_t\phi=c^2\partial_x^2\phi-m^2_{\rm eff}\phi$",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(os.path.dirname(__file__), "..", "figures", "condensate_regimes.png")
fig.savefig(out, dpi=140)
fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)
