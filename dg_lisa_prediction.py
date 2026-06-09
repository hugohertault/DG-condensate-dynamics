#!/usr/bin/env python3
"""
dg_lisa_prediction.py
=====================

Falsifiable prediction of the Dark Geometry inheritance scenario:
the decollapse echo frequency and whether it lands in the LISA band.

This is the part of the framework that can be WRONG. Unlike the
numerical solver checks (which only verify that the code integrates the
chosen equation correctly), this computes a physical observable from the
conformal exponent and confronts it with the published LISA sensitivity
band. If the predicted frequency fell outside [1e-4, 1e-1] Hz, the
scenario would be excluded by LISA.

Chain of reasoning (condensate_dynamics.tex):
  1. The inheritance exponent is fixed by the conformal weight of the
     axiom, NOT fitted:
         p = (d+1)/d = 4/3      (d = 3)
  2. Lambda ~ S_parent^{-p}, with our cosmological constant
         Lambda / M_Pl^4 ~ 1e-122.
  3. Solve for the parent entropy S_parent, then the parent mass via the
     Bekenstein-Hawking relation S = 4 pi (M/M_Pl)^2  (in natural units,
     for a Schwarzschild hole S = A/4 = 4 pi G^2 M^2 / (hbar c ...)).
  4. The decollapse echo resonant frequency is the light-crossing /
     ringdown scale
         f_res = c^3 / (2 pi G M_parent).
  5. Compare to the LISA band.

NOTHING below is hard-coded to 1.7e7 Msun or 2 mHz; those are the
paper's quoted values and we check whether an independent computation
reproduces them.
"""

import numpy as np

# ---- constants (SI) ----
c    = 2.99792458e8       # m/s
G    = 6.67430e-11        # m^3 kg^-1 s^-2
hbar = 1.054571817e-34    # J s
kB   = 1.380649e-23       # J/K
Msun = 1.98892e30         # kg

# ---- framework input ----
D_SPACE = 3
P_EXPONENT = (D_SPACE + 1) / D_SPACE        # 4/3, conformal weight
LAMBDA_OVER_MPL4 = 1.0e-122                  # observed cosmological constant

# Planck units
M_Pl_kg = np.sqrt(hbar * c / G)              # Planck mass ~2.176e-8 kg

# LISA sensitivity band (Hz), standard quoted range
LISA_LO, LISA_HI = 1.0e-4, 1.0e-1
LIGO_LO, LIGO_HI = 1.0e1, 1.0e3


def parent_entropy_from_lambda(lam=LAMBDA_OVER_MPL4, p=P_EXPONENT):
    """Lambda ~ S^{-p}  =>  S ~ Lambda^{-1/p}."""
    return lam ** (-1.0 / p)


def mass_from_entropy(S):
    """Bekenstein-Hawking: S = 4 pi (M/M_Pl)^2  =>  M = M_Pl sqrt(S/(4 pi))."""
    return M_Pl_kg * np.sqrt(S / (4.0 * np.pi))


def resonant_frequency(M_kg):
    """f_res = c^3 / (2 pi G M)  (inverse light-crossing of r_s, in Hz)."""
    return c**3 / (2.0 * np.pi * G * M_kg)


def hawking_temperature(M_kg):
    return hbar * c**3 / (8.0 * np.pi * G * M_kg * kB)


def in_band(f, lo, hi):
    return lo <= f <= hi


def main():
    print("=" * 68)
    print(" Dark Geometry -- falsifiable decollapse-echo prediction")
    print("=" * 68)
    print(f" d = {D_SPACE}   =>   inheritance exponent p = (d+1)/d = {P_EXPONENT:.6f}")
    print(f" Lambda / M_Pl^4 = {LAMBDA_OVER_MPL4:.1e}")
    print()

    S = parent_entropy_from_lambda()
    M = mass_from_entropy(S)
    f = resonant_frequency(M)
    TH = hawking_temperature(M)

    print(" Derived (no fitting):")
    print(f"   parent entropy   S_parent = {S:.3e}")
    print(f"   parent mass      M_parent = {M:.3e} kg = {M/Msun:.3e} M_sun")
    print(f"   resonant freq    f_res    = {f:.3e} Hz")
    print(f"   Hawking temp     T_H      = {TH:.3e} K")
    print()

    print(" Confrontation with detectors:")
    lisa = in_band(f, LISA_LO, LISA_HI)
    ligo = in_band(f, LIGO_LO, LIGO_HI)
    print(f"   LISA band [{LISA_LO:.0e}, {LISA_HI:.0e}] Hz : "
          f"{'INSIDE  <-- testable by LISA' if lisa else 'outside'}")
    print(f"   LIGO band [{LIGO_LO:.0e}, {LIGO_HI:.0e}] Hz : "
          f"{'inside' if ligo else 'outside'}")
    print()

    # cross-check against the paper's quoted numbers (are they reproduced?)
    print(" Cross-check vs paper's quoted values:")
    M_paper, f_paper = 1.7e7, 2e-3
    print(f"   M_parent : computed {M/Msun:.2e} vs paper {M_paper:.1e} M_sun"
          f"   (ratio {M/Msun/M_paper:.2f})")
    print(f"   f_res    : computed {f:.2e} vs paper {f_paper:.1e} Hz"
          f"   (ratio {f/f_paper:.2f})")
    print()

    # sensitivity: how does the exponent choice move the prediction?
    print(" Sensitivity to the exponent (why p=4/3 is the physical one):")
    print(f"   {'reading':<22}{'p':>8}{'M_parent [Msun]':>20}{'f_res [Hz]':>14}{'band':>10}")
    for name, p in [("boundary (area)", 1.0),
                    ("conformal (axiom)", 4.0/3.0),
                    ("naive bulk", 1.5)]:
        Sp = parent_entropy_from_lambda(p=p)
        Mp = mass_from_entropy(Sp)
        fp = resonant_frequency(Mp)
        band = "LISA" if in_band(fp, LISA_LO, LISA_HI) else (
               "LIGO" if in_band(fp, LIGO_LO, LIGO_HI) else "neither")
        print(f"   {name:<22}{p:>8.4f}{Mp/Msun:>20.3e}{fp:>14.3e}{band:>10}")
    print()
    print(" The conformal exponent p=4/3 is the one selected by the bulk-")
    print(" filling dynamics (Criterion 2). This script shows it is also the")
    print(" one whose echo lands in LISA -- a falsifiable, not fitted, result.")
    print("=" * 68)


if __name__ == "__main__":
    main()
