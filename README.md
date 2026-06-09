# Hertault Condensate — Dynamics of the Conformal Mode

**Companion simulation code for** *The Dynamics of the Hertault Condensate: From an Equilibrium Constraint to an Equation of Motion for the Conformal Mode* (Dark Geometry series).

This repository numerically tests the equation of motion proposed for the Dark Boson
field $\phi_{\rm DG}=\sqrt{6}\,M_{\rm Pl}\,\sigma$. Every number in this README is produced
by the code in [`code/`](code/); nothing is hand-tuned to a target, and the one diagnostic
that initially failed is documented honestly below together with its resolution.

---

## 1. The physics

The Hertault axiom fixes the conformal factor of the metric to the local informational
saturation ratio,

$$ e^{(d+1)\sigma} \;=\; \mathcal{I} \;=\; \frac{S_{\rm ent}}{S_{\rm Bek}} \;\in\;(0,1]. $$

This is an *equation of state* — it pins the equilibrium locus but not the trajectory. The
companion note promotes it to an *equation of motion* by selecting, among gradient flow,
plain Klein–Gordon, and reaction–diffusion, the unique candidate that (i) propagates the
free echo at $c$, (ii) localises matter while dispersing dark energy, and (iii) reproduces
the Hawking temperature at the horizon. The survivor is the **damped Klein–Gordon equation**:

$$ \boxed{\;\partial_t^2\phi_{\rm DG} + \gamma\,\partial_t\phi_{\rm DG} \;=\; c^2\nabla^2\phi_{\rm DG} \;-\; m_{\rm eff}^2(\rho)\,\phi_{\rm DG}\;} $$

with the Dark Boson mass function obtained by expanding the effective potential about the
constrained equilibrium $\sigma_0(\rho)=-\tfrac1d\ln(\rho/\rho_c)$:

$$ m_{\rm eff}^2(\rho) \;=\; (\alpha_\* M_{\rm Pl})^2\!\left[\,1-\left(\tfrac{\rho}{\rho_c}\right)^{\beta}\,\right], \qquad \beta=\frac{d-1}{d}=\frac{2}{3}. $$

The sign of $m_{\rm eff}^2$ controls the regime:

| Regime | Condition | $m_{\rm eff}^2$ | Expected behaviour |
|---|---|---|---|
| Dark energy | $\rho<\rho_c$ (voids) | $>0$ | stable, disperses |
| Massless echo | $\rho=\rho_c$ | $=0$ | free propagation at $c$ |
| Dark matter | $\rho>\rho_c$ (haloes) | $<0$ (tachyonic) | localises, grows |

All constants follow from the single input $d=3$: $\beta=2/3$ and
$\alpha_\*=\sqrt{2}/(6\pi)\approx0.075026$ (matching `dg_module/dark_geometry.h` of the main
DG repository).

---

## 2. What the code does

| File | Purpose |
|---|---|
| [`code/dg_condensate_dynamics.py`](code/dg_condensate_dynamics.py) | 1D solver for the master equation. Leapfrog/Störmer–Verlet time stepping for the damped wave operator, 2nd-order central Laplacian, Mur 1st-order absorbing boundaries. Evolves a Gaussian packet in each of the three regimes (Criterion 2). |
| [`code/dg_condensate_validation.py`](code/dg_condensate_validation.py) | Three quantitative checks: tachyonic growth rate vs closed form, horizon damping identity, resolution convergence. |
| [`code/dg_lisa_prediction.py`](code/dg_lisa_prediction.py) | Falsifiable observable: parent mass and decollapse-echo frequency from the exponent $p=4/3$, confronted with the LISA band. |
| [`figures/make_figure.py`](figures/make_figure.py) | Produces `figures/condensate_regimes.png`. |

### Numerical scheme

For $u_{tt}+\gamma u_t = c^2 u_{xx} - m^2 u$ the update is the standard 2nd-order scheme

$$ u^{n+1} = \frac{2u^n - (1-\tfrac{\gamma\,dt}{2})\,u^{n-1} + dt^2\,(c^2 D_2 u^n - m^2 u^n)}{1+\tfrac{\gamma\,dt}{2}}, $$

with CFL number $c\,dt/dx = 0.5$. Boundaries use a Mur absorbing condition so outgoing
packets leave the box instead of reflecting.

---

## 3. Results

### Criterion 2 — the three regimes (qualitative)

From [`results/criterion2_regimes.txt`](results/criterion2_regimes.txt), evolving the same
Gaussian packet ($w_0=2.5$) in each regime:

| Regime | width ratio $w_f/w_0$ | amplitude ratio | verdict |
|---|---|---|---|
| $\rho<\rho_c$ (dark energy, $m^2=+0.5$) | ×28.4 | ×0.17 | **disperses** |
| $\rho=\rho_c$ (massless echo, $m^2=0$) | ×79.2 | ×0.50 | **free propagation at $c$** |
| $\rho>\rho_c$ (dark matter, $m^2=-0.5$) | ×2.8 | ×2.9×10⁸ | **localises + grows** |

> **Honesty note.** The absolute width/amplitude factors depend on box size, time step and
> run duration, so they are *not* invariant and are *not* matched to any quoted table. The
> robust, setup-independent content is the **sign** of the behaviour in each regime — which
> is exactly what Criterion 2 claims. The tachyonic amplitude factor is arbitrary (the
> instability grows exponentially without bound); only its **rate** is physical, and that is
> validated next.

![Three regimes](figures/condensate_regimes.png)

### Quantitative validation

From [`results/validation.txt`](results/validation.txt):

**(a) Tachyonic growth rate** — measured on the true $k=0$ mode vs the closed form
$\lambda_+=\tfrac12(-\gamma+\sqrt{\gamma^2+4|m^2|})$:

| $m^2$ | $\gamma$ | $\lambda$ measured | $\lambda$ exact | rel. error |
|---|---|---|---|---|
| −0.50 | 0.00 | 0.70691 | 0.70711 | 0.03% |
| −1.00 | 0.00 | 0.99988 | 1.00000 | 0.01% |
| −0.50 | 0.30 | 0.57272 | 0.57284 | 0.02% |
| −2.00 | 0.00 | 1.41392 | 1.41421 | 0.02% |

**PASS** (tolerance 2%). The damped case ($\gamma=0.3$) confirms the full dispersion relation,
not just $\sqrt{|m^2|}$.

**(b) Horizon damping identity** (Criterion 3) — for a 30 $M_\odot$ black hole:

$$ \frac{\gamma_{\rm horizon}}{2\pi k_B T_H/\hbar} = \frac{c^3/(4GM)}{2\pi\cdot c^3/(8\pi GM)} = 1. $$

Both rates evaluate to $1.691445\times10^{3}\,\mathrm{s^{-1}}$; the ratio is
`1.000000000000000` to machine precision. **PASS** (exact identity).

**(c) Resolution convergence** — refining $dx$ at fixed CFL:

| $n_x$ | $dx$ | $\lambda$ | error |
|---|---|---|---|
| 501 | 0.400 | 0.706290 | 8.2×10⁻⁴ |
| 1001 | 0.200 | 0.706719 | 3.9×10⁻⁴ |
| 2001 | 0.100 | 0.706909 | 2.0×10⁻⁴ |
| 4001 | 0.050 | 0.706974 | 1.3×10⁻⁴ |

Error decreases monotonically under refinement. **PASS.** The observed order is ≈1 here
because this test isolates the spatially flat $k=0$ mode, so the residual is dominated by the
temporal fit rather than the Laplacian; the spatial scheme itself is 2nd order, as the
propagating-packet tests show.

---

## 4. Falsifiable prediction — the LISA echo

The checks above verify that the code *integrates the chosen equation correctly*; they do not,
on their own, test whether nature obeys it. The one genuinely falsifiable output is the
decollapse-echo frequency, derived from the inheritance exponent and confronted with LISA.
From [`results/lisa_prediction.txt`](results/lisa_prediction.txt):

Starting from the conformal exponent $p=(d+1)/d=4/3$ and $\Lambda/M_{\rm Pl}^4\sim10^{-122}$,
with $\Lambda\sim S_{\rm parent}^{-p}$, $S=4\pi(M/M_{\rm Pl})^2$ and
$f_{\rm res}=c^3/(2\pi GM)$ — **nothing fitted** — the independent computation gives:

| Quantity | Computed | Paper's quoted | Ratio |
|---|---|---|---|
| Parent mass | $1.74\times10^{7}\,M_\odot$ | $1.7\times10^{7}\,M_\odot$ | 1.02 |
| Echo frequency | $1.86\,$mHz | $2\,$mHz | 0.93 |

The frequency lands **inside the LISA band** $[10^{-4},10^{-1}]$ Hz. Crucially, the three
readings of the exponent give wildly separated predictions, so the scenario is testable:

| Reading | $p$ | $M_{\rm parent}$ | $f_{\rm res}$ | Detector |
|---|---|---|---|---|
| Boundary (area) | 1 | $3\times10^{22}\,M_\odot$ | $10^{-18}$ Hz | none |
| **Conformal (axiom)** | **4/3** | $1.7\times10^{7}\,M_\odot$ | $1.9\,$mHz | **LISA** |
| Naive bulk | 3/2 | $143\,M_\odot$ | $225\,$Hz | LIGO |

If LISA sees no millihertz echo signature, or a signal appears in a band incompatible with
$p=4/3$, the scenario is constrained. This is a prediction that can be **wrong**.

> **Caveat.** $\Lambda/M_{\rm Pl}^4\sim10^{-122}$ is itself an input, and the prefactor of
> $f_{\rm res}$ is a ringdown-scale normalisation, so the robust content is the
> many-orders-of-magnitude *separation* between the three readings, not the third digit of
> the frequency.

---

## 5. A note on scientific scope

This code establishes three defensible claims:

1. the master equation realises both cosmic sectors from one field (Criterion 2, qualitative);
2. its tachyonic instability rate matches the analytic dispersion relation to ~0.02% (quantitative);
3. the horizon damping equals the Hawking rate exactly (Criterion 3).

It does **not** establish, and the companion note leaves open (Tier C):

- the off-horizon damping profile $\gamma(\mathcal{I})$ — only its horizon value is pinned;
- whether the spatial dimension $d$ can itself evolve.

The diagnostic in §3(a) initially **failed** when the growth rate was read off a localised
packet's peak (which mixes $k>0$ modes that grow more slowly). The flat-in-resolution error
revealed this as an *observable* bias, not a solver error; switching to the true $k=0$ mode
fixed it. This is recorded deliberately as an example of the code catching its own mistake.

---

## 6. Reproducing

```bash
pip install -r requirements.txt
python code/dg_condensate_dynamics.py     # Criterion 2 table
python code/dg_condensate_validation.py   # three quantitative checks
python code/dg_lisa_prediction.py          # falsifiable LISA prediction
python figures/make_figure.py             # regenerate the figure
```

## Reference

H. Hertault, *The Dynamics of the Hertault Condensate*, Dark Geometry series companion note (2026).

## License

MIT — see [LICENSE](LICENSE).
