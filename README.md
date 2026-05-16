# ODE Model of an Inflammatory Feedback Network

A computational implementation and dynamical-systems analysis of a coupled
cell–cytokine model of inflammation. The model reproduces the qualitative
**bistability** between a *resolved* and a *chronic* inflammatory state described
in published ODE models of acute inflammation, and provides a reusable framework
for bifurcation analysis and global sensitivity analysis of cytokine-driven
disease dynamics.

**Stack:** Python · SciPy · SALib · NumPy · Matplotlib

---

## Overview

The model couples an activated immune-cell population to three cytokines —
TNF-α, IL-6, and IL-10. The pro-inflammatory cytokines (TNF-α, IL-6) recruit and
activate immune cells, which in turn produce more cytokines (**positive
feedback**). IL-10 provides anti-inflammatory **negative feedback** by
suppressing both recruitment and pro-inflammatory cytokine production.

The competition between these feedback loops produces bistability: over a range
of parameters the system has two coexisting stable steady states, separated by
an unstable saddle equilibrium. Which state the system settles into is set by
the strength of the initial immune challenge.

The notebook performs six analyses:

1. **Trajectory simulation** — resolved vs. chronic outcomes from identical parameters.
2. **Equilibrium analysis** — equilibria located by root-finding, classified by Jacobian eigenvalues.
3. **Bifurcation diagram** — numerical continuation revealing a hysteresis loop bounded by two saddle-node bifurcations.
4. **Basins of attraction** — the separatrix mapped in the plane of initial conditions.
5. **Global sensitivity analysis** — Sobol indices identifying the parameters that most strongly drive the chronic outcome.
6. **Two-parameter regime map** — delineating monostable and bistable regions.

---

## Model

State vector `y = [M, T, I6, I10]`:

| Symbol | Meaning |
|--------|---------|
| `M`    | activated immune-cell population (e.g. macrophages) |
| `T`    | TNF-α concentration (pro-inflammatory) |
| `I6`   | IL-6 concentration (pro-inflammatory) |
| `I10`  | IL-10 concentration (anti-inflammatory) |

Governing equations, with activating Hill function
`h⁺(x,k) = xⁿ / (kⁿ + xⁿ)` and inhibitory Hill function
`h⁻(x,k) = kⁿ / (kⁿ + xⁿ)`:

```
dM/dt   = [ s_M + k_MT·h⁺(T,K_MT) + k_MI6·h⁺(I6,K_MI6) ]·h⁻(I10,K_M10)·(1 − M/M_max) − d_M·M
dT/dt   = p_T·M·h⁻(I10,K_T10)   − d_T·T
dI6/dt  = p_I6·M·h⁻(I10,K_I610) − d_I6·I6
dI10/dt = p_I10·M·h⁺(I6,K_106)  − d_I10·I10
```

The shared Hill coefficient `n` sets the cooperativity of the feedback switches
and is the key structural parameter enabling bistability.

### Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `s_M`     | 0.01  | baseline (cytokine-independent) cell recruitment |
| `k_MT`    | 1.20  | recruitment driven by TNF-α |
| `k_MI6`   | 0.50  | recruitment driven by IL-6 |
| `K_MT`    | 1.00  | half-saturation, TNF-α recruitment |
| `K_MI6`   | 1.50  | half-saturation, IL-6 recruitment |
| `d_M`     | 0.15  | cell clearance / apoptosis rate |
| `M_max`   | 8.00  | carrying capacity of activated cells |
| `p_T`     | 1.20  | TNF-α production per activated cell |
| `K_T10`   | 0.60  | IL-10 half-saturation suppressing TNF-α |
| `d_T`     | 1.10  | TNF-α decay rate |
| `p_I6`    | 0.90  | IL-6 production per activated cell |
| `K_I610`  | 0.70  | IL-10 half-saturation suppressing IL-6 |
| `d_I6`    | 0.80  | IL-6 decay rate |
| `p_I10`   | 0.80  | IL-10 production per activated cell |
| `K_106`   | 1.20  | IL-6 half-saturation driving IL-10 production |
| `d_I10`   | 0.90  | IL-10 decay rate |
| `K_M10`   | 0.80  | IL-10 half-saturation suppressing recruitment |
| `n_hill`  | 3.0   | Hill coefficient (cooperativity; structural) |

Parameter values are illustrative and chosen to place the system in the
biologically relevant bistable regime.

---

## Key results

**Bistability.** At the default parameters the system has three coexisting
equilibria:

| Equilibrium | M | T | I6 | I10 | Stability |
|-------------|------|------|------|------|-----------|
| Resolved    | 0.070 | 0.076 | 0.079 | 0.000 | stable node |
| Saddle      | 0.259 | 0.283 | 0.292 | 0.003 | unstable (separatrix) |
| Chronic     | 1.819 | 0.815 | 1.075 | 0.676 | stable focus |

Stability is confirmed by the eigenvalues of the numerically-evaluated Jacobian.
The chronic state is a stable **focus** — its complex eigenvalues produce the
damped cytokine oscillations seen before the system locks into chronic
inflammation.

**Hysteresis.** Sweeping the TNF-α production rate `p_T` reveals a bistable
window at **`p_T ∈ [0.457, 1.760]`**, bounded by two saddle-node (fold)
bifurcations. Outside this window the system is monostable — resolved for low
`p_T`, chronic for high `p_T`.

**Sensitivity / candidate therapeutic targets.** Variance-based global
sensitivity analysis (Sobol indices, Saltelli sampling, all 17 kinetic
parameters varied ±40%) ranks the IL-10-mediated suppression parameters as the
dominant drivers of chronic-inflammation severity:

| Rank | Parameter | S₁ | Sₜ | Interpretation |
|------|-----------|------|------|----------------|
| 1 | `K_T10`  | 0.256 | 0.311 | IL-10 potency against TNF-α |
| 2 | `K_I610` | 0.103 | 0.271 | IL-10 potency against IL-6 (large interaction term) |
| 3 | `p_T`    | 0.092 | 0.137 | TNF-α production rate |
| 4 | `d_T`    | 0.096 | 0.135 | TNF-α decay rate |

The large gap between total-order (Sₜ) and first-order (S₁) indices for `K_I610`
indicates it acts predominantly through interactions with other parameters.

---

## Usage

The notebook is designed for Google Colab (no local setup required):

1. Open [Google Colab](https://colab.research.google.com/) → **File → Upload notebook**.
2. Upload `Inflammatory_Network_ODE_Model.ipynb`.
3. **Runtime → Run all.**

`SALib` is installed by the first cell; all other dependencies are pre-installed
in Colab. Use the **standard CPU runtime** — the workload is sequential ODE
integration and does not benefit from a GPU. Full execution takes roughly
3–4 minutes.

To run locally instead:

```bash
pip install numpy scipy matplotlib SALib
jupyter notebook Inflammatory_Network_ODE_Model.ipynb
```

### Performance note

The global sensitivity analysis (Cell 7) runs at `N = 128`
(~4,600 model evaluations, ~75 s). Increase to `N = 256` for tighter bootstrap
confidence intervals at the cost of roughly 5 minutes of runtime.

---

## Scope and limitations

- The model reproduces the **qualitative** bistability of published inflammation
  models; it is **not** fitted to a specific experimental dataset, and the
  parameter values are illustrative rather than calibrated. Quantitative claims
  about a particular biological system would require parameter estimation
  against measured cytokine time-courses.
- Bistability depends on sufficient feedback cooperativity. Reducing the Hill
  coefficient `n_hill` narrows or eliminates the bistable window.
- The bifurcation diagram is traced by integration- and root-finding-based
  continuation rather than a dedicated continuation package (e.g. AUTO, MatCont);
  this is adequate for the folds present here but does not perform formal
  branch-point detection.

The framework is modular — the cell–cytokine structure, Hill-type feedback,
continuation routine, and Sobol analysis can be re-parameterised for other
systems in which cytokine dynamics drive disease progression, including
organ-on-chip models.

---

## References

1. Reynolds A. et al. (2006). A reduced mathematical model of the acute inflammatory response. *Journal of Theoretical Biology* 242(1): 220–236.
2. Day J. et al. (2006). A reduced mathematical model of the acute inflammatory response II. *Journal of Theoretical Biology* 242(1): 237–256.
3. Herald M.C. (2010). General model of inflammation. *Journal of Theoretical Biology* 264(3): 935–945.
4. Saltelli A. et al. (2010). Variance based sensitivity analysis of model output. *Computer Physics Communications* 181(2): 259–270.
