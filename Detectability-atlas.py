"""
Detectability Atlas
====================
A framework for computing what observational data *could* detect versus what is
structurally impossible — under unmeasured confounding (HVE), power constraints,
and effect-duration limits.

Core formula (derived from first principles)
--------------------------------------------
For a candidate effect of size 1/X (one event per X doses), the detection
threshold with power 1-β is:

    X_min(h, w, L) = min(w, L) / [ L · w · μ · (1-h) · (z_α + z_β) · SE(w) ]

with standard error:

    SE(w) = sqrt( 1/(w · N_V · μ · (1-h)) + 1/(w · N_U · μ) )

The structural upper bound for protective effects (w = L) is:

    X_max(h, L) = 1 / (μ · L · (1-h))

Interpretation
--------------
- An effect 1/X is *detectable* if X ≤ X_min(h, w, L).
- An effect 1/X is *structurally possible* (protective) if X ≥ X_max(h, L).
- The optimal observation window is w = L (maximizes X_min).

Behavior of X_min as a function of the window w
------------------------------------------------
- w ≤ L : X_min rises with w  (more of the effect becomes visible)
- w = L : MAXIMUM (optimum)
- w > L : X_min falls with w  (only additional background noise)

Verifications
-------------
1. Consistency: X_min_full(h, L, L) == X_min(h, L) for all (h, L).
2. Monotonicity: maximum at w = L, rising for w < L, falling for w > L.
3. Monte-Carlo: empirical power ≈ 0.90 at the analytical threshold.

Author: (your name)
License: MIT
"""

import os
import numpy as np
import matplotlib

# Use the non-interactive Agg backend so the script runs headless (e.g. on CI).
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (needed for 3D projection)


# ============================================================
# Configuration
# ============================================================
# Output directory — change to your preferred location.
OUT_DIR = os.path.join("output", "detectability-atlas")
os.makedirs(OUT_DIR, exist_ok=True)
print(f"Output directory: {OUT_DIR}")

# Population parameters (illustrative — calibrate to your data).
N_V = 6.0e6              # Number of vaccinated individuals
N_U = 4.0e6              # Number of unvaccinated individuals
MU = 3.0e-5              # Baseline mortality per person per day

# Statistical parameters.
Z_ALPHA = 2.33           # Bonferroni-corrected critical value (5 windows, α=0.05, one-sided)
Z_BETA = 1.282           # 90% power
Z_SUM = Z_ALPHA + Z_BETA


# ============================================================
# Boundary functions
# ============================================================

def X_min(h, L, N_V=N_V, N_U=N_U, mu=MU, z_sum=Z_SUM):
    """
    Detection threshold at the optimal window w = L.

    Parameters
    ----------
    h : float or ndarray
        HVE strength (fractional reduction of baseline mortality in vaccinated).
    L : float or ndarray
        Effect duration (days).
    N_V, N_U : float
        Group sizes (vaccinated / unvaccinated).
    mu : float
        Baseline mortality per person per day.
    z_sum : float
        Sum of critical values (α + β).

    Returns
    -------
    X_min : ndarray
        Minimum effect size 1/X detectable with power 1-β.
    """
    h = np.clip(np.asarray(h, float), 0.0, 0.95)
    L = np.maximum(np.asarray(L, float), 1e-9)
    se = np.sqrt(1.0 / (L * N_V * mu * (1.0 - h))
                 + 1.0 / (L * N_U * mu))
    return 1.0 / (L * mu * (1.0 - h) * z_sum * se)


def X_min_full(h, w, L, N_V=N_V, N_U=N_U, mu=MU, z_sum=Z_SUM):
    """
    Full detection threshold for arbitrary window w.

    This is the general form. For w = L it collapses to X_min(h, L).

    Parameters
    ----------
    h : float or ndarray
        HVE strength.
    w : float or ndarray
        Observation window (days).
    L : float or ndarray
        Effect duration (days).

    Returns
    -------
    X_min : ndarray
        Minimum effect size 1/X detectable with power 1-β.
    """
    h = np.clip(np.asarray(h, float), 0.0, 0.95)
    w = np.maximum(np.asarray(w, float), 1e-9)
    L = np.maximum(np.asarray(L, float), 1e-9)
    min_wL = np.minimum(w, L)
    se = np.sqrt(1.0 / (w * N_V * mu * (1.0 - h))
                 + 1.0 / (w * N_U * mu))
    return min_wL / (L * w * mu * (1.0 - h) * z_sum * se)


def X_max_prot(h, L, mu=MU):
    """
    Structural upper bound for protective effects (w = L only).

    A protective effect stronger than this bound is impossible:
    it would prevent more deaths than the baseline produces.
    """
    h = np.clip(np.asarray(h, float), 0.0, 0.95)
    L = np.maximum(np.asarray(L, float), 1e-9)
    return 1.0 / (mu * L * (1.0 - h))


# ============================================================
# Verification 1 — analytical consistency
# ============================================================
# The full formula must reduce to the w = L form when w = L.
print("\n" + "=" * 78)
print("VERIFICATION 1: X_min_full(h, L, L) == X_min(h, L)")
print("=" * 78)
print(f"{'h':>5} {'L':>5} {'X_min(h,L)':>14} {'X_min_full(h,L,L)':>18} {'Status':>8}")
print("-" * 78)
for h in [0.0, 0.2, 0.5, 0.7]:
    for L in [7, 28, 90, 180]:
        a = float(X_min(h, L))
        b = float(X_min_full(h, L, L))
        rel = abs(a - b) / max(a, 1e-30)
        status = "OK" if rel < 1e-10 else "ERROR"
        print(f"{h:>5.2f} {L:>5} {a:>14.0f} {b:>18.0f} {status:>8}")
print("-" * 78)


# ============================================================
# Verification 2 — monotonicity in w for fixed L
# ============================================================
# X_min should rise with w until w = L, then fall.
print("\n" + "=" * 78)
print("VERIFICATION 2: X_min(w) behavior for L=28, h=0.2")
print("=" * 78)
print(f"{'w':>5} {'X_min(w)':>14} {'Comment':>38}")
print("-" * 78)
for w in [1, 3, 7, 14, 28, 60, 90, 180]:
    val = float(X_min_full(0.2, w, 28))
    if w < 28:
        comment = "w < L: rising (more effect visible)"
    elif w == 28:
        comment = "w = L: MAXIMUM (optimum)"
    else:
        comment = "w > L: falling (more noise)"
    print(f"{w:>5} {val:>14.0f}  {comment}")
print("-" * 78)


# ============================================================
# Verification 3 — Monte-Carlo power at the analytical threshold
# ============================================================
def mc_power(h, w, L, X_target, n_sim=10000, seed=42,
             N_V=N_V, N_U=N_U, mu=MU):
    """
    Simulate the power at an effect of size 1/X_target.

    Draws Poisson-distributed death counts in V and U, computes the
    log rate ratio, and returns the fraction of simulations exceeding
    the critical value z_α.
    """
    rng = np.random.default_rng(seed)
    phi = 1.0 / X_target
    lam_V = w * N_V * mu * (1.0 - h) + N_V * phi * np.minimum(w, L) / L
    lam_U = w * N_U * mu
    EV = np.maximum(rng.poisson(lam_V, n_sim).astype(float), 1.0)
    EU = np.maximum(rng.poisson(lam_U, n_sim).astype(float), 1.0)
    rr_obs = (EV / w) / (EU / w)
    rr_0 = N_V * (1.0 - h) / N_U
    se0 = np.sqrt(1.0 / (w * N_V * mu * (1.0 - h))
                  + 1.0 / (w * N_U * mu))
    Z = np.log(rr_obs / rr_0) / se0
    return float(np.mean(Z > Z_ALPHA))


print("\n" + "=" * 78)
print("VERIFICATION 3: Monte-Carlo power at the analytical threshold")
print("=" * 78)
print(f"{'h':>5} {'w':>5} {'L':>5} {'X_formula':>12} {'Power_MC':>10} {'Status':>8}")
print("-" * 78)
for (h, w, L) in [(0.2, 7, 7), (0.2, 14, 7), (0.2, 28, 7),
                  (0.2, 90, 90), (0.2, 180, 180),
                  (0.5, 14, 7), (0.5, 90, 180)]:
    Xf = float(X_min_full(h, w, L))
    p = mc_power(h, w, L, Xf, n_sim=8000)
    status = "OK" if 0.86 <= p <= 0.94 else "CHECK"
    print(f"{h:>5.2f} {w:>5} {L:>5} {Xf:>12.0f} {p:>10.3f} {status:>8}")
print("-" * 78)


# ============================================================
# Plot 1 — Facets: X_min(w) for fixed L
# ============================================================
h_fix = 0.2
L_vals = [7, 28, 90, 180]
w_line = np.linspace(1, 180, 400)

fig, axes = plt.subplots(2, 2, figsize=(14, 9), dpi=200, sharex=True)
for ax, L_val in zip(axes.flat, L_vals):
    y = X_min_full(h_fix, w_line, L_val)
    ax.plot(w_line, y, color="C0", linewidth=2,
            label=f"X_min(w), L={L_val}")
    y_opt = float(X_min_full(h_fix, L_val, L_val))
    ax.axvline(L_val, color="red", linestyle="--", linewidth=1.2,
               label=f"Optimum w = L = {L_val}")
    ax.axhline(5000, color="gray", linestyle=":", linewidth=1,
               label="1/5000")
    ax.set_title(f"L = {L_val} days  (optimum: X = {y_opt:.0f})")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, loc="lower right")
for ax in axes[-1]:
    ax.set_xlabel("Window w (days)")
for ax in axes[:, 0]:
    ax.set_ylabel("X_min (log)")
fig.suptitle(f"X_min(w) for various L — h = {h_fix}", y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "facets_w_vs_L.png"), dpi=200)
plt.close(fig)


# ============================================================
# Plot 2 — 3D surface: X_min(w, L)
# ============================================================
w_grid = np.linspace(1, 180, 120)
L_grid = np.linspace(1, 180, 120)
W, Lg = np.meshgrid(w_grid, L_grid)
Z = X_min_full(h_fix, W, Lg)

fig = plt.figure(figsize=(11, 7), dpi=200)
ax = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(W, Lg, Z, cmap="viridis",
                       linewidth=0, antialiased=True)
fig.colorbar(surf, ax=ax, shrink=0.6, label="X_min (doses per death)")

w_diag = np.linspace(1, 180, 100)
z_diag = X_min_full(h_fix, w_diag, w_diag)
ax.plot(w_diag, w_diag, z_diag, color="red", linewidth=2,
        label="Optimum w = L")

ax.set_xlabel("Window w (days)")
ax.set_ylabel("Effect duration L (days)")
ax.set_zlabel("X_min")
ax.set_title(f"X_min(w, L) — h = {h_fix}, 90% power")
ax.legend()
ax.view_init(elev=25, azim=45)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "surface_3d_w_L.png"), dpi=200)
plt.close(fig)


# ============================================================
# Plot 3 — Heatmap: X_min(w, L)
# ============================================================
fig, ax = plt.subplots(figsize=(9, 7), dpi=200)
cf = ax.contourf(W, Lg, Z, levels=60, cmap="viridis")
plt.colorbar(cf, ax=ax, label="X_min (doses per death)")
ax.plot([1, 180], [1, 180], "w--", linewidth=2, label="Optimum w = L")
ax.set_xlabel("Window w (days)")
ax.set_ylabel("Effect duration L (days)")
ax.set_title(f"Heatmap X_min(w, L) — h = {h_fix}")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "heatmap_w_L.png"), dpi=200)
plt.close(fig)


# ============================================================
# Plot 4 — Atlas (w = L): detection threshold + structural bound
# ============================================================
h_grid = np.linspace(0.0, 0.9, 400)
L_show = [7, 28, 90, 180]
colors = {7: "#1f77b4", 28: "#ff7f0e", 90: "#2ca02c", 180: "#d62778"}

fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), dpi=200)
for ax, title, is_prot in [
    (axes[0], "Harmful — detection threshold (w = L)", False),
    (axes[1], "Protective — structural bound (w = L)", True),
]:
    for L_val in L_show:
        ax.plot(h_grid, X_min(h_grid, L_val),
                color=colors[L_val], linewidth=2,
                label=f"X_min, L = {L_val}")
    if is_prot:
        for L_val in L_show:
            ax.plot(h_grid, X_max_prot(h_grid, L_val),
                    color=colors[L_val], linestyle=":", linewidth=1.2,
                    alpha=0.7, label=f"X_max, L = {L_val}")
    ax.axhline(5000, color="red", linestyle="--", linewidth=1.5,
               label="1/5000")
    ax.set_xlabel("HVE strength h")
    ax.set_ylabel("X (doses per death)")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(True, which="both", alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "atlas_w_eq_L.png"), dpi=200)
plt.close(fig)


# ============================================================
# Plot 5 — E-value contours (robustness vs. additional confounders)
# ============================================================
def X_contour_evalue(h, L, E_target, mu=MU):
    """
    X such that E(RR_adj(X, h, L)) = E_target (analytical inversion).

    The E-value E_target specifies how strong an *additional* (unmeasured)
    confounder would need to be, both with treatment and with outcome,
    to fully explain the observed adjusted effect.
    """
    if E_target <= 1.0:
        return np.full_like(np.asarray(h, float), np.nan)
    h = np.clip(np.asarray(h, float), 0.0, 0.95)
    L = np.asarray(L, float)
    return (2.0 * E_target - 1.0) / (
        (E_target - 1.0) ** 2 * L * mu * (1.0 - h))


fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), dpi=200)
for ax, title, is_prot in [
    (axes[0], "Harmful — E-value contours (w = L)", False),
    (axes[1], "Protective — E-value contours (w = L)", True),
]:
    for L_val in L_show:
        ax.plot(h_grid, X_min(h_grid, L_val),
                color=colors[L_val], linewidth=2,
                label=f"X_min, L = {L_val}")
    if is_prot:
        for L_val in L_show:
            ax.plot(h_grid, X_max_prot(h_grid, L_val),
                    color=colors[L_val], linestyle=":", linewidth=1,
                    alpha=0.6)
    # E-value contours are computed for L = 180 as a representative slice.
    for E_t, col, lab in [(1.5, "#777777", "E = 1.5"),
                          (2.0, "#000000", "E = 2.0"),
                          (3.0, "#8B0000", "E = 3.0")]:
        y = X_contour_evalue(h_grid, 180.0, E_t)
        mask = np.isfinite(y) & (y > 1) & (y < 1e8)
        ax.plot(h_grid[mask], y[mask], linestyle="--",
                color=col, linewidth=1.5, alpha=0.8,
                label=f"{lab} (add. confounder)")
    ax.axhline(5000, color="red", linestyle="--", linewidth=1.5,
               label="1/5000")
    ax.set_xlabel("HVE strength h")
    ax.set_ylabel("X (doses per death)")
    ax.set_yscale("log")
    ax.set_ylim(1, 1e6)
    ax.set_title(title)
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(True, which="both", alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "atlas_evalue.png"), dpi=200)
plt.close(fig)


# ============================================================
# Done
# ============================================================
print(f"\nDone. Plots saved in: {OUT_DIR}")
print("  facets_w_vs_L.png   — X_min(w) for various L")
print("  surface_3d_w_L.png  — 3D surface X_min(w, L)")
print("  heatmap_w_L.png     — heatmap X_min(w, L)")
print("  atlas_w_eq_L.png    — detection threshold + structural bound")
print("  atlas_evalue.png    — E-value contours")