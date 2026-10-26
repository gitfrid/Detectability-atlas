# Detectability Atlas

**Before you analyze observational data for a small effect, ask: could you even see it?**

This is a small, focused planning tool. It does not perform causal inference. It does not adjust for confounding. It answers one question, honestly and quantitatively:

> Given your sample size, baseline event rate, presumed confounding strength, observation window, and effect duration — **what is the smallest effect you could detect with 90% power?**

If the answer is "an effect larger than the one you're looking for", no amount of clever statistics will save the analysis. Better to know that before you start.

## What it computes

For an effect of size **1/X** (one additional event per X doses):

| Quantity | Meaning |
|---|---|
| **X_min(h, w, L)** | Detection threshold. If X ≤ X_min, the effect is detectable with 90% power. |
| **X_max(h, L)** | Structural bound. Protective effects stronger than this are impossible (you cannot prevent more deaths than the baseline produces). |
| **E-value contour** | Robustness bound. If an additional unmeasured confounder of this strength existed, it would fully explain your observed effect. |

The optimal observation window is **w = L** (matching window to effect duration).

## Why this matters

Most power calculations answer: "how large must N be to detect an effect of size θ?" This tool answers the inverse question, which is often more useful *before* an analysis: **"given my N, what θ can I detect — and is my target θ even in range?"**

For rare-event data with unmeasured confounding (vaccine safety, drug safety, policy evaluation), the answer is frequently: *no, not without controlling the confounding*. This tool makes that visible, in advance, in one figure.

## What this tool does **not** do

- It does **not** estimate causal effects.
- It does **not** adjust for confounding — it *parameterizes* it (`h`), because confounding strength is what you don't know.
- It does **not** replace Target Trial Emulation, negative controls, instrumental variables, or randomization.
- It is **not** a novel statistical method. It uses textbook power analysis for Poisson rates, with an HVE-adjusted baseline and a structural upper bound added.

## Honest limitations

| Limitation | Effect |
|---|---|
| **Uniform effect distribution over `L`** | Real effects often wane. The uniform assumption is *conservative* for waning (formula understates detectability) and *anti-conservative* for delayed-onset effects. |
| **Normal approximation** | Accurate for expected counts ≳ 10. For very small windows (E_V ≈ 1000), empirical power is ~5% below the analytical 90%. |
| **`h` must be specified** | The HVE strength is exactly what you don't know. The tool shows you *how detectability depends on it* — it cannot tell you which `h` is correct. |
| **`X_max` bound valid only for w ≤ L** | For w > L, the structural bound becomes `1/(μ w (1-h))`. |
| **Illustrative parameters** | `N_V = 6e6`, `N_U = 4e6`, `μ = 3e-5` are placeholder values. Calibrate to your data before drawing conclusions. |

## Who should use it

- Anyone planning an observational analysis with **rare events** and **likely unmeasured confounding** (e.g. Healthy Vaccinee Effect, frailty selection).
- want to check whether a claimed "no effect" is a real null or a power failure.
- If you want explain why "no statistically significant effect" is not the same as "no effect".


