# Detectability Atlas

Before you analyze observational data for a small effect, ask three questions:
Could you even see it? Is it structurally possible? And if it is harmful, would the HVE mask it?

This tool answers all three, symmetrically, in two directions. It does not
estimate causal effects. It does not adjust for confounding. It quantifies what
observational data could detect, and what is impossible regardless of sample size.

---

## 1/5000 ON ALL-CAUSE MORTALITY: CAN THE DATA SHOW IT?

    Harmful:     No.  Masked by HVE for any plausible h.

    Protective:  No.  Structurally impossible for p_C <= 7%.

    Two different reasons. Both say: the data cannot resolve at this size.

(illustrative: h = 0.5, L = 180, p_C = 0.05, p_T = 0.20)

---

## The figure

![Three bounds on harmful and protective effects](output/detectability-atlas/atlas_w_eq_L.png)

Left (harmful):    X_min (detection), X_vis (masking), X_max_tox (structural).
Right (protective): X_min (detection), X_max_prot (structural).
Red dashed line:   the 1/5000 benchmark.

Assumes optimal window w = L. For w != L, thresholds fall; see next figure.
N_V, N_U, mu, p_C, p_T are illustrative. Calibrate to your data.

---

## Three bounds, two directions

For an effect of size 1/X (one additional or prevented event per X doses):

    X_min(h, w, L):        Detection threshold. Detectable at 90% power if
                           X <= X_min. Both directions.

    X_max_tox(p_T):        Structural bound. Harmful effects stronger than
                           1/p_T are impossible. Harmful only.

    X_vis(h, L; p_T):      Visibility bound. Harmful effects weaker than X_vis
                           are masked by HVE. Harmful only.

    X_max_prot(h, L; p_C): Structural bound. Protective effects stronger than
                           this are impossible. Protective only.

The bounds, explicit:

    Detection:       X_min = min(w,L) / [L*w*mu*(1-h)*(z_a+z_b)*SE(w)]
                     SE(w) = sqrt(1/(w*N_V*mu*(1-h)) + 1/(w*N_U*mu))

    Harmful:         X_max_tox ~ 1/p_T
                     X_vis     = p_T / (mu*h*L)

    Protective:      X_max_prot = 1 / (p_C*mu*(1-h)*L)

where:
    p_C = COVID-attributable fraction of ACM (~2-10%, observable)
    p_T = fraction of ACM causes a toxic effect could affect
          (~10-30%, biological: cardiovascular, thrombotic, infections;
          not accidents, suicide, most cancers)

---

## The visible range (illustrative: h=0.5, L=180, p_C=0.05, p_T=0.20)

    Protective:    7,407 <= X <= 9,866     (a 25% window in X)
    Harmful:           5 <= X <  74        (a factor of 15 in X)

These ranges are DISJOINT. Everything in the gap cannot be reached in
either direction — the data is blind to it.

Why p_C and p_T matter more than h:
    p_C is observable from the data (COVID death fraction).
    p_T is bounded by biology (which causes can a vaccine affect?).
    The structural bounds are more trustworthy than the detection
    threshold, which depends on the unknown h.

---

## Why w = L matters

![X_min(w) for various L](output/detectability-atlas/facets_w_vs_L.png)

X_min is maximized when the observation window w matches the effect
duration L. In practice L is unknown, so real-world detectability
lies BELOW the w = L curve. The hero figure is a best case.

Behavior:
    w < L:   X_min rises with w    (more effect visible)
    w = L:   MAXIMUM               (optimum)
    w > L:   X_min falls with w    (only background noise)

---

## Worked example: Czech FOI dataset

The tool was built to ask: could the Czech FOI data (~11M individuals,
2021-2022) distinguish a causal vaccine effect from HVE on ACM?

    Harmful:     1/5000 masked by any h >= 0.01.
                 Absence of signal is not evidence of absence.

    Protective:  1/5000 structurally impossible for p_C <= 7%.
                 Presence of signal is not evidence of cause.

Symmetric conclusion: for 1/5000 on ACM, observational data alone
cannot resolve the question, in either direction. The tool quantifies
why: the protective range is blocked by p_C, the harmful range by X_vis.

---

## What this tool does NOT do

- Does not estimate causal effects (neither direction).
- Does not adjust for confounding — it parameterizes it (h).
- Does not replace Target Trial Emulation, negative controls, IVs,
  or randomization.
- Is not a novel statistical method. Textbook power analysis for
  Poisson rates, with an HVE-adjusted baseline, a visibility bound,
  and two structural bounds added.

---

## Limitations (honest)

Uniform effect over L:
    Real effects often wane. Uniform assumption is conservative for
    waning, anti-conservative for delayed onset.

Normal approximation:
    Accurate for expected counts >= 10. For very small windows
    (E_V ~ 1000), empirical power is ~5% below target. Script prints
    CHECK in this case; expected, not a bug.

h is unknown:
    Tool shows how detectability depends on h. It cannot tell you
    which h is correct.

p_C, p_T must be specified:
    Observable (p_C) or biologically bounded (p_T). The choice matters;
    README shows sensitivity.

X_max_tox ~ 1/p_T is nearly constant in h:
    By design — it is a biology bound, not a confounding bound.

X_vis depends on h:
    The only bound that depends critically on h. Fundamental
    sensitivity of the analysis.

Parameters are illustrative:
    N_V=6e6, N_U=4e6, mu=3e-5, p_C=0.05, p_T=0.20.
    Calibrate before drawing conclusions.

---

## Who should use it

    - Planning observational analysis with rare events, aggregated
      endpoints (ACM), and likely unmeasured confounding (vaccine
      safety, drug safety, policy evaluation).
    - Anyone checking whether a claimed "no harmful effect" or
      "protective effect" is real or a power/confounding/structural
      artefact.
    - Anyone explaining that "no significant effect" is not the same
      as "no effect", and that "not enough power", "structurally
      impossible", "masked by confounding" are three different reasons.

---

## Status

Methodological note, not a framework.
Mathematics is correct and Monte-Carlo-verified in both directions.
Cite as a planning tool, not as a solution to confounding.
