# MRD ctDNA Calculator — Adult Lymphomas (CXJ1)

Web calculator implementing the **V15.130 Cox proportional hazards model**
for prediction of progression-free survival (PFS) at 12 and 24 months
after CXJ1 follow-up timepoint, based on circulating tumor DNA (ctDNA)
minimal residual disease (MRD) markers.

🔬 **Live calculator** : https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/

## Pipelines supported

The calculator implements two parallel pipelines, both calibrated on
the same lymphoma cohort:

| Pipeline | N | Events | C-index (apparent) | C-index CV-5 | AIC |
|---|---|---|---|---|---|
| **IDES** (MOABI panel) | 143 | 31 | 0.859 | 0.877 | 238.59 |
| **PV** (variants phasés) | 91 | 16 | 0.928 | 0.956 | 93.05 |

## Model variables (3 predictors)

1. **Score (joint bi-exponential)** — trajectory fit between diagnostic
   VAF and CXJ1 VAF, scored against good vs. bad responder reference
   trajectories. Personalized by histology (Hodgkin classique vs. rest).
2. **v5A_gated** — driver gene signal :
   `log1p(Σ reads_drv) × (n_drv_pos / n_total_pos) × quali_patient`
   with drivers panel D11 = {BCL2, BCL6, BCL7A, BTG2, CIITA, CXCR4,
   IRF8, MYC, PAX5, S1PR2, TP53}.
3. **TEP2_pos** — Deauville score at CXJ1 timepoint (binary positive).

## Risk stratification (V15.132 — Tern 20/50)

Three clinical zones based on predicted r12 :

| Zone | Range | IDES PFS@12 | PV PFS@12 |
|---|---|---|---|
| 🟢 Faible | r12 ≤ 20% | 92.7% | 95.7% |
| 🟡 Modéré | 20% < r12 ≤ 50% | 84.2% | 90.0% |
| 🔴 Élevé | r12 > 50% | 8.3% | 0% |

## Technical VPN curve (V15.50 LIS + Hill smoothing)

The calculator also displays a depth-dependent technical Negative
Predictive Value (NPV) curve, indicating the reliability of a negative
MRD result given the sequencing depth used. Methodology :
Longest Increasing Subsequence (LIS) on empirical cohort data, then
smooth Hill-function fit with high-end envol to the 100% asymptote.

## Usage

Open `index.html` (or visit the GitHub Pages URL) in any modern browser.
The calculator is **100% client-side**, no data is transmitted, no
external dependencies. Inputs are validated locally.

## Citation

> Calculateur MRD CXJ1 V15.102 — Modèle pronostique pour la maladie
> résiduelle ctDNA dans les lymphomes adultes après CXJ1. Cohort
> development : MOABI panel + variants phasés (V15.130).

## Disclaimer

This calculator is a **research prototype** intended for academic and
methodological discussion. It is **not** a validated clinical decision
support tool and must not be used for individual patient management
without independent clinical assessment by a qualified hematologist.

## License

MIT — see LICENSE file.
