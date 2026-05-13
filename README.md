# MRD ctDNA Calculator — Adult Lymphomas (CXJ1)

Web calculator implementing **Cox proportional hazards models** (V15.130 + V15.143)
for prediction of progression-free survival (PFS) at 12 and 24 months
after CXJ1 follow-up timepoint, based on circulating tumor DNA (ctDNA)
minimal residual disease (MRD) markers.

🔬 **Live calculator** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/>
📚 **Methodology & validation** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/methodology.html>

---

## Two Cox models (toggle inside calculator)

### M2 (V15.130) — with TEP (3 variables) — **default**
| Pipeline | N | Events | C-index (apparent) | C-index CV-5 | AIC |
|---|---|---|---|---|---|
| **IDES** (MOABI panel) | 143 | 31 | **0.859** | 0.877 | 238.59 |
| **PV** (phased variants) | 91 | 16 | **0.928** | 0.956 | 93.05 |

### M1 (V15.143) — without TEP (2 variables) — fallback when TEP unavailable
| Pipeline | N | Events | C-index (apparent) | AIC |
|---|---|---|---|---|
| **IDES** | 175 | 44 | 0.786 | 388.93 |
| **PV** | 108 | 22 | 0.826 | 161.19 |

Toggle "TEP non disponible" in the calculator switches between models.

---

## Model variables

1. **Score (joint bi-exponential)** — trajectory fit between diagnostic
   VAF and CXJ1 VAF, scored against good vs. bad responder reference
   trajectories. Personalized by histology (Hodgkin classique vs. rest→DLBCL).
2. **v5A_gated** — driver gene signal :
   `log1p(Σ reads_drv) × (n_drv_pos / n_total_pos) × quali_patient`
   with drivers panel D11 = {BCL2, BCL6, BCL7A, BTG2, CIITA, CXCR4,
   IRF8, MYC, PAX5, S1PR2, TP53}.
3. **TEP2_pos** — Deauville score at TEP interim cycle 2 (binary positive ≥ 4).
   *Only in M2.*

---

## Risk stratification (V15.132 — Tern 20/50)

Three clinical zones based on predicted r12 (12-month PFS risk):

| Zone | Range | IDES PFS@12 / @24 | PV PFS@12 / @24 |
|---|---|---|---|
| 🟢 **Faible** | r12 ≤ 20% | 92.7% / 92.7% | 95.7% / 95.7% |
| 🟡 **Modéré** | 20% < r12 ≤ 50% | 84.2% / **31.6%** | 90.0% / **45.0%** |
| 🔴 **Élevé** | r12 > 50% | 8.3% / 8.3% | 0% / 0% |

Note : the Modéré zone shows late relapses between 12 and 24 months → justifies close
surveillance even if PFS@12 looks reassuring.

---

## Key finding : baseline scores don't add value

Likelihood-ratio tests of adding **IPI** (DLBCL) or **Hasenclever IPS** (Hodgkin) to the M2 model:

| Score | Pipeline | LR test | Verdict |
|---|---|---|---|
| IPI | IDES (N=55) | χ²=1.97, **p=0.16** | NS |
| IPI | PV (N=40) | χ²=0.98, **p=0.32** | NS |
| Hasenclever | IDES (N=77 Hodgkin) | χ²=0.07, **p=0.79** | NS |
| Hasenclever | PV (N=43 Hodgkin) | χ²=1.25, **p=0.26** | NS |

→ **Mid-treatment molecular (ctDNA) + metabolic (TEP) response overrides baseline prognostic scores.**

---

## Files

| File | Description |
|---|---|
| `index.html` | Main calculator (standalone, no dependencies) |
| `methodology.html` | Methodology & validation page (performance, KM, calibration, baseline-score tests) |
| `figs/` | Embedded PNG figures for the methodology page |
| `LICENSE` | MIT |
| `.gitignore` | Excludes patient data (xlsx, csv, input/) |

---

## Usage

Open `index.html` (or visit the GitHub Pages URL) in any modern browser.
The calculator is **100% client-side**, no data is transmitted, no external dependencies.

---

## Citation

> Calculateur MRD CXJ1 (V15.130 + V15.143) — Cox proportional hazards models for ctDNA MRD
> in adult lymphomas after CXJ1 follow-up. Cohort development : MOABI panel (IDES) and
> phased variants (PV). Laboratoire d'immunologie biologique GHU Mondor — secteur biologie
> moléculaire (hmn-immuno-bm), 2026.

---

## Disclaimer

This calculator is a **research prototype** intended for academic and methodological
discussion. It is **not** a validated clinical decision support tool and must not be
used for individual patient management without independent clinical assessment by a
qualified hematologist.

---

## License

MIT — see [LICENSE](LICENSE).
