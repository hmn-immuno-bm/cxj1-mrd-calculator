# MRD ctDNA Calculator — Adult Lymphomas (CXJ1)

Web calculator implementing **Cox proportional hazards models** for prediction of progression-free survival (PFS) at 12 and 24 months after mid-treatment timepoint (CXJ1), using circulating tumor DNA (ctDNA) minimal residual disease (MRD) markers in adult B-cell lymphomas and classical Hodgkin.

🔬 **Live calculator** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/>
📚 **Methodology** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/methodology.html>

## Overview

The calculator predicts post-treatment relapse risk from four covariates measured around the CXJ1 timepoint (typically end of cycle 2 or 3 of first-line chemo-immunotherapy):

1. **Kinetic score** — distance from "good responder" vs "bad responder" mono-exponential ctDNA decay trajectories, computed from VAF and cfDNA at diagnosis and CXJ1.
2. **Baseline tumor burden** — log(1 + ctDNA at diagnosis), capturing the prognostic weight of initial tumor mass.
3. **Driver gene signal (v5A_gated)** — residual NGS signal on the 11 oncogenic drivers (D11: BCL2, BCL6, BCL7A, BTG2, CIITA, CXCR4, IRF8, MYC, PAX5, S1PR2, TP53) at CXJ1, gated by a pipeline-specific quality filter.
4. **Interim FDG-PET** — Deauville score ≥ 4 at end of cycle 2 (binary).

## 8 model variants (auto-selected from inputs)

The calculator switches automatically between 8 variants based on data availability:

| Pipeline | Quantification | TEP | Variant | N | events | C-index |
|---|---|---|---|---|---|---|
| **IDES** | hEq (VAF × cfDNA) | ✅ | **M2_hEq** ⭐ | 134 | 24 | **0.860** |
| IDES | hEq | ❌ | M1_hEq | 158 | 30 | 0.807 |
| IDES | VAF % only | ✅ | M2_VAF | 135 | 24 | 0.846 |
| IDES | VAF % only | ❌ | M1_VAF | 160 | 30 | 0.775 |
| **PV** | hEq (VAF × cfDNA) | ✅ | **M2_hEq** ⭐ | 88 | 15 | **0.921** |
| PV | hEq | ❌ | M1_hEq | 103 | 20 | 0.841 |
| PV | VAF % only | ✅ | M2_VAF | 89 | 15 | 0.918 |
| PV | VAF % only | ❌ | M1_VAF | 105 | 20 | 0.822 |

⭐ = optimal configuration when all data is available.

**Measured impact of missing data:**
- Missing cfDNA Qubit (hEq → VAF): IDES ΔC-index −0.013 / PV −0.003 with TEP; −0.032 / −0.019 without TEP
- Missing PET (M2 → M1): IDES ΔC-index −0.053 / PV −0.080 with cfDNA; −0.071 / −0.096 without cfDNA

**Clinical baseline comparison** (single-covariate Cox on patient union, N=183, 34 events):

| Model | Covariates | C-index |
|---|---|---|
| TEP only (baseline) | 1 (TEP2_pos only) | **0.646** |
| M1_hEq MRD ctDNA only (no TEP) | 3 (score, burden, v5A) | **0.807 / 0.841** |
| M2_hEq MRD + TEP (full) | 4 | **0.860 / 0.921** |

**ctDNA MRD is the dominant prognostic signal**: adding ctDNA to TEP improves C-index by **+0.18** on average, while adding TEP to ctDNA improves it by only **+0.07**. The interim PET retains incremental value when combined with MRD, but cannot replace it.

## Pipelines

- **IDES** — MOABI hybrid-capture panel sequencing with Watch List approach; quality gate = Monte-Carlo p-value ≤ 1/100001 (10⁵ permutations).
- **PV** — Phased-variant UMI-based sequencing; quality gate = ≥ 3 positive doublets AND ≥ 3 total UMIs at CXJ1.

## Model specification

For each variant, the Cox model is fitted with partial pooling stratified by pipeline:
- `β_score`, `β_burden`, `β_TEP` are **shared** between IDES and PV (LR test p > 0.97 — supports biological universality of these markers).
- `β_v5A` is **pipeline-specific** (Reads_alt scale for IDES vs UMI scale for PV — LR test rejects sharing p < 0.005).
- Baseline hazard `S₀(t)` is **stratified by pipeline**.

### Linear predictor and risk

```
LP_raw = β_score · score_kin + β_burden · log_burden + β_v5A · v5A_gated + β_TEP · TEP2_pos
LP    = LP_raw − LP_offset            ← centering on cohort mean
S(t)  = S₀(t)^exp(LP)
r(t)  = 1 − S(t)
```

The **LP_offset** is the cohort-mean β·X̄ stored per variant in the data file. This centering matches lifelines' default behavior where `baseline_survival_` is evaluated at the mean of training covariates (not at zero).

### Production coefficients (M2_hEq, all data available)

Shared coefficients (identical for IDES and PV):
- β_score_kin = +0.470
- β_log_burden_diag = +0.298
- β_TEP2_pos = +0.911

Pipeline-specific:
- β_v5A_gated IDES = +0.590, PV = +1.491
- LP_offset IDES = +3.402, PV = +3.559
- S₀(12) IDES = 0.9048, PV = 0.9352

## Mono-exponential trajectories (4 parameters per mode)

ctDNA decay rates fitted by histology × responder status:

| Pipeline | Histology | λ_good (j⁻¹) | t½_good | λ_bad (j⁻¹) | t½_bad |
|---|---|---|---|---|---|
| IDES | Hodgkin | 0.649 | 1.1 d | 0.371 | 1.9 d |
| IDES | DLBCL | 0.388 | 1.8 d | 0.217 | 3.2 d |
| PV | Hodgkin | 0.558 | 1.2 d | 0.241 | 2.9 d |
| PV | DLBCL | 0.372 | 1.9 d | 0.241 | 2.9 d |

Histology is binarized: **Hodgkin classical** uses the Hodgkin trajectory; all other B-cell lymphomas (DLBCL, HGBL, PMBL, Burkitt, transformed indolent, etc.) use the DLBCL trajectory.

## Clinical risk stratification (Tern 10/40)

| Zone | r12 range | Action |
|---|---|---|
| 🟢 **Low** | r12 ≤ 10% | Standard surveillance |
| 🟡 **Intermediate** | 10% < r12 ≤ 40% | Closer surveillance (additional PET, regular MRD) |
| 🔴 **High** | r12 > 40% | Multidisciplinary discussion for additional therapy (consolidation, CAR-T, trial) |

Within the development cohort, the High zone reaches **0% PFS at 24 months** in both pipelines (universal relapse).

## Local calibration (adaptive KNN)

Alongside the Cox prediction, the calculator displays a **non-parametric KM estimate** computed on neighbors of the patient by LP. The KNN strategy is empirically validated:

1. Take all neighbors with `|LP_i − LP_query| ≤ 1.0`
2. If fewer than 8, expand to 8 nearest (sparse-tail safety)
3. If more than 40, truncate to 40 nearest

This reduces calibration error in the prediction tail (ACE_tail) by ~10× compared to a fixed K=30 neighborhood. Empty / sparse regions are flagged with a visual warning so the clinician knows the Cox model is extrapolating.

The KM 95% Greenwood log-log confidence interval is reported numerically below the survival plot.

## Validation

- **Linearity**: martingale residual lowess amplitude = 0.13 (< 0.3 threshold); Grambsch-Therneau PH test p > 0.48 for all covariates; LR test of adding RCS splines on log_burden p = 0.80 (non-significant).
- **Bootstrap optimism-corrected C-index** (B=500, Harrell): 0.860 apparent → 0.852 corrected (optimism = 0.009).
- **β_log_burden stability**: 95% percentile CI [0.10, 0.47], 100% of bootstraps positive.
- **Partial pooling validation**: LR test of sharing β_score + β_burden + β_TEP between IDES and PV: χ² = 0.20, df = 3, p = 0.978 (supports sharing). Sharing v5A: rejected, p = 0.005.
- **Predictive contribution per covariate** (LR test): log_burden p = 0.002 (M2_hEq), v5A_gated p < 0.001, TEP p = 0.001, score_kin p = 0.005.

## Confidence intervals (Wilson)

The calculator displays **95% Wilson confidence intervals** on observed NPV/PPV for the patient's risk group, addressing calibration uncertainty at extreme risk values (where N < 10 in development cohort).

## Privacy

All computations run client-side in the browser. No patient data is transmitted to any server. The model coefficients are loaded at startup and remain local.

## Citation

> Calculateur MRD CXJ1 — Cox proportional hazards models for ctDNA MRD in adult B-cell and Hodgkin lymphomas after mid-treatment. Partial pooling across IDES (MOABI hybrid-capture) and PV (phased-variant UMI) pipelines, with mono-exponential decay trajectories per histology × responder strata and integrated baseline tumor burden. Laboratoire d'immunologie biologique GHU Mondor — secteur biologie moléculaire (hmn-immuno-bm), 2026.

## License

MIT

## Disclaimer

The calculator is a **research prototype**. It is not validated for individual clinical care. Any therapeutic decision must rely on the independent evaluation of a qualified hematologist. Calibration at extreme risk values (r12 > 60% or ctDNA_diag < 500 hEq/mL) is limited by development cohort size. External validation on independent cohort is pending.
