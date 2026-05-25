# MRD ctDNA Calculator — Adult Lymphomas (CXJ1)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version V230.4](https://img.shields.io/badge/version-V230.4-blue.svg)](https://github.com/hmn-immuno-bm/cxj1-mrd-calculator/releases)
[![DOI (Zenodo)](https://img.shields.io/badge/DOI-pending%20(Zenodo)-lightgrey.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Cite this software](https://img.shields.io/badge/cite-CITATION.cff-blueviolet.svg)](./CITATION.cff)

Web calculator implementing **Cox proportional hazards models** for prediction of **Event-Free Survival (EFS = relapse OR progression)** at 12 and 24 months after the mid-treatment timepoint (CXJ1), using circulating tumor DNA (ctDNA) minimal residual disease (MRD) markers in adult B-cell lymphomas and classical Hodgkin.

🔬 **Live calculator** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/>
📚 **Methodology** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/methodology.html>
🛠️ **Pipeline source** : <https://github.com/alessiocg/cxj1-mrd-lymphoma-pipeline>
📜 **TRIPOD+AI checklist** : [TRIPOD_AI_checklist.md](./TRIPOD_AI_checklist.md)
📑 **Citation file** : [CITATION.cff](./CITATION.cff)

## Overview

The calculator predicts post-CXJ1 EFS risk from four covariates measured around the CXJ1 timepoint (typically end of cycle 2 or 3 of first-line chemo-immunotherapy):

1. **Kinetic score** — Gaussian log-likelihood distance between observed and predicted ctDNA decay ratios under "good responder" vs "bad responder" mono-exponential trajectories (fitted per histology × responder strata).
2. **Baseline tumor burden** — log(1 + ctDNA at diagnosis) where ctDNA = VAF × cfDNA in hEq mode, or log(1 + VAF at diagnosis) in VAF mode.
3. **Driver gene signal (v5A_gated) — HISTOLOGY-SPECIFIC PANEL since V230.4 (May 2026)** — residual NGS signal on the appropriate driver panel at CXJ1, gated by a pipeline-specific quality filter. The panel is selected automatically based on the patient's histology:
   - **Hodgkin classique** : 16 genes — `BTG2, BZRAP1, CD83, CIITA, CXCR4, DTX1, IRF8, ITPKB, PAX5, RHOH, S1PR2, SOCS1, STAT6, TP53, ZCCHC7-GRHPR, ZFP36L1`. Axis CXCR4 + JAK-STAT (Reed-Sternberg signature) + B-cell signaling.
   - **Non-Hodgkin (DLBCL, FL, MZL, MCL…)** : 16 genes — `BCL2, BIRC3, BTG1, BZRAP1, CIITA, DTX1, FOXO1, HIST1H1E, KLF2, LTB, MYC, PAX5, PIM1, POU2AF1, S1PR2, TMSB4X`. Axis B-cell signaling + chromatin + anti-apoptosis (DLBCL signature).
   - **Overlap**: 5 shared genes (`BZRAP1, CIITA, DTX1, PAX5, S1PR2`) — 11 distinct genes per histology.
   - **Methodology** : exhaustive backward elimination from whitelist of non-rotten genes (Ig V(D)J genes excluded ; HR<1 or NS univariate genes blacklisted). Reveals that the legacy D11 panel (BCL2, BCL6, …) was an NH-biased signature, with several D11 genes (BCL2, BCL6, BIRC3) actually anti-prognostic in Hodgkin.
4. **Interim FDG-PET** — Deauville score (linear, 1-5) at end of cycle 2, with **ΔSUVmax% fallback** if Deauville unavailable.

**Endpoint**: EFS (Event-Free Survival = relapse OR progression). Non-lymphoma deaths are censored (3 patients : 2 COVID, 1 hemorrhagic ulcer under transfusion contraindication). Justification in methodology §XI.10.

## 6 model variants (auto-selected from inputs)

The calculator switches automatically between 6 variants based on data availability:

| Variant | Quantification | TEP | N stacked (events) | **C-corrected V230.4** | Δ vs V230.2 |
|---|---|---|---|---|---|
| **M2_hEq Deauv** ⭐ | hEq (VAF × cfDNA) | Deauville 1-5 | 230 (42) | **0.862** | +0.009 |
| M2_hEq Delta | hEq (VAF × cfDNA) | ΔSUVmax% (fallback) | 200 (40) | **0.859** | +0.012 |
| M1_hEq | hEq (VAF × cfDNA) | — absent | 273 (53) | **0.845** | +0.018 |
| M2_VAF Deauv | VAF % only | Deauville 1-5 | 230 (42) | **0.832** | +0.024 |
| M2_VAF Delta | VAF % only | ΔSUVmax% (fallback) | 200 (40) | **0.834** | +0.020 |
| M1_VAF | VAF % only | — absent | 273 (53) | **0.801** | +0.020 |

⭐ = optimal configuration when all data is available. **C-index optimism-corrected (bootstrap V230.4 B=500, Harrell), optimism médian +0.007 (max +0.013) — excellent generalization.** **V230.4 update (May 2026)** : histology-specific v5A driver panels (16 H + 16 NH genes, backward elimination from whitelist of non-rotten genes). V230.4 dominates V230.2 on all 6 variants (ΔC_corrected +0.009 to +0.024).

**Measured impact of missing data** (Δ C-index, M2_hEq Deauv ⭐ vs alternatives):

- Missing cfDNA Qubit (hEq → VAF): IDES −0.023 / PV −0.071 in M2 ; IDES −0.032 / PV −0.046 in M1.
- Missing Deauville (M2 → M1): IDES −0.032 / PV −0.034 in hEq ; IDES −0.041 / PV −0.009 in VAF.
- Substituting Deauville by ΔSUVmax% (M2_Deauv → M2_Delta): equivalent (Δ < 0.02 in all variants); cohort shrinks by 13% (loss of patients with Deauville but no SUVmax baseline).

**Clinical baseline comparison** (TEP-only Cox, on union of patients with TEP available, N=183, 34 events):

| Model | Covariates | C-index (corr.) |
|---|---|---|
| TEP only (clinical baseline) | 1 (TEP2_pos binary) | **0.65** |
| M1_hEq (MRD only, no TEP) | 3 (score, burden, v5A) | **0.81 / 0.85** |
| M2_hEq Deauv ⭐ (MRD + TEP) | 4 (score, burden, v5A, Deauville) | **0.84 / 0.87** |

**ctDNA MRD is the dominant prognostic signal**: adding MRD to TEP improves C-index by **+0.20**, while adding TEP to MRD improves it only by **+0.03 to +0.04**. Interim PET retains incremental value when combined with MRD but cannot replace it.

## Pipelines

- **IDES** — MOABI hybrid-capture panel sequencing with Watch List approach; quality gate = Monte-Carlo p-value ≤ 1/100001 (10⁵ permutations).
- **PV** — Phased-variant UMI-based sequencing; quality gate = ≥ 3 positive doublets AND ≥ 3 total UMIs at CXJ1.

The PV cohort is **strictly nested in the IDES cohort** (PV ⊂ IDES, 92/138 in M2, 109/164 in M1) — when both pipelines are available, PV is treated as an internal sensitivity analysis on a subset filtered by stricter PV-side QC.

## Model specification

For each variant, the Cox model is fitted with **partial pooling B′** stratified by `(pipeline × histology)` since V230.2 (4 baselines h0 per variant: IDES×Hodgkin, IDES×Non-Hodgkin, PV×Hodgkin, PV×Non-Hodgkin; `cluster=NOM`, L2 penalizer = 0.05):

- `β_score`, `β_burden`, `β_TEP` are **shared** between IDES and PV (LR test against fully separated model : p > 0.97).
- `β_v5A` is **pipeline-specific** (Reads_alt scale for IDES vs UMI scale for PV — LR test rejects pooling p < 0.01).
- Baseline hazard `S₀(t)` is **stratified by pipeline**.

Validation of B′ on the V230.1 stacked cohort (IDES+PV, N=230 in M2_hEq Deauv): AIC = 305.79, lowest among all 16 possible pooling combinations (cf. methodology §VIII).

### Linear predictor and risk

```
LP_raw = β_score · score_kin + β_burden · log_burden + β_v5A · v5A_gated  [+ β_TEP · TEP_var]
LP     = LP_raw − LP_offset            ← centering on cohort mean
S(t)   = S₀(t)^exp(LP)
r(t)   = 1 − S(t)
```

The **LP_offset** is the cohort-mean `β · X̄` stored per variant × pipeline in the data file. This centering matches `lifelines.predict_survival_function` (which evaluates baseline at the mean of training covariates, not at zero).

### Production coefficients (M2_hEq Deauv ⭐)

Shared (identical for IDES and PV):
- β_score_kin = **+0.374**
- β_log_burden_hEq = **+0.264**
- β_Deauville_num = **+0.365** (HR per unit = 1.44 ; HR Deauville 5 vs 1 = 4.31)

Pipeline-specific:
- β_v5A_gated IDES = **+0.597**, PV = **+1.500**
- LP_offset IDES = **+3.920**, PV = **+4.094**
- S₀(12) IDES = 0.890, PV = 0.924

## Mono-exponential trajectories (4 parameters per mode × pipeline)

ctDNA decay rates fitted by histology × responder status (good = no EFS event @ 12 mo ; bad = EFS event @ 12 mo). Fitted by `differential_evolution` on log-ratio observed vs predicted.

**hEq mode** (ratio ctDNA absolute = VAF × cfDNA):

| Pipeline | Histology | λ_good (j⁻¹) | t½_good | λ_bad (j⁻¹) | t½_bad |
|---|---|---|---|---|---|
| IDES | Hodgkin | 0.646 | 1.1 d | 0.371 | 1.9 d |
| IDES | DLBCL | 0.369 | 1.9 d | 0.211 | 3.3 d |
| PV | Hodgkin | 0.415 | 1.7 d | 0.364 *(fallback DLBCL/good)* | 1.9 d |
| PV | DLBCL | 0.364 | 1.9 d | 0.237 | 2.9 d |

Histology is binarized: **Hodgkin classical** uses the Hodgkin trajectory; all other lymphomas (DLBCL, HGBL, PMBL, Burkitt, transformed indolent, NLPHL, EBV+, etc.) use the DLBCL trajectory.

**Fallback for Hodgkin/bad PV → DLBCL/good** (V229.1 decision) : the PV cohort has 0 calibration patients in the Hodgkin × bad-responder cell (7 Hodgkin classical patients with EFS=1, all with VAF_diag undetectable in PV). The chosen fallback is `λ_bad_Hodgkin_PV = λ_good_DLBCL_PV` (≈ 0.36 j⁻¹), consistent with IDES where empirically `λ_bad_Hodgkin ≈ λ_good_DLBCL ≈ 0.37`.

## Clinical risk stratification (Tern 15/75, V230.1)

| Zone | r12 range | Suggested action |
|---|---|---|
| 🟢 **Low** | r12 ≤ 15% | Standard surveillance |
| 🟡 **Intermediate** | 15% < r12 ≤ 75% | Closer surveillance (additional PET, regular MRD) |
| 🔴 **High** | r12 > 75% | Mandatory MDT discussion (consolidation, CAR-T, trial) |

Thresholds **15/75** selected by exhaustive scan with **round-number candidates only** (multiples of 5%) on the V230.1 EFS cohort (V15.332). 15/75 is the **only round threshold** that simultaneously :
- Maintains the correct ordering of the 3 KM curves at **BOTH 12 months AND 24 months** (Green > Orange > Red) across all 12 KMs (6 variants × 2 pipelines),
- Provides a Faible-Modéré gap of **18 percentage points at 12 months** and **32 percentage points at 24 months**,
- Keeps the High zone universally marked (EFS@24 = 0 % in every variant ; N_E ≥ 3-7 per cohort).

## Local calibration (adaptive KNN)

Alongside the Cox prediction, the calculator displays a **non-parametric KM estimate** computed on the patient's LP-neighbors:

1. Take all neighbors with `|LP_i − LP_query| ≤ 1.0`
2. If fewer than 8, expand to the 8 nearest (sparse-tail safety)
3. If more than 40, truncate to the 40 nearest

Validated empirically (V208 leave-one-out): ACE_tail divided by ~10× vs a fixed K=30. Empty/sparse neighborhoods are flagged with a visual warning so the clinician knows the Cox model is extrapolating.

The KM 95% Greenwood log-log confidence interval is reported numerically below the survival plot.

## PV vs IDES on overlap cohort N=92

When comparing both pipelines on the **same 92 patients** (M2_hEq Deauv overlap, PV ⊂ IDES strict), the C-index is essentially identical (IDES 0.893 vs PV 0.890, Δ = −0.002). The per-variant C-index gap (0.844 vs 0.890) reported in the variants table comes mostly from a **cohort-selection effect** (PV excludes 46 IDES-positive patients without PV-detectable variants — typically harder-to-predict cases), not from a discriminatory edge of the PV pipeline.

PV's real advantages on the same cohort :

| Metric (overlap N=92) | IDES | PV | Δ |
|---|---|---|---|
| Calibration slope (target = 1.0) | 1.43 | **1.09** | −0.34 |
| Brier R² @ 12mo (higher = better) | 0.47 | **0.57** | +0.10 |
| Continuous NRI (PV vs IDES) | — | +0.97 | ⭐ |
| IDI | — | +0.025 | + |
| HR per 1 SD of LP | 3.98 | 5.39 | +1.41 |
| **PPV zone Élevé (r12 ≥ 75%)** | 67 % | **100 %** | +33 pt |
| **Specificity zone Élevé** | 92 % | **100 %** | +8 pt |

→ PV is more parsimonious (9 vs 17 patients in the High zone) but each flag is certified by an EFS event. See methodology §VI.bis.

## Validation (V230.2, May 2026)

**V230.2** — minimal architectural update over V230.1, preserving all coefficients and adding only :
1. **Enlarged v5A panel** : D11 (11 genes) → D11 + STAT6 + SOCS1 (13 genes, +2 Hodgkin-specific JAK-STAT drivers)
2. **Strata h0 by `(pipeline × histo)`** : 4 baselines per variant instead of 2

Performance vs V230.1 (cohort M2_hEq Deauv stacked N=230, 42 events):
- ΔAIC = −42.6 (large fit gain from histo-stratified baseline)
- ΔC_all ≈ 0 (neutral discrimination)
- VPP zone Élevé Hodgkin (r12 ≥ 10%) : **16% → 29%** (clinical gain)
- All coefficients remain significant (β_score p=0.046, β_burden p=0.017, β_v5A_IDES p<0.001, β_v5A_PV p<0.001, β_Deauv p=0.006)



- **Bootstrap optimism-corrected C-index** (Harrell B=500): M2_hEq Deauv IDES ⭐ apparent 0.845 → corrected **0.835** (optimism = +0.010, IC95% percentile [0.76 – 0.92]) ; PV apparent 0.893 → corrected **0.873** (optimism = +0.020, IC95% [0.80 – 0.97]). Max optimism across 12 fits = +0.027, median +0.015.
- **Linearity**: martingale residual lowess amplitude = 0.13 (< 0.3 threshold); Grambsch-Therneau PH test p > 0.18 for all covariates; LR test for adding RCS splines on log_burden p = 0.80 NS.
- **β_log_burden stability**: 95% percentile CI [0.10, 0.47], 100% of bootstraps positive.
- **Partial pooling B′ validation**: LR test of sharing β_score + β_burden + β_TEP between IDES and PV: χ² ≈ 0, df = 3, p = 1.00 ✓. Sharing v5A rejected (LR p=0.041).
- **Predictive contribution per covariate** (LR test, M2_hEq Deauv): log_burden p = 0.002 ; v5A_gated p < 0.001 ; Deauville_num p < 0.001 ; score_kin p ≈ 0.005.
- **Independence vs IPI / Hasenclever** (V230.1, V15.296): on the 6 variants, IPI and Hasenclever are **absorbed by LP V230.1** (all p_score+ > 0.09 in DLBCL, > 0.56 in Hodgkin). Reciprocally, LP V230.1 adds significant information beyond clinical scores (DLBCL p_LP+ < 0.001 on 6/6 variants; Hodgkin p_LP+ ≤ 0.012 on 6/6).

## V230.1 patches (technical reintegration of 2 patients)

V230.1 = V230 + 2 patients re-integrated after audit (V15.289-291):

- **Patient A** (DLBCL bad responder): inconsistent MOABI WL prefix between C2J1 and C3J1 (different 3-letter prefixes due to a compound surname). The default parser only matched one prefix and silently lost the other timepoint file. Patch: an internal alias module documents the prefix mapping (NAS files untouched). Audit confirmed this was the only such case among 191 WL_NGS files.
- **Patient B** (Hodgkin bad responder): C1J1 date error in the clinical database (a copy-paste typo placed the C3J1+1d date in the C1J1 cell ; the true C1J1 was confirmed by clinical re-review). With the corrected date, `delai_real(C2J1) = +21 d` (instead of negative) and the patient correctly enters the M1 IDES cohort. Patch: an internal clinical-patches module applies the correction in memory (NAS clinical file untouched).

**Impact**: +2 M1 IDES (162→164), +1 M1 PV (108→109), +2 M2 IDES (136→138), +1 M2 PV (91→92). ΔC-index < 0.005 ; model coefficients changed by < 2%. The architecture is unchanged from V230 — only the cohort is enriched by 2 real-world cases.

## Confidence intervals (Wilson)

The calculator displays **95% Wilson confidence intervals** on observed NPV/PPV for the patient's risk group, addressing calibration uncertainty at extreme risk values (where N < 10 in development cohort).

## Privacy

All computations run client-side in the browser. No patient data is transmitted to any server. The model coefficients are loaded at startup and remain local.

## Citation

If you use this calculator, methodology, or any derived artifact (code, JSON parameters, panel definitions), please cite it via the [CITATION.cff](./CITATION.cff) file (GitHub displays a "Cite this repository" button in the right sidebar that exports BibTeX / APA / etc.).

**Bibtex (placeholder — replace DOI once Zenodo archive is set up)** :

```bibtex
@software{cxj1_mrd_calculator_v230_4,
  author       = {Caulier, Alexis and {hmn-immuno-bm team}},
  title        = {{MRD ctDNA Calculator — Adult Lymphomas (CXJ1, V230.4)}},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {V230.4},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/},
  note         = {Methodology: https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/methodology.html}
}
```

**Plain-text** :

> Caulier A, hmn-immuno-bm team. MRD ctDNA Calculator — Adult Lymphomas (CXJ1, V230.4): Cox proportional hazards models for ctDNA MRD in adult B-cell and Hodgkin lymphomas after the mid-treatment timepoint. Partial pooling across IDES (MOABI hybrid-capture) and PV (phased-variant UMI) pipelines, mono-exponential decay trajectories per histology × responder strata, and integrated baseline tumor burden. Endpoint EFS (relapse OR progression). Laboratoire d'immunologie biologique GHU Mondor — secteur biologie moléculaire (hmn-immuno-bm), AP-HP, Créteil, France ; 2026. Available from: https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/. DOI: 10.5281/zenodo.XXXXXXX.

## Data availability (FAIR)

This project follows the **FAIR principles** (Findable, Accessible, Interoperable, Reusable):

- **Findable** : DOI (Zenodo, pending), GitHub repository indexed by Google Scholar / OpenAlex / ORCID claims, semantic title and abstract.
- **Accessible** : code MIT-licensed, web calculator publicly hosted (HTTPS), CITATION.cff machine-readable, methodology HTML5 standards-compliant.
- **Interoperable** : model parameters in JSON (`_calculator_data_v2.json`, `_calculator_data_v2_VAF_variants.json`) with documented schema in methodology §XII.1 and §XII.20.
- **Reusable** : MIT license, full Python pipeline (135+ scripts `_v15_*.py`) under [alessiocg/cxj1-mrd-lymphoma-pipeline](https://github.com/alessiocg/cxj1-mrd-lymphoma-pipeline), TRIPOD+AI checklist provided, version-tagged releases on GitHub.

**Patient-level data**: not shareable due to GDPR (single-center retrospective cohort of identifiable lymphoma patients). **Aggregated, de-identified cohort arrays** (predicted r12, observed T/E, predicted LP) embedded in the production JSONs enable external reproduction of all reported C-index, calibration, NRI/IDI, and DCA computations. **External validation kits** (anonymized cohort arrays + reference scripts) available upon reasonable request to hmn-immuno-bm.

## License

MIT

## Disclaimer

The calculator is a **research prototype**. It is not validated for individual clinical care. Any therapeutic decision must rely on the independent evaluation of a qualified hematologist. Calibration at extreme risk values (r12 > 60% or ctDNA_diag < 500 hEq/mL) is limited by development cohort size. **External validation on an independent cohort is in progress.**
