# MRD ctDNA Calculator — Adult Lymphomas (CXJ1)

Web calculator implementing **Cox proportional hazards models** for prediction of **Event-Free Survival (EFS = relapse OR progression)** at 12 and 24 months after the mid-treatment timepoint (CXJ1), using circulating tumor DNA (ctDNA) minimal residual disease (MRD) markers in adult B-cell lymphomas and classical Hodgkin.

🔬 **Live calculator** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/>
📚 **Methodology** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/methodology.html>

## Overview

The calculator predicts post-CXJ1 EFS risk from four covariates measured around the CXJ1 timepoint (typically end of cycle 2 or 3 of first-line chemo-immunotherapy):

1. **Kinetic score** — Gaussian log-likelihood distance between observed and predicted ctDNA decay ratios under "good responder" vs "bad responder" mono-exponential trajectories (fitted per histology × responder strata).
2. **Baseline tumor burden** — log(1 + ctDNA at diagnosis) where ctDNA = VAF × cfDNA in hEq mode, or log(1 + VAF at diagnosis) in VAF mode.
3. **Driver gene signal (v5A_gated)** — residual NGS signal on the 11 oncogenic drivers (D11: BCL2, BCL6, BCL7A, BTG2, CIITA, CXCR4, IRF8, MYC, PAX5, S1PR2, TP53) at CXJ1, gated by a pipeline-specific quality filter.
4. **Interim FDG-PET** — Deauville score (linear, 1-5) at end of cycle 2, with **ΔSUVmax% fallback** if Deauville unavailable.

**Endpoint**: EFS (Event-Free Survival = relapse OR progression). Non-lymphoma deaths are censored (3 patients : 2 COVID, 1 hemorrhagic ulcer under transfusion contraindication). Justification in methodology §XI.10.

## 6 model variants (auto-selected from inputs)

The calculator switches automatically between 6 variants based on data availability:

| Variant | Quantification | TEP | IDES N (events) | C-IDES (corr.) | PV N (events) | C-PV (corr.) |
|---|---|---|---|---|---|---|
| **M2_hEq Deauv** ⭐ | hEq (VAF × cfDNA) | Deauville 1-5 | **138 (26)** | **0.84** | **92 (16)** | **0.87** |
| M2_hEq Delta | hEq (VAF × cfDNA) | ΔSUVmax% (fallback) | 120 (25) | 0.82 | 80 (15) | 0.88 |
| M1_hEq | hEq (VAF × cfDNA) | — absent | 164 (32) | 0.81 | 109 (21) | 0.85 |
| M2_VAF Deauv | VAF % only | Deauville 1-5 | 138 (26) | 0.81 | 92 (16) | 0.81 |
| M2_VAF Delta | VAF % only | ΔSUVmax% (fallback) | 120 (25) | 0.81 | 80 (15) | 0.83 |
| M1_VAF | VAF % only | — absent | 164 (32) | 0.77 | 109 (21) | 0.79 |

⭐ = optimal configuration when all data is available. C-index optimism-corrected (bootstrap V230.1 B=500, Harrell).

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

For each variant, the Cox model is fitted with **partial pooling B′** stratified by pipeline (`cluster=NOM`, L2 penalizer = 0.05):

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

## Clinical risk stratification (Tern 5/35, V230)

| Zone | r12 range | Suggested action |
|---|---|---|
| 🟢 **Low** | r12 ≤ 5% | Standard surveillance |
| 🟡 **Intermediate** | 5% < r12 ≤ 35% | Closer surveillance (additional PET, regular MRD) |
| 🔴 **High** | r12 > 35% | Mandatory MDT discussion (consolidation, CAR-T, trial) |

Thresholds **5/35** selected by exhaustive scan of 117 candidate (lo, hi) couples on the V230 EFS cohort (V15.283/284). 5/35 is Pareto-optimal: it maximizes the χ² separation Low ↔ Intermediate on both IDES (χ²=4.10, p=0.043) and PV (χ²=5.00, p=0.025) while keeping the High zone strongly marked (N=20 IDES with EFS@24=21%, N=9 PV with EFS@24=0%).

## Local calibration (adaptive KNN)

Alongside the Cox prediction, the calculator displays a **non-parametric KM estimate** computed on the patient's LP-neighbors:

1. Take all neighbors with `|LP_i − LP_query| ≤ 1.0`
2. If fewer than 8, expand to the 8 nearest (sparse-tail safety)
3. If more than 40, truncate to the 40 nearest

Validated empirically (V208 leave-one-out): ACE_tail divided by ~10× vs a fixed K=30. Empty/sparse neighborhoods are flagged with a visual warning so the clinician knows the Cox model is extrapolating.

The KM 95% Greenwood log-log confidence interval is reported numerically below the survival plot.

## Validation (V230.1, May 2026)

- **Bootstrap optimism-corrected C-index** (Harrell B=500): M2_hEq Deauv IDES ⭐ apparent 0.845 → corrected **0.835** (optimism = +0.010, IC95% percentile [0.76 – 0.92]) ; PV apparent 0.893 → corrected **0.873** (optimism = +0.020, IC95% [0.80 – 0.97]). Max optimism across 12 fits = +0.027, median +0.015.
- **Linearity**: martingale residual lowess amplitude = 0.13 (< 0.3 threshold); Grambsch-Therneau PH test p > 0.18 for all covariates; LR test for adding RCS splines on log_burden p = 0.80 NS.
- **β_log_burden stability**: 95% percentile CI [0.10, 0.47], 100% of bootstraps positive.
- **Partial pooling B′ validation**: LR test of sharing β_score + β_burden + β_TEP between IDES and PV: χ² ≈ 0, df = 3, p = 1.00 ✓. Sharing v5A rejected (LR p=0.041).
- **Predictive contribution per covariate** (LR test, M2_hEq Deauv): log_burden p = 0.002 ; v5A_gated p < 0.001 ; Deauville_num p < 0.001 ; score_kin p ≈ 0.005.
- **Independence vs IPI / Hasenclever** (V230.1, V15.296): on the 6 variants, IPI and Hasenclever are **absorbed by LP V230.1** (all p_score+ > 0.09 in DLBCL, > 0.56 in Hodgkin). Reciprocally, LP V230.1 adds significant information beyond clinical scores (DLBCL p_LP+ < 0.001 on 6/6 variants; Hodgkin p_LP+ ≤ 0.012 on 6/6).

## V230.1 patches (technical reintegration of 2 patients)

V230.1 = V230 + 2 patients re-integrated after audit (V15.289-291):

- **PATIENT_A** (DLBCL bad responder): inconsistent MOABI WL prefix (`NGA_*.xlsx` at C2J1 vs `YOU_*.xlsx` at C3J1, due to compound surname). The default parser only matched `YOU` and silently lost the C2J1 file. Patch: `scripts/_patient_aliases.py` documents the `NGA → PATIENT_A` alias (NAS files untouched). Audit V15.289 confirmed this was the only such case among 191 WL_NGS files.
- **PATIENT_B** (Hodgkin bad responder): C1J1 date error in the clinical database (`2022-11-03` was actually C3J1+1d due to a copy-paste in the wrong cell ; the true C1J1 is `2022-09-21`, confirmed by user 2026-05-20). With the corrected date, `delai_real(C2J1) = +21 d` and the patient correctly enters the M1 IDES cohort. Patch: `scripts/_clinical_patches.py` applies the correction in memory (NAS clinical file untouched).

**Impact**: +2 M1 IDES (162→164), +1 M1 PV (108→109), +2 M2 IDES (136→138), +1 M2 PV (91→92). ΔC-index < 0.005 ; model coefficients changed by < 2%. The architecture is unchanged from V230 — only the cohort is enriched by 2 real-world cases.

## Confidence intervals (Wilson)

The calculator displays **95% Wilson confidence intervals** on observed NPV/PPV for the patient's risk group, addressing calibration uncertainty at extreme risk values (where N < 10 in development cohort).

## Privacy

All computations run client-side in the browser. No patient data is transmitted to any server. The model coefficients are loaded at startup and remain local.

## Citation

> Calculateur MRD CXJ1 — Cox proportional hazards models for ctDNA MRD in adult B-cell and Hodgkin lymphomas after the mid-treatment timepoint. Partial pooling across IDES (MOABI hybrid-capture) and PV (phased-variant UMI) pipelines, with mono-exponential decay trajectories per histology × responder strata and integrated baseline tumor burden. Endpoint EFS (relapse OR progression). Laboratoire d'immunologie biologique GHU Mondor — secteur biologie moléculaire (hmn-immuno-bm), 2026.

## License

MIT

## Disclaimer

The calculator is a **research prototype**. It is not validated for individual clinical care. Any therapeutic decision must rely on the independent evaluation of a qualified hematologist. Calibration at extreme risk values (r12 > 60% or ctDNA_diag < 500 hEq/mL) is limited by development cohort size. **External validation on an independent cohort is in progress.**
