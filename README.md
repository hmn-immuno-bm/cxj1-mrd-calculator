# MRD ctDNA Calculator — Adult Lymphomas (CXJ1)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version V232](https://img.shields.io/badge/version-V232-16a34a.svg)](https://github.com/hmn-immuno-bm/cxj1-mrd-calculator/releases)
[![DOI (Zenodo)](https://img.shields.io/badge/DOI-pending%20(Zenodo)-lightgrey.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Cite](https://img.shields.io/badge/cite-CITATION.cff-blueviolet.svg)](./CITATION.cff)

Web calculator implementing a **Cox proportional hazards model** for prediction of **Event-Free Survival (EFS = relapse OR progression)** at 12 and 24 months after the mid-treatment timepoint (CXJ1), using **circulating tumor DNA (ctDNA) minimal residual disease** in adult B-cell lymphomas and classical Hodgkin.

🔬 **Live calculator** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/>
📚 **Methodology** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/methodology.html>
🛠️ **Pipeline source** : <https://github.com/alessiocg/cxj1-mrd-lymphoma-pipeline>
📜 **TRIPOD+AI checklist** : [TRIPOD_AI_checklist.md](./TRIPOD_AI_checklist.md)
📑 **Citation file** : [CITATION.cff](./CITATION.cff)

---

## 1. What the calculator does

From a mid-treatment plasma sample (C2J1 or C3J1) and standard clinical data, it predicts the probability of relapse or progression at 12 and 24 months. The model combines, in a single linear predictor:

1. **A time-dependent ctDNA decay score** (`F15_strat`) — log of the observed CXJ1/diag ctDNA ratio, corrected for elapsed time with a histology-specific apparent half-life.
2. **Baseline tumor burden** — log of the ctDNA load at diagnosis (or VAF if cfDNA quantification unavailable).
3. **Driver gene signal at CXJ1** (`v5A_gated`) — residual NGS signal on a histology-specific 16-gene driver panel, gated by a pipeline-specific quality filter.
4. **Interim FDG-PET** — continuous Deauville (1–5) at end of cycle 2, with ΔSUVmax% fallback if Deauville unavailable.

A patient routing layer switches automatically between **6 variants** depending on which inputs the site has available (with/without TEP × with/without Qubit cfDNA × Deauville-or-Delta-or-absent).

**Endpoint** : EFS (relapse OR progression). Non-lymphoma deaths (3 patients: 2 COVID, 1 hemorrhagic ulcer on Jehovah's Witness) are censored. Detail in methodology §XII.10.

---

## 2. The model — V232 (May 2026)

### 2.1 F15_strat (time-dependent ctDNA decay)

```
F15_strat = log(ratio_eff) + (ln(2) / t½_histo) × delai_real

ratio_eff = ratio_obs           if MRD+ (quality gate passes)
          = rc_floor             if MRD−  (polish v56 active)
ratio_obs = ctDNA_CXJ1 / ctDNA_diag    (hEq mode)
          = VAF_CXJ1   / VAF_diag      (VAF mode, no cfDNA quantification)

t½_HODGKIN    = 14 days      (justified by Tobit MLE, p < 10⁻⁹⁹)
t½_NotHodgkin = 28 days      (median KM time-to-MRD-neg 29d vs 42d)
delai_real    = days between C1J1 and CXJ1 sampling
rc_floor      = patient-specific detection threshold (1/depth for IDES, 3/PCU for PV)
```

**Property** : if a patient follows exactly an exponential decay with half-life `t½_histo`, then `F15_strat = log(ratio_initial)` — invariant in time. Slower-than-expected decay → `F15_strat` rises (bad responder). The fitted Cox coefficient `β_F15 = +0.156` confirms the expected sign (higher F15 → higher hazard).

**Portability**: by construction, this formula is robust to changes in sampling timing (J7 vs J21) and to changes in sequencing depth (the `rc_floor` is exposed as an optional input in the calculator UI). No re-fit needed for external cohorts with different protocols, only re-validation of the calibration.

### 2.2 Linear predictor and risk

```
LP_raw = β_F15 · F15_strat
       + β_burden · log_burden
       + β_v5A_pipe · v5A_gated
       + β_TEP · TEP_var          (M2 variants only)
LP     = LP_raw − LP_offset[histo]      ← centering on cohort mean, by stratum
S(t)   = S₀[pipe×histo](t)^exp(LP)
r(t)   = 1 − S(t)
```

The `LP_offset` and the baseline survival `S₀(t)` are **stratified by `pipeline × histology`** (4 strata per variant). This absorbs the differential baseline hazard between Hodgkin and DLBCL on each pipeline without spending degrees of freedom on extra covariates.

### 2.3 Partial pooling B′

`β_F15`, `β_burden`, `β_TEP` are **shared** between IDES and PV (LR test against fully separated model : `p > 0.97`). `β_v5A` is **pipeline-specific** because of the very different scales (raw reads for IDES, UMIs for PV ; LR test rejects pooling `p < 0.01`). This gives a single set of "biological" coefficients estimated jointly on the stacked cohort, while preserving the technical specificity of each pipeline. Patients sequenced on both pipelines (92 patients, doubled rows) are handled by a `cluster_col = patient_id` robust sandwich estimator.

### 2.4 Production coefficients (M2_hEq Deauv ⭐)

Shared (identical for IDES and PV pipelines):
- `β_F15_strat` = **+0.156** (p = 3.9×10⁻³)
- `β_log_burden_hEq` = **+0.231** (p = 5.5×10⁻³)
- `β_Deauville_num` = **+0.336** (p = 1.7×10⁻³) — HR per unit Deauville = 1.40; HR Deauville 5 vs 1 = 3.83

Pipeline-specific:
- `β_v5A` IDES = **+0.599** (p = 1.6×10⁻⁷), PV = **+0.682** (p = 4.5×10⁻⁶)
- `LP_offset` HODGKIN: IDES +2.26 / PV +2.07; DLBCL: IDES +3.14 / PV +3.04
- `S₀(12)` HODGKIN: 0.932 (both pipelines, IDES.HODGKIN fallback fix); DLBCL: IDES 0.878 / PV 0.866

---

## 3. The 6 variants — auto-routing by data availability

| Variant | Quantif. | TEP | N stacked (events) | C-index corr. ¹ | Notes |
|---|---|---|---|---|---|
| **M2_hEq Deauv** ⭐ | hEq (VAF×cfDNA) | Deauville 1–5 | 230 (42) | **0.852** | optimum: cascade ✓/✓ ALL_SIG ✓ |
| M1_hEq | hEq (VAF×cfDNA) | — | 273 (53) | 0.815 | fallback no PET |
| M2_VAF Deauv | VAF % | Deauville 1–5 | 230 (42) | 0.843 | fallback no Qubit cfDNA |
| M2_hEq Delta | hEq (VAF×cfDNA) | ΔSUVmax% | 200 (40) | 0.827 | fallback no Deauville (uses SUVmax shift) |
| M2_VAF Delta | VAF % | ΔSUVmax% | 200 (40) | 0.809 | dual fallback |
| M1_VAF | VAF % | — | 273 (53) | 0.793 | minimal model |

¹ Bootstrap optimism-corrected B=500 (cluster bootstrap by patient), on the stacked IDES+PV cohort.

⭐ = optimal configuration when all data is available.

---

## 4. Performance — variant M2_hEq Deauv ⭐ (V232 May 2026)

Stacked IDES+PV cohort, N=230, 42 events.

| Metric | Apparent | Optimism-corrected (B=500) |
|---|---|---|
| C-index global | 0.862 | **0.852** |
| Slope Hodgkin (target 1.0) | +1.43 | **+1.21 [+0.53 ; +1.86]** ✓ |
| Slope Non-Hodgkin (target 1.0) | +1.20 | **+1.14 [+0.98 ; +1.48]** ✓ |
| Optimism (Harrell) | +0.009 | (negligible — strong generalization) |

**Cascade validation** : both slopes' 95% CI contain 1.0 AND lower bounds > 0 → calibration validated on both histology subgroups simultaneously. All 5 covariates significant (`ALL_SIG ✓`).

**Independence from clinical scores** (V232 inherits the V230.1 result, architecture-stable) : on the 6 variants, **IPI** (DLBCL) and **Hasenclever** (Hodgkin) are absorbed by the model's LP (all `p > 0.09` in DLBCL, `p > 0.56` in Hodgkin). Conversely, LP adds significant prognostic information beyond the clinical scores (DLBCL `p_LP+ < 0.001` on 6/6 variants; Hodgkin `p_LP+ ≤ 0.012` on 6/6).

---

## 5. Pipelines

- **IDES** — MOABI hybrid-capture panel sequencing with Watch List approach; quality gate = Monte-Carlo p-value ≤ 1/100001 (10⁵ permutations).
- **PV** — Phased-variant UMI-based sequencing; quality gate = ≥ 3 positive doublets AND ≥ 3 total UMIs at CXJ1.

The PV cohort is **strictly nested in the IDES cohort** (PV ⊂ IDES, 92/138 in M2, 109/164 in M1). When both pipelines are available, PV is treated as an internal sensitivity analysis on a subset filtered by stricter PV-side QC. The two pipelines reach essentially the same C-index on the overlap (`Δ ≈ 0`), but PV is better calibrated (slope 1.09 vs 1.43 on V230.5b overlap data — V232 update pending).

---

## 6. Clinical risk stratification (Tern)

| Zone | r12 range | Suggested action |
|---|---|---|
| 🟢 **Low** | r12 ≤ cutoff_LO | Standard surveillance |
| 🟡 **Intermediate** | between LO and HI | Closer surveillance (additional PET, regular MRD) |
| 🔴 **High** | r12 > cutoff_HI | Mandatory MDT discussion (consolidation, CAR-T, trial) |

Cutoffs are **histology-adapted in the live calculator UI** :
- **Non-Hodgkin** : tricolor Tern 5 % / 45 % (3 zones)
- **Hodgkin classical** : binary cutoff 15 % (no high-risk zone — PPV plateaus at ~25-50 % even at high r12 with 9 events only)

Justification : the V230.x history (V230.2 slope 0.70 → V230.4 slope 2.14 → V230.5b slope 1.19 → V232 slope 1.21) shows that the Hodgkin calibration has been a moving target. The conservative binary 15 % cutoff is maintained until external validation, despite V232's now-validated bootstrap slope IC95% containing 1.0.

---

## 7. Local calibration (adaptive KNN)

Alongside the Cox prediction, the calculator displays a **non-parametric KM estimate** computed on the patient's LP-neighbors:

1. Take all neighbors with `|LP_i − LP_query| ≤ 1.0`
2. If fewer than 8, expand to the 8 nearest (sparse-tail safety)
3. If more than 40, truncate to the 40 nearest

Empty/sparse neighborhoods are flagged with a visual warning so the clinician knows the Cox model is extrapolating. The KM 95% Greenwood log-log confidence interval is reported numerically below the survival plot.

---

## 8. Privacy

All computations run **client-side** in the browser. No patient data is transmitted to any server. The model coefficients are loaded at startup and remain local.

---

## 9. Version history (concise)

- **V232 (May 28, 2026 — current)** ⭐ — `score_kin` (4 fitted mono-exp trajectories per pipeline) replaced by `F15_strat` (analytical formula, 2 biological constants). Same 5 covariates, same partial pooling B′, same strata pipe×histo, same v5A panels (V230.5b). **Time-dependent by construction**, portable to external cohorts with different timing or depth without re-fit. C-index optimism-corrected = 0.852 (equivalent to V230.5b at 0.852 — no discrimination cost). Cascade slope_H/slope_NH ✓/✓ validated bootstrap B=500.
- **V230.5b (May 2026 — backup)** — v5A per-GENE counting + C2-priority bi-tp fusion + Ig exclusion + Panel H swap (LTB ↔ ZCCHC7-GRHPR). Slope Hodgkin 2.14 → 1.19. Backup files: `_calculator_data_v2_V230_5b_backup.json`, `index_V230_5b_backup.html`.
- Earlier versions (V229.3, V230.1-V230.4) — see methodology §XII for full chronology of methodological iterations.

---

## 10. Citation

If you use this calculator, methodology, or any derived artifact (code, JSON parameters, panel definitions), please cite it via the [CITATION.cff](./CITATION.cff) file (GitHub displays a "Cite this repository" button in the right sidebar that exports BibTeX / APA / etc.).

**BibTeX** (placeholder — replace DOI once Zenodo archive is set up):

```bibtex
@software{cxj1_mrd_calculator_v232,
  author       = {Caulier, Alexis and {hmn-immuno-bm team}},
  title        = {{MRD ctDNA Calculator — Adult Lymphomas (CXJ1, V232)}},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {V232},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/},
  note         = {Methodology: https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/methodology.html}
}
```

**Plain text**:

> Caulier A, hmn-immuno-bm team. *MRD ctDNA Calculator — Adult Lymphomas (CXJ1, V232)*: Cox proportional hazards model for ctDNA MRD in adult B-cell and Hodgkin lymphomas after the mid-treatment timepoint. Time-dependent decay score `F15_strat` with histology-stratified apparent half-life (14d Hodgkin, 28d Non-Hodgkin), histology-specific v5A driver panels (16 H + 16 NH genes, per-gene counting), partial pooling B′ across IDES (MOABI hybrid-capture) and PV (phased-variant UMI) pipelines, strata baseline by `(pipeline × histology)`. Endpoint EFS (relapse OR progression). Laboratoire d'immunologie biologique GHU Mondor — secteur biologie moléculaire (hmn-immuno-bm), AP-HP, Créteil, France; 2026. Available from: https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/. DOI: 10.5281/zenodo.XXXXXXX.

---

## 11. Data availability (FAIR)

This project follows the **FAIR principles**:

- **Findable** — DOI (Zenodo, pending), GitHub repository indexed by Google Scholar / OpenAlex.
- **Accessible** — code MIT-licensed, web calculator publicly hosted (HTTPS), CITATION.cff machine-readable, methodology HTML5-compliant.
- **Interoperable** — model parameters in JSON (`_calculator_data_v2.json`, `_calculator_data_v2_VAF_variants.json`) with documented schema in methodology §IX and §XII.27.
- **Reusable** — MIT license; full Python pipeline (180+ scripts `_v15_*.py`) under [alessiocg/cxj1-mrd-lymphoma-pipeline](https://github.com/alessiocg/cxj1-mrd-lymphoma-pipeline); TRIPOD+AI checklist provided; version-tagged releases on GitHub.

**Patient-level data**: not shareable due to GDPR (single-center retrospective cohort of identifiable lymphoma patients). **Aggregated, de-identified cohort arrays** (predicted r12, observed T/E, predicted LP, isHodgkin) are embedded in the production JSONs (`cohort_arrays` field) and enable external reproduction of all reported C-index, calibration, NRI/IDI, and DCA computations. **External validation kits** (anonymized cohort arrays + reference scripts) available upon reasonable request to hmn-immuno-bm.

---

## 12. License & disclaimer

MIT.

The calculator is a **research prototype**. It is not validated for individual clinical care. Any therapeutic decision must rely on the independent evaluation of a qualified hematologist. Calibration at extreme risk values (`r12 > 60%` or `ctDNA_diag < 500 hEq/mL`) is limited by development cohort size. **External validation on an independent cohort is in progress.**
