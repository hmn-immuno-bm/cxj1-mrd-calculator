# MRD ctDNA Calculator — Adult Lymphomas (CXJ1) — v2 NEW

Web calculator implementing **Cox proportional hazards models** (V15.185, M2 v2)
for prediction of progression-free survival (PFS) at 12 and 24 months
after CXJ1 mid-treatment timepoint, using circulating tumor DNA (ctDNA)
minimal residual disease (MRD) markers.

🔬 **Live calculator (v2)** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/index_v2.html>
🔬 **Live calculator (v1, production)** : <https://hmn-immuno-bm.github.io/cxj1-mrd-calculator/>

## What's new in v2

| Aspect | v1 (V15.130, production) | **v2 (V15.185, NEW)** |
|---|---|---|
| **Trajectory model** | Bi-exponential (3 params) | **Mono-exponential (1 param)** |
| **VAF representation** | VAF % (relative) | **Absolute ctDNA quantity (hEq) = VAF × cfDNA** |
| **Mid-treatment timing** | Theoretical (J21/J42) | **Real Glims sampling date** (fallback to theoretical for 18 aberrant cases) |
| **Baseline burden** | Not used | **Integrated into score** : score = score_kinetic + log(1 + ctDNA_diag) |
| **Cox covariates** | 3 (score + v5A + TEP) | **3 (same)** |
| **Number of trajectory params** | 12 (3 × 4 trajectories) | **4 (1 × 4 trajectories)** |

## Two Cox models (toggle inside calculator)

### M2 v2 — with TEP (3 variables)

| Pipeline | N | Events | C-index | AIC |
|---|---|---|---|---|
| IDES (MOABI panel) | 135 | 24 | **0.859** | 183.36 |
| PV (phased variants) | 89 | 15 | **0.928** | 87.56 |

### M1 v2 — without TEP (2 variables, fallback)

| Pipeline | N | Events | C-index | AIC |
|---|---|---|---|---|
| IDES | 159 | 30 | **0.811** | 256.78 |
| PV | 104 | 20 | **0.847** | 143.95 |

→ Significant improvement vs v1 M1 (IDES: 0.752 → **0.811**, +0.06).

## Score formula (M2 v2)

### Step 1 — inputs
- `VAF_diag` (%) — mean VAF at diagnosis
- `VAF_CXJ1` (%) — mean VAF at mid-treatment, **polished** (set to 0 if quality gate fails)
- `cfDNA_diag` (hEq/mL) — plasma cell-free DNA concentration at diagnosis
- `cfDNA_CXJ1` (hEq/mL) — at mid-treatment
- `t` (days) — real interval from C1J1 to CXJ1 sampling
- `histo` — Hodgkin classical OR DLBCL/other B-cell

### Step 2 — ctDNA quantities (haploid genome equivalents)

```
ctDNA_diag = VAF_diag × cfDNA_diag
ctDNA_CXJ1 = VAF_CXJ1_polished × cfDNA_CXJ1
```

### Step 3 — Mono-exponential trajectories (4 parameters total !)

```
mono_exp(t, λ) = exp(−λ·t)
```

Fitted parameters per (histology × responder class) :

| Histology | λ_good (j⁻¹) | t½ good | λ_bad (j⁻¹) | t½ bad | σ² |
|---|---|---|---|---|---|
| **IDES Hodgkin** | 0.651 | 1.1 d | 0.371 | 1.9 d | 30.48 |
| **IDES DLBCL** | 0.388 | 1.8 d | 0.217 | 3.2 d | 33.23 |
| **PV Hodgkin** | 0.565 | 1.2 d | 0.241 (DLBCL fallback) | — | 29.51 |
| **PV DLBCL** | 0.372 | 1.9 d | 0.241 | 2.9 d | 25.22 |

### Step 4 — Kinetic score (log-likelihood ratio)

```
ratio_obs = ctDNA_CXJ1 / ctDNA_diag   (clamp 1e-7 if ctDNA_CXJ1 = 0)
pred_good = exp(−λ_good · t)
pred_bad  = exp(−λ_bad  · t)

log L_good = −(log(ratio_obs) − log(pred_good))² / (2σ²)
log L_bad  = −(log(ratio_obs) − log(pred_bad))²  / (2σ²)

score_kinetic = log L_bad − log L_good
```

### Step 5 — Total score (with baseline burden, coefficient 1)

```
score_NEW = score_kinetic + log(1 + ctDNA_diag)
```

This is a key innovation: the **diagnostic ctDNA burden is integrated directly into the score** (no separate covariate), reflecting that high initial tumor load is prognostic independently of the relative response dynamics.

### Step 6 — Cox prediction

```
h(t | x) = h₀(t) · exp(β₁·score_NEW + β₂·v5A_gated + β₃·TEP2_pos)
```

| Covariate | coef | HR | p (IDES M2) |
|---|---|---|---|
| score_NEW | 0.335 | 1.40 | 0.001 |
| v5A_gated | 0.592 | 1.81 | <10⁻⁶ |
| TEP2_pos | 0.952 | 2.59 | 0.010 |

## Risk stratification (Tern 10/40)

| Zone | r12 range | IDES PFS@12 / @24 | PV PFS@12 / @24 |
|---|---|---|---|
| 🟢 **Low** | r12 ≤ 10% | 95.1% / 95.1% | 98.3% / 98.3% |
| 🟡 **Intermediate** | 10% < r12 ≤ 40% | 85.5% / 82.1% | 83.1% / 68.2% |
| 🔴 **High** | r12 > 40% | 18.8% / **0.0%** | 0.0% / 0.0% |

Note: the **High zone PFS@24 is 0% in both pipelines** — these patients all relapse within 24 months. The Intermediate zone shows continued late relapses (vs Low which is stable) → close surveillance justified.

## Methodological validation history

The v2 model was selected after **systematic comparison of 11+ alternative formulations** including :
- Pooled VAF (Σ Reads/Σ Depth) with multiple imputation strategies (V155-162) — all degraded
- Different VAF aggregations (POS only, with zeros, IDES vs MOABI source) (V163-173)
- Bi-temporal aggregation strategies (V164: MIN, MEAN, MAX, C3-priority) — all degraded
- Longitudinal joint models (V166: slope, NLME f_patient, time-varying Cox) — all NS
- Real dates (delai_real) — equivalent to delai_corrige (V165)
- Multiple Tern thresholds with cross-validation (V176) — 10/40 confirmed optimal
- Polish strengthening with K-thresholds (V171) and soft weighting (V172) — all degraded
- LCMM stratified trajectories with baseline ctDNA (V181) — overfit
- Histo-specific α coefficients (V180) — power-limited
- ctDNA score formulations (V178) — V184 mono-exp + log_burden optimal

→ The v2 model represents the **best compromise** of model parsimony (3× fewer parameters) and prognostic accuracy (+0.013 C-index) validated on the development cohort.

## Confidence intervals (Wilson)

The v2 calculator displays **95% Wilson confidence intervals** on observed VPP/VPN for the patient's risk group, addressing calibration uncertainty at extreme risk values (where N < 10 in development cohort).

## Citation

> Calculateur MRD CXJ1 v2 (V15.185) — Cox proportional hazards models for ctDNA MRD
> in adult lymphomas after CXJ1 mid-treatment. Mono-exponential trajectories on
> absolute ctDNA quantity (hEq) with integrated baseline burden term. Cohort
> development : MOABI panel (IDES) and phased variants (PV). Laboratoire d'immunologie
> biologique GHU Mondor — secteur biologie moléculaire (hmn-immuno-bm), 2026.

## License

MIT
