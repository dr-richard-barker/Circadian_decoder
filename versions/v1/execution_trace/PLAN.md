# Plan: Circadian Time-of-Day Decoding of Plant Spaceflight Transcriptomes

## Summary

Build a reproducible analysis pipeline that applies the ChronoGauge circadian-time (CT) ensemble model — the successor to the time-of-day decoder introduced in Gardiner et al. 2021 PNAS (PMC8364196) — to NASA GeneLab Arabidopsis spaceflight RNA-seq data, to test whether spaceflight alters plant circadian phase. The analysis runs as a **forked design**: (A) the 15-accession set from Barker et al. 2023 npj Microgravity, and (B) all available Arabidopsis spaceflight RNA-seq studies in GeneLab, so we can compare how dataset scope affects conclusions. Statistical testing uses per-study flight-vs-ground comparison of CT residuals, followed by random-effects meta-analysis across studies, stratified by tissue, genotype, and hardware. The deliverable is a full submission-ready package: all figures/tables, npj Microgravity-format manuscript, LaTeX, GitHub repo, and Zenodo-ready archive.

---

## Key decisions (locked)

| Decision | Choice |
|---|---|
| Model | ChronoGauge (2025, Reynolds et al., Nature Comms) — bagging ensemble, pre-trained HuggingFace models, MIT license |
| Dataset scope | **Fork**: (A) 15 accessions from 2023 meta-analysis; (B) all Arabidopsis spaceflight RNA-seq in GeneLab |
| Phase-shift test | Per-study flight-vs-ground test, then random-effects meta-analysis across studies |
| Deliverable | Analysis + figures first → then full submission-ready manuscript/LaTeX/GitHub/Zenodo |

---

## Subsystem 1: Environment and dependencies

**Goal**: A single reproducible Python + R environment.

- Python 3.9+ environment via `uv`: tensorflow, scikit-learn, numpy, pandas, tqdm, matplotlib, seaborn, huggingface-hub, requests, pyyaml
- R environment for meta-analysis and ggplot2 figures: metafor, ggplot2, ggprism, ComplexHeatmap, svglite, dplyr, tidyr, readr
- Clone ChronoGauge repo (`ConnorReynoldsUK/ChronoGauge`) and download pre-trained RNA-seq ensemble from HuggingFace (`conjr94/ChronoGauge_RNAseq`)
- All code runs on a right-sized sandbox machine (8 CPU / 32 GB); no HPC needed — ChronoGauge inference is <1 min per sample for 100 sub-models

---

## Subsystem 2: Data acquisition (GeneLab OSDR)

**Goal**: Download processed RNA-seq + metadata for both forks.

### Fork A — 15 accessions (Barker 2023)
GLDS-7, 17, 37, 38, 44, 46, 120, 121, 136, 147, 205, 208, 213, 218, 251

### Fork B — all Arabidopsis spaceflight RNA-seq
Query OSDR API (`https://osdr.nasa.gov/bio/api/spaceflight/studies`) filtering by organism = Arabidopsis thaliana and assay type = RNA-seq. Cross-reference with the 15-accession list; Fork B is a superset.

### Per-study files to download
1. **`*Normalized_Counts.csv`** — DESeq2 median-of-ratios normalized gene counts (gene × sample matrix)
2. **ISA-Tab metadata zip** (`OSD-*_metadata_OSD-*-ISA.zip`) → parse `s_*.txt` sample table for: tissue, genotype/ecotype, spaceflight/ground, hardware, light regime, harvest age, photoperiod
3. **`*differential_expression*.csv`** — pre-computed DEG tables (for cross-referencing clock genes)

### Metadata curation (critical step)
Extract and harmonize into a single sample-level table:
- `study_id` (GLDS/OSD accession)
- `sample_name`, `organism`, `ecotype`/`genotype` (Col-0, WS, Ws-2, Ler-2, Cvi-0, phyD, met1-7, elp2-5, etc.)
- `tissue` (root, leaf, shoot, hypocotyl, whole seedling, callus)
- `condition` (flight / ground_control)
- `hardware` (Veggie, BRIC, APH, EMCS, APEX, TAGES)
- `light_regime` (continuous light / LL, photoperiod LD, dark)
- `photoperiod_hours` (e.g., 24:0, 16:8, 12:12)
- `harvest_age_days`
- `zeitgeber_time_ZT` (hours after lights-on, if recoverable from metadata; many studies harvest at a single timepoint)
- `replicate_group`

**ZT/CT handling**: ChronoGauge was trained under constant light (LL → CT labels). For GeneLab samples under LL, predicted CT is directly comparable to experimental ZT. For photoperiod (LD) samples, we treat predicted CT as an estimate of internal circadian phase and compare flight vs. ground *within the same light regime*. Studies where ZT cannot be recovered are flagged and handled as "unknown-phase" — included in the flight-vs-ground CT comparison but excluded from CT−ZT residual analysis.

---

## Subsystem 3: ChronoGauge application

**Goal**: Predict circadian time (CT, 0–24h) for every GeneLab sample.

### Gene ID mapping
- GeneLab Normalized_Counts uses Ensembl/TAIR10 gene IDs (AGI codes like `AT1G01010`)
- ChronoGauge feature sets use AGI codes — direct mapping expected; verify and log any missing features
- For each of the 100 sub-models: extract that model's feature gene set from the expression matrix, MinMaxScale (fit on ChronoGauge training data stats), predict CT
- Aggregate 100 predictions via **circular mean** → single CT estimate per sample
- Record circular variance across sub-models as an **uncertainty/confidence metric** ("circadian fingerprint")

### Validation on ground controls
- For studies with recoverable ZT: plot predicted CT vs. known ZT for ground controls → expected ~y=x line if the clock is functioning
- Report MAE (mean absolute error, circular) and R² on ground controls as model validation
- This is the key "does the model work on spaceflight data" sanity check

### Output
Per-sample table: `study_id, sample_name, tissue, genotype, hardware, condition, light_regime, ZT_known, predicted_CT, CT_circular_variance, n_submodels_used`

---

## Subsystem 4: Statistical analysis (per-study then meta)

**Goal**: Test whether spaceflight shifts circadian phase, stratified by tissue/genotype/hardware.

### Step 1 — Per-study analysis
For each study with both flight and ground samples:
- **Primary metric**: CT residual = circular_distance(predicted_CT, known_ZT) for studies with ZT; or raw predicted_CT for studies without ZT (compare flight vs. ground distributions)
- **Test**: Watson-Williams multi-sample test for circular means (flight vs. ground), or linear mixed model on CT residuals if ZT known
- **Effect size**: mean circular phase shift (flight − ground) in hours, with 95% CI
- Require ≥3 biological replicates per condition per study; flag studies with fewer

### Step 2 — Random-effects meta-analysis
- Pool per-study effect sizes via `metafor::rma()` (random-effects, REML estimator)
- Primary forest plot: phase shift (hours) per study + pooled estimate
- Heterogeneity: Cochran's Q, I², τ²
- **Stratified meta-analyses** by:
  - Tissue (root vs. leaf vs. shoot vs. whole)
  - Genotype/ecotype (Col-0 vs. WS vs. others)
  - Hardware (Veggie vs. BRIC vs. APH vs. EMCS)
  - Light regime (LL vs. LD)
- Sensitivity: leave-one-study-out, funnel plot asymmetry (Egger's test)

### Step 3 — Fork comparison
- Run the full pipeline identically on Fork A (15 accessions) and Fork B (all)
- Compare pooled effect sizes, heterogeneity, and stratified conclusions
- This directly answers: "does expanding the dataset change the conclusion about spaceflight and circadian phase?"

---

## Subsystem 5: Supplementary ML and data visualization

**Suggested additions beyond the core analysis:**

1. **Circadian fingerprint analysis** — ChronoGauge's per-sub-model variance as a "clock robustness" metric. Test if flight samples show higher sub-model disagreement (less coherent clock) than ground. ANOVA or Kruskal-Wallis by condition × tissue.

2. **Core clock gene expression** — Extract CCA1, LHY, TOC1, PRR9/7/5, GI, ELF4, ELF3, LUX, TIC expression from Normalized_Counts. Heatmap of flight/ground log2FC per study, annotated by predicted CT shift. Tests whether phase shift is driven by canonical clock genes.

3. **Circadian gene set enrichment** — Among each study's DEGs (flight vs. ground), test enrichment of the MetaCycle-defined circadian gene set (from PMC8364196, ~9,394 genes) via ORA (clusterProfiler). Are spaceflight DEGs enriched for clock-controlled genes?

4. **PCA/UMAP colored by predicted CT** — Dimensionality reduction of normalized counts, colored by predicted CT, shaped by flight/ground. Visual check for whether flight samples occupy a different region of circadian phase space.

5. **Polar/circular phase plots** — Rose diagrams of predicted CT for flight vs. ground, per tissue. Intuitive visualization of phase distribution shifts.

6. **Bayesian phase estimation** (optional, if time permits) — Wrap ChronoGauge predictions in a simple Bayesian model to get posterior CT distributions per sample, enabling uncertainty-aware phase-shift testing.

---

## Subsystem 6: Figures and tables

### Main figures (target: 5–6, npj Microgravity allows up to 8)

| # | Figure | Content |
|---|---|---|
| 1 | Dataset overview | Sankey/alluvial: studies × tissue × genotype × hardware × light regime, for both forks. Table 1 embedded or adjacent. |
| 2 | ChronoGauge validation | Predicted CT vs. known ZT for ground controls (scatter + y=x line), MAE annotation, per-tissue facet. Proves model works on GeneLab data. |
| 3 | Primary meta-analysis forest plot | Phase shift (hours) per study + pooled estimate, Fork A and Fork B side-by-side. The headline result. |
| 4 | Stratified phase shifts | Panel grid: phase shift by tissue (a), genotype (b), hardware (c), light regime (d). Box/violin or forest subplots. |
| 5 | Circadian fingerprint / clock robustness | Sub-model circular variance, flight vs. ground, per tissue. Tests clock coherence. |
| 6 | Core clock gene response | Heatmap of log2FC (flight/ground) for ~15 core clock genes × studies, annotated by phase-shift direction. |

### Supplementary figures
- S1: Per-study predicted CT distributions (flight vs. ground histograms/polar plots)
- S2: PCA/UMAP colored by predicted CT
- S3: Circadian gene set enrichment results per study
- S4: Fork A vs. Fork B comparison summary (all stratified results side-by-side)
- S5: Sensitivity analyses (leave-one-out, funnel plots)
- S6: ChronoGauge feature gene coverage per study

### Tables
- **Table 1**: Dataset characteristics (study, OSD/GLDS ID, tissue, genotype, hardware, light regime, n flight, n ground, ZT recoverable?)
- **Table 2**: Per-study phase-shift results (effect size, 95% CI, p-value, n)
- **Table 3**: Meta-analysis pooled results by stratum (tissue/genotype/hardware/light)
- **Table S1**: Full sample metadata table
- **Table S2**: ChronoGauge prediction outputs per sample

### Format
- All figures: SVG (primary) + PNG (300 DPI), colorblind-friendly palettes, Liberation Sans font
- R/ggplot2 for forest plots, heatmaps (ComplexHeatmap), polar plots; Python/seaborn for PCA/UMAP and validation scatter
- Media output check on every figure before finalizing

---

## Subsystem 7: Manuscript (npj Microgravity format)

**Structure** (Nature Portfolio format):

1. **Abstract** (~200 words, unstructured)
2. **Introduction** — circadian clock in plants, spaceflight stress, why time-of-day decoding matters, ChronoGauge as successor to PMC8364196 approach, study aims
3. **Results** — dataset overview → model validation → primary meta-analysis → stratified results → fork comparison → clock gene/fingerprint analyses
4. **Discussion** — interpretation by tissue/genotype/hardware, comparison to prior spaceflight circadian literature, limitations (single-timepoint, light regime confounds, cross-ecotype generalization), implications for space agriculture
5. **Methods** — data acquisition, ChronoGauge application, statistical framework, fork design, software/versions
6. **Data availability** — GeneLab accessions, ChronoGauge repo/HuggingFace, this study's GitHub/Zenodo
7. **Code availability** — GitHub repo URL, Zenodo DOI
8. **References** — inline numbered, Vancouver style

**LaTeX document**: `manuscript.tex` using a Nature-compatible template (or `elsarticle`/`sn-jnl` class), pulling in all figures via `\includegraphics`, tables via `booktabs`, with separate figure legends file. Compiled to PDF.

---

## Subsystem 8: GitHub and Zenodo packaging

### GitHub repo structure
```
spaceflight-circadian-decoder/
├── README.md                    # Project overview, quickstart, citation
├── LICENSE                       # MIT (ChronoGauge) + CC-BY-4.0 (our code)
├── Dockerfile                    # Reproducible environment
├── environment.yml / pyproject.toml  # Python deps
├── requirements_R.txt            # R deps
├── .gitignore
├── config/
│   ├── accessions_forkA.yaml     # 15-accession list
│   ├── accessions_forkB.yaml     # all-Arabidopsis list (generated)
│   └── analysis_params.yaml      # thresholds, stratification vars
├── src/
│   ├── data_acquisition.py       # OSDR API download + ISA-Tab parsing
│   ├── metadata_curation.py      # harmonize sample metadata
│   ├── chronogauge_apply.py      # load models, predict CT, circular mean
│   ├── statistical_analysis.py   # per-study tests + metafor call
│   ├── figures/                  # one script per figure
│   └── tables/                   # one script per table
├── notebooks/
│   ├── 01_data_acquisition.ipynb
│   ├── 02_chronogauge_prediction.ipynb
│   ├── 03_statistical_analysis.ipynb
│   └── 04_figures_tables.ipynb
├── manuscript/
│   ├── manuscript.tex
│   ├── references.bib
│   ├── figures/                  # final SVG/PNG
│   └── tables/                   # final CSV/LaTeX tables
├── results/                      # generated outputs (gitignored, in Zenodo)
└── tests/
    └── test_pipeline.py          # smoke test on 1 study
```

### Zenodo archive
- Tagged GitHub release → Zenodo DOI
- Archive includes: code + results/ (all figures, tables, prediction outputs, meta-analysis objects) + manuscript PDF + supplementary
- `CITATION.cff` file for Zenodo metadata
- `zenodo_upload_manifest.json` listing all files to include

---

## Compute / resource estimate

| Component | Estimate | Basis |
|---|---|---|
| Data download (Fork A: 15 studies) | ~30 min, ~2 GB | Normalized_Counts ~50–150 MB each + metadata zips |
| Data download (Fork B: ~25–35 studies) | ~60 min, ~5 GB | ~2× Fork A |
| Metadata curation | ~1–2 h (mostly automated + manual review) | ISA-Tab parsing + ZT extraction |
| ChronoGauge inference | ~15–30 min for ~1000 samples | 100 sub-models × <1s per sample |
| Statistical analysis | ~5 min | metafor is fast |
| Figure generation | ~10 min | R + Python plotting |
| Manuscript writing | LLM generation + review | Not compute-bound |
| **Total compute** | ~2–3 h on one 8-CPU/32-GB machine | No HPC needed |

**Execution target**: Single right-sized sandbox machine (8 CPU / 32 GB) via ManageMachine. Data staged in `/workspace/`, checkpoints to `/mnt/shared-workspace/`, deliverables to `/mnt/results/`.

---

## Execution order

1. **Environment setup** — install deps, clone ChronoGauge, download HuggingFace models
2. **Data acquisition** — download Fork A then Fork B, parse metadata
3. **Metadata curation** — harmonize sample table, extract ZT/light regime
4. **ChronoGauge application** — predict CT for all samples, both forks
5. **Model validation** — predicted vs. known ZT on ground controls
6. **Statistical analysis** — per-study tests, meta-analysis, stratification, fork comparison
7. **Supplementary ML** — fingerprint, clock genes, enrichment, PCA, polar plots
8. **Figure generation** — all main + supplementary figures (with media output checks)
9. **Table generation** — all main + supplementary tables
10. **Manuscript** — write full npj Microgravity text
11. **LaTeX** — compile manuscript + figures + tables
12. **GitHub packaging** — repo structure, README, Docker, tests
13. **Zenodo packaging** — release tag, archive manifest, CITATION.cff

Each step surfaces outputs to `/mnt/results/` as they complete for user review.

---

## Assumptions

1. GeneLab Normalized_Counts.csv files use AGI gene codes (AT*G*) compatible with ChronoGauge features — will verify on first download; if Ensembl IDs, add mapping step.
2. ZT/harvest-time metadata is recoverable for a meaningful subset of studies from ISA-Tab or associated publications — if not, the CT−ZT residual analysis shrinks but the flight-vs-ground CT comparison remains valid.
3. ChronoGauge pre-trained RNA-seq ensemble on HuggingFace contains 100 saved Keras/TF models + feature gene lists — will verify file structure on download.
4. Fork B (all Arabidopsis RNA-seq) will yield ~25–35 studies — will confirm via OSDR API query.
5. npj Microgravity accepts Nature Portfolio LaTeX format — will use `sn-jnl` or equivalent class.
6. Studies with <3 replicates per condition are flagged but retained in descriptive analyses; excluded from per-study inferential tests.
