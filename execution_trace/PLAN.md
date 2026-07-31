# Plan: Figure Fixes + Three Additional ML/Visualization Approaches

## Summary

Fix readability issues across all 12 existing figures (legends overlapping data, labels clipped, broken panels), then implement three new analyses: (1) circadian trajectory analysis via t-SNE/PCA on the 100-sub-model prediction space, (2) gene set enrichment (limma-voom + fgsea) on the 4 studies with significant phase shifts, and (3) an interactive Plotly dashboard deployable to GitHub Pages. All new outputs integrate into the manuscript as supplementary figures S7–S9, new tables, updated text, and a regenerated PDF. The GitHub repo archive is updated.

## Part 1: Figure Readability Fixes

### Problems identified from code analysis

| Figure | Problem | Fix |
|--------|---------|-----|
| Fig 1 (all panels) | Legends placed inside plot area with default `loc='best'` — overlap bars in stacked bar charts, especially Panel B (10+ hardware categories) | Move legends outside plots using `bbox_to_anchor=(1.02, 1)` with `loc='upper left'`; increase figure width to accommodate |
| Fig 3 (polar) | Legend at `bbox_to_anchor=(1.3, 1.1)` — may be clipped despite `bbox_inches='tight'` | Adjust to `(1.15, 1.05)`, reduce marker size in legend |
| Fig 4 (forest) | Study labels hardcoded at `x=-13` — off-chart since x-axis auto-scales from data range | Compute label x-position from axis left limit minus padding; use `ax.annotate` with `xycoords='axes fraction'` for label column, or use a two-panel layout (text column + plot column) |
| Fig 5 (stratified) | `suptitle` at `y=1.02` with `tight_layout` — title may be clipped | Use `constrained_layout=True` or `fig.suptitle(..., y=0.98)` with `plt.subplots_adjust(top=0.9)` |
| Fig S1 Panel B | 23-study legend with `ncol=2` inside plot — overlaps scatter | Move legend outside plot to the right, use `ncol=1` with small font, or use a colorbar instead of legend for 23 categories |
| Fig S3/S4 | Long hardware/light regime labels rotated 30° — clipped at bottom | Rotate to 45°, use `ha='right'`, add `plt.subplots_adjust(bottom=0.25)` or increase figure height |
| Fig S5 | Same `x=-13` label issue as Fig 4 | Same fix as Fig 4 |
| Fig S6 Panel A | References `mean_circular_variance` column that doesn't exist in `per_study_results` — silently plots zeros | Compute mean circular variance per study from `predictions_df` and merge into per_study before plotting |

### Implementation

Rewrite the affected functions in `/workspace/spaceflight-circadian-decoder/src/figures/generate_figures.py`:
- `fig1_study_overview`: legends outside, wider figure (figsize 18×6)
- `fig3_phase_shift_polar`: legend repositioned
- `fig4_forest_plot`: label column using axes-fraction coordinates
- `fig5_stratified_analysis`: constrained_layout, suptitle y adjustment
- `figS1_pca_umap`: Panel B legend outside or colorbar
- `figS3_phase_by_hardware`: 45° rotation, bottom margin
- `figS4_phase_by_light_regime`: 45° rotation, bottom margin
- `figS5_fork_comparison`: same label fix as fig4
- `figS6_model_uncertainty`: compute mean_circular_variance from predictions_df

Re-run figure generation from cached data (predictions + per_study_results + metadata). No need to re-run ChronoGauge or statistical analysis.

## Part 2: Circadian Trajectory Analysis (t-SNE + PCA)

### Approach
Use the 100 ChronoGauge sub-model predictions per sample (`all_model_predictions.csv`, 603 samples × ~100 features) as a circadian fingerprint feature space. Apply:
- **PCA** (linear): shows global variance structure, explained variance per component
- **t-SNE** (non-linear, sklearn): preserves local neighborhoods, reveals trajectory-like structure

### Analysis steps
1. Load `all_model_predictions.csv` (603 × ~100), fill NaN with 12 (noon)
2. Standardize features (StandardScaler)
3. PCA → 2 components, record explained variance
4. t-SNE → 2 components (perplexity=30, n_iter=1000, random_state=42)
5. Merge with metadata (condition, osd_id, tissue, light_regime, ecotype)
6. Visualize:
   - Panel A: t-SNE colored by condition (flight/ground)
   - Panel B: t-SNE colored by study
   - Panel C: PCA colored by condition with explained variance
   - Panel D: t-SNE colored by predicted CT (continuous color scale)
7. Quantitative test: compute per-study centroid distance between flight and ground clusters in t-SNE space; correlate with phase shift magnitude

### Output
- `figS7_circadian_trajectory.png/.svg` — 4-panel figure
- `data/trajectory_analysis.csv` — per-study centroid distances and correlation with phase shift

## Part 3: Gene Set Enrichment on Phase-Shifted Studies

### Studies
The 4 studies with significant phase shifts (p < 0.05) that have RNA-seq normalized counts:
- OSD-321 (22F/22G, whole_seedling, dark, phase shift −0.43h, p=1.8e-7)
- OSD-193 (16F/16G, root, unknown, phase shift −0.24h, p=0.002)
- OSD-38 (3F/9G, whole_seedling, dark, phase shift −0.52h, p=0.026)
- OSD-281 (16F/16G, root, unknown, phase shift −0.17h, p=0.038)

### DEG method: limma-voom + fgsea
1. For each study: load normalized counts, map GSM sample IDs to flight/ground using harmonized metadata
2. Filter low-expression genes (keep genes with CPM > 1 in at least min(group sizes) samples)
3. Apply TMM normalization (edgeR::calcNormFactors) on the count-like matrix
4. voom transform (limma::voom) with design matrix ~condition
5. Fit linear model (limma::eBayes), extract t-statistics per gene
6. Rank genes by t-statistic → ranked list for fgsea
7. Build Arabidopsis GO BP and KEGG pathway gene sets from org.At.tair.db
8. Run fgsea (fgsea::fgsea) with 10000 permutations, minSize=15, maxSize=500
9. Focus on circadian-related pathways: GO:0007623 (circadian rhythm), GO:0042752 (circadian regulation of gene expression), GO:0042754 (negative regulation of circadian rhythm), KEGG ath04710 (circadian rhythm)

### Output
- `figS8_enrichment.png/.svg` — enrichment results: Panel A = circadian pathway GSEA running enrichment scores for OSD-321; Panel B = top 10 enriched GO BP terms across studies (dot plot); Panel C = circadian gene set NES across 4 studies
- `tables/tableS3_enrichment_results.csv` — full fgsea results for all studies
- `data/deg_results/` — per-study DEG tables (gene, logFC, t, p)

### R implementation
Write an R script (`/workspace/spaceflight-circadian-decoder/src/enrichment_analysis.R`) that:
- Loads each study's normalized counts
- Runs limma-voom DEG
- Builds gene sets from org.At.tair.db
- Runs fgsea
- Saves results as CSV
- Generates the enrichment figure using ggplot2/enrichplot

## Part 4: Interactive Dashboard for GitHub Pages

### Approach
GitHub Pages serves static files only — no server-side processing. Create a self-contained HTML file with:
- Embedded Plotly.js (from CDN)
- All prediction/metadata data embedded as JSON (compressed if large)
- Multiple linked views with client-side filtering via Plotly's `updatemenus` and JavaScript callbacks

### Dashboard views
1. **Circadian phase overview**: Scatter plot of predicted CT vs circular variance, colored by condition. Dropdown filters: study, tissue, light regime.
2. **Phase shift forest plot**: Per-study phase shifts with CIs and pooled estimate. Hover shows full details.
3. **Circadian fingerprint**: Box plot of circular variance by condition, with study-level breakdown.
4. **Clock gene expression**: Heatmap of 11 core clock genes across samples (subset for performance).
5. **t-SNE trajectory**: The new trajectory analysis embedded as an interactive scatter.

### Deployment structure
```
docs/                          # GitHub Pages root (docs/ folder convention)
├── index.html                 # Main dashboard
├── data.js                    # Embedded JSON data
├── plotly.min.js              # Plotly.js (local copy, no CDN dependency)
└── README.md                  # Deployment instructions
```

### Output
- `docs/index.html` — self-contained interactive dashboard
- `docs/data.js` — embedded data as JavaScript variable
- `docs/plotly.min.js` — local Plotly.js copy
- `docs/README.md` — GitHub Pages deployment instructions
- Copy to `/mnt/results/` for user access

## Part 5: Manuscript Integration

### New content
1. **Supplementary Methods section**: "Circadian trajectory analysis" and "Gene set enrichment analysis"
2. **Supplementary Results section**: New findings from trajectory and enrichment analyses
3. **New figures**: figS7 (trajectory), figS8 (enrichment)
4. **New table**: tableS3 (enrichment results)
5. **Updated supplementary figure/table references** in main text

### LaTeX update
- Add new figures to `manuscript.tex` with proper captions
- Add new table
- Update supplementary methods and results sections
- Regenerate PDF (compile in /workspace/, copy to /mnt/results/)

### Repo update
- Add `src/enrichment_analysis.R` to repo
- Add `src/trajectory_analysis.py` to repo
- Add `docs/` folder for GitHub Pages dashboard
- Update README with new analysis descriptions
- Regenerate MANIFEST.txt
- Rebuild tar.gz archive

## Compute/Resource Estimate

| Task | Runtime | Memory | Disk |
|------|---------|--------|------|
| Figure fixes (re-run from cache) | ~2 min | <1 GB | negligible |
| Trajectory analysis (t-SNE on 603×100) | <1 min | <1 GB | negligible |
| DEG + enrichment (4 studies, limma-voom + fgsea) | ~10 min | ~2 GB | ~50 MB |
| Dashboard HTML generation | ~1 min | <1 GB | ~5 MB |
| Manuscript LaTeX compilation | ~2 min | <1 GB | ~5 MB |
| **Total** | **~16 min** | **<2 GB peak** | **~60 MB** |

All work runs on the default worker-0 machine. No HPC or additional machines needed.

## Execution Order

1. Fix all 12 figures (rewrite generate_figures.py, re-run from cached data)
2. Run trajectory analysis (new Python script)
3. Run enrichment analysis (new R script)
4. Generate interactive dashboard (new Python script generating HTML)
5. Update manuscript text + LaTeX, compile PDF
6. Update GitHub repo (add new files, regenerate archive)
7. Final verification of all deliverables

## Assumptions

- Cached prediction data (`/mnt/results/data/`) is valid and complete from prior analysis
- GeneLab normalized count files are available for the 4 significant RNA-seq studies
- org.At.tair.db GO/KEGG annotations are sufficient for Arabidopsis enrichment (installed and verified)
- GitHub Pages static HTML approach is acceptable (no server-side Dash)
- New supplementary figures are numbered S7 (trajectory) and S8 (enrichment)
- The existing manuscript structure can accommodate new supplementary sections without major restructuring
