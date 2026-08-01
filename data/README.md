# Data Dictionary — Spaceflight Circadian Decoder

This directory contains the curated metadata, predictions, and results of the circadian clock phase analysis for spaceflight-grown *Arabidopsis thaliana* samples.

## Files

### `all_predictions.csv`
Per-sample aggregated circadian time predictions.
- **`sample_name`**: NASA GeneLab sample identifier (matching ISA-Tab)
- **`predicted_CT`**: Circular mean of the 100 sub-model predictions (hours, range [0, 24])
- **`circular_variance`**: Circular variance of predictions ($1-R$, range [0, 1]). Lower values indicate higher ensemble consensus (clock robustness).
- **`n_models_used`**: Number of ChronoGauge sub-models utilized for prediction (sub-models are skipped if <80% of their feature genes are present in the expression matrix).
- **`osd_id`**: NASA Open Science Data Repository study accession (e.g., `OSD-321`)

### `all_model_predictions.csv`
Raw prediction outputs from each of the 100 individual neural network sub-models in the ChronoGauge ensemble across all 603 samples. Columns `0` to `99` represent the sub-model IDs.

### `clock_gene_expression.csv`
Z-scored expression values of 11 core circadian clock genes used for visualization and validation.
- **`gene`**: AGI identifier (e.g., `AT2G46830`)
- **`sample_name`**: Sample identifier
- **`expression`**: Z-scored expression value (normalized per gene across samples within each study)
- **`gene_name`**: Common gene symbol (e.g., `CCA1`, `LHY`, `TOC1`)
- **`osd_id`**: Study accession

### `per_study_results.csv`
Circadian phase shift statistics and Watson-Williams test results for the 18 studies containing both flight and ground controls.
- **`osd_id`**: Study accession
- **`n_flight`**: Number of spaceflight replicates
- **`n_ground`**: Number of ground control replicates
- **`phase_shift_hours`**: Phase shift between conditions ($Flight - Ground$, circular difference in range $[-12, 12]$ hours). Negative values indicate a phase advance.
- **`ci_lower`**: Lower bound of the 95% bootstrap confidence interval
- **`ci_upper`**: Upper bound of the 95% bootstrap confidence interval
- **`p_value`**: Watson-Williams test p-value (test for equality of circular means)
- **`tissue`**: Curated tissue type
- **`light_regime`**: Curated growth light regime (dark, continuous light, photoperiod, or unknown)
- **`hardware`**: Curated hardware configuration
- **`included_in_meta`**: Boolean indicating if the study was included in the meta-analysis (requires at least 3 replicates in both conditions)

### `trajectory_analysis.csv`
Centroid distances in t-SNE/PCA space and Spearman correlation results.
- **`osd_id`**: Study accession
- **`tsne_centroid_distance`**: Euclidean distance between flight and ground control centroids in the 2D t-SNE embedding
- **`pca_centroid_distance`**: Euclidean distance in the 2D PCA embedding
- **`absolute_phase_shift`**: Absolute value of the circular phase shift (hours)

### `meta_analysis_results.json`
Structured JSON output of the random-effects meta-analysis performed using R `metafor`. Contains effect sizes (pooled phase shifts), 95% confidence intervals, p-values, Cochran's Q heterogeneity test statistics, and $I^2$ percentages, overall and stratified by tissue, light regime, and hardware.

### `deg_results/`
Directory containing differential expression analysis results (flight vs ground) generated using `limma-voom` for the 4 studies showing statistically significant phase shifts:
- `OSD-38_deg.csv`
- `OSD-193_deg.csv`
- `OSD-281_deg.csv`
- `OSD-321_deg.csv`
Columns include gene identifier, log2 fold change (`logFC`), average expression (`AveExpr`), t-statistic (`t`), raw p-value (`P.Value`), and false discovery rate (`adj.P.Val`).
