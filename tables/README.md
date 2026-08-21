# Tables Index — Spaceflight Circadian Decoder

This directory contains the main and supplementary tables corresponding to the manuscript and dashboard.

## Table List

| Table File | Type | Description |
|------------|------|-------------|
| `table1_dataset_characteristics.csv` | Main | Summary characteristics of the 24 Arabidopsis spaceflight studies included in ChronoGauge analysis (and 26 total catalogued in OSDR). Details assay platform, tissue type, ecotype, hardware configuration, light regime, and replicate counts for flight and ground. |
| `table2_per_study_results.csv` | Main | Per-study circadian phase shifts (hours), bootstrap 95% confidence intervals, and Watson-Williams test p-values comparing flight and ground control sample distributions. |
| `table3_meta_analysis_overall.csv` | Main | Pooled random-effects meta-analysis (REML) estimates overall and for Fork A/B comparison. |
| `table3_meta_analysis_tissue.csv` | Main | Pooled estimates stratified by tissue type (root, whole seedling, leaf, shoot, hypocotyl). |
| `table3_meta_analysis_light_regime.csv` | Main | Pooled estimates stratified by growth light regime (dark, continuous light, photoperiod). |
| `table3_meta_analysis_hardware.csv` | Main | Pooled estimates stratified by hardware configuration. |
| `tableS1_sample_metadata.csv` | Supp | Complete harmonized sample-level metadata for all 1,169 samples parsed from ISA-Tab. Includes ground-only studies (OSD-46, OSD-136, OSD-208) and flight-only studies (OSD-251, OSD-346) which were excluded from phase-shift comparisons. |
| `tableS2_chronogauge_validation.csv` | Supp | ChronoGauge time-prediction validation statistics (Mean Absolute Error, Median Error, RMSE, circular variance) evaluated on held-out constant-light RNA-seq test data. |
| `tableS3_enrichment_results.csv` | Supp | Full fast Gene Set Enrichment Analysis (fgsea) output across 816 GO Biological Process terms and KEGG pathways for the 4 studies with significant phase shifts. |
