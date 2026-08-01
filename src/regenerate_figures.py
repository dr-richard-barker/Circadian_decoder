"""
Regenerate all figures from cached data.
Does NOT re-run ChronoGauge or statistical analysis — uses saved CSVs.
"""
import os
import sys
import json
import numpy as np
import pandas as pd

src_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, src_dir)
sys.path.insert(0, os.path.join(src_dir, 'figures'))
sys.path.insert(0, os.path.join(src_dir, 'tables'))

if os.path.exists("/mnt/results"):
    RESULTS_DIR = "/mnt/results"
elif os.path.exists("/results"):
    RESULTS_DIR = "/results"
else:
    RESULTS_DIR = os.path.abspath(os.path.join(src_dir, ".."))

FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
DATA_OUT_DIR = os.path.join(RESULTS_DIR, "data")

if os.path.exists("/workspace/genelab_data"):
    DATA_DIR = "/workspace/genelab_data"
else:
    DATA_DIR = os.path.abspath(os.path.join(src_dir, "..", "genelab_data"))

from generate_figures import (
    fig1_study_overview, fig3_phase_shift_polar, fig4_forest_plot,
    fig5_stratified_analysis, fig6_circadian_fingerprint,
    figS1_pca_umap, figS2_core_clock_genes, figS3_phase_by_hardware,
    figS4_phase_by_light_regime, figS5_fork_comparison, figS6_model_uncertainty
)

# Fork definitions
FORK_A_IDS = [7, 17, 37, 38, 44, 46, 120, 121, 136, 147, 205, 208, 218, 251]

def main():
    print("=" * 60)
    print("REGENERATING ALL FIGURES FROM CACHED DATA")
    print("=" * 60)

    # Load cached data
    metadata_path = os.path.join(DATA_DIR, 'harmonized_metadata.csv')
    if not os.path.exists(metadata_path):
        metadata_path = os.path.join(RESULTS_DIR, 'tables', 'tableS1_sample_metadata.csv')
    metadata = pd.read_csv(metadata_path)
    predictions_df = pd.read_csv(os.path.join(DATA_OUT_DIR, 'all_predictions.csv'))
    model_preds_df = pd.read_csv(os.path.join(DATA_OUT_DIR, 'all_model_predictions.csv'))
    per_study = pd.read_csv(os.path.join(DATA_OUT_DIR, 'per_study_results.csv'))
    clock_expr_df = pd.read_csv(os.path.join(DATA_OUT_DIR, 'clock_gene_expression.csv'))

    with open(os.path.join(DATA_OUT_DIR, 'meta_analysis_results.json')) as f:
        meta_results = json.load(f)

    print(f"Loaded: {len(metadata)} metadata, {len(predictions_df)} predictions, {len(per_study)} per-study")

    # Main figures
    print("\n--- Main figures ---")
    fig1_study_overview(metadata, per_study, FIGURES_DIR)
    fig3_phase_shift_polar(per_study, FIGURES_DIR)
    fig4_forest_plot(per_study, meta_results, FIGURES_DIR)
    fig5_stratified_analysis(per_study, meta_results['tissue'], FIGURES_DIR, stratify_by='tissue')
    fig6_circadian_fingerprint(predictions_df, metadata, FIGURES_DIR)

    # Supplementary figures
    print("\n--- Supplementary figures ---")
    # figS1: PCA on model predictions
    figS1_pca_umap(model_preds_df, metadata, FIGURES_DIR)

    # figS2: clock gene heatmap
    if len(clock_expr_df) > 0:
        figS2_core_clock_genes(clock_expr_df, metadata, FIGURES_DIR)

    figS3_phase_by_hardware(per_study, FIGURES_DIR)
    figS4_phase_by_light_regime(per_study, FIGURES_DIR)

    # figS5: fork comparison
    fork_a_set = {f'OSD-{x}' for x in FORK_A_IDS}
    fork_a_results = per_study[per_study['osd_id'].isin(fork_a_set)]
    fork_b_results = per_study.copy()  # Fork B = all
    figS5_fork_comparison(fork_a_results, fork_b_results, FIGURES_DIR)

    # figS6: model uncertainty
    figS6_model_uncertainty(per_study, predictions_df, FIGURES_DIR)

    print("\n" + "=" * 60)
    print("ALL FIGURES REGENERATED")
    print("=" * 60)
    print(f"Output: {FIGURES_DIR}")

    # Verify all files exist
    expected = [
        'fig1_study_overview', 'fig3_phase_shift_polar', 'fig4_forest_plot',
        'fig5_stratified_tissue', 'fig6_circadian_fingerprint',
        'figS1_pca', 'figS2_clock_genes', 'figS3_phase_by_hardware',
        'figS4_phase_by_light_regime', 'figS5_fork_comparison', 'figS6_model_uncertainty'
    ]
    for name in expected:
        for ext in ['.png', '.svg']:
            path = os.path.join(FIGURES_DIR, name + ext)
            if os.path.exists(path):
                size = os.path.getsize(path)
                print(f"  {name}{ext}: {size:,} bytes")
            else:
                print(f"  MISSING: {name}{ext}")


if __name__ == "__main__":
    main()
