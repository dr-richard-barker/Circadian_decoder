"""
Generate validation figure (fig2) and table (tableS2) using ChronoGauge's own test data.
Also generate figS1 (PCA/UMAP) from the actual predictions.
"""
import os
import sys
import json
import numpy as np
import pandas as pd

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

src_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, src_dir)
sys.path.insert(0, os.path.join(src_dir, 'figures'))
sys.path.insert(0, os.path.join(src_dir, 'tables'))

from chronogauge_apply import predict_ct, circular_error
from generate_figures import fig2_chronogauge_validation, figS1_pca_umap, save_figure
from generate_tables import tableS2_chronogauge_validation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.decomposition import PCA

rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
rcParams['svg.fonttype'] = 'none'

if os.path.exists("/workspace/ChronoGauge/data"):
    CG_DATA = "/workspace/ChronoGauge/data"
else:
    CG_DATA = os.path.abspath(os.path.join(src_dir, "..", "ChronoGauge", "data"))

if os.path.exists("/mnt/results"):
    RESULTS_DIR = "/mnt/results"
elif os.path.exists("/results"):
    RESULTS_DIR = "/results"
else:
    RESULTS_DIR = os.path.abspath(os.path.join(src_dir, ".."))

FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
TABLES_DIR = os.path.join(RESULTS_DIR, "tables")
DATA_OUT_DIR = os.path.join(RESULTS_DIR, "data")

if not os.path.exists(os.path.join(CG_DATA, 'expression_matrices/x_test_rna.csv')):
    print(f"ChronoGauge validation data not found at {CG_DATA}.")
    print("Skipping validation generation. Using pre-generated validation figures/tables.")
    sys.exit(0)

# ============================================================
# Fig 2: ChronoGauge validation on held-out test data
# ============================================================
print("Generating Fig 2: ChronoGauge validation...")

# Load test RNA-seq data
X_test = pd.read_csv(os.path.join(CG_DATA, 'expression_matrices/x_test_rna.csv'), index_col=0)
Y_test = pd.read_csv(os.path.join(CG_DATA, 'targets/target_test_rna.csv'), index_col=0)

# Run ChronoGauge with all 100 models
results, all_preds = predict_ct(X_test, platform='rnaseq', n_models=100)

# Compute errors
true_times = Y_test.iloc[:, 0].values % 24
pred_times = results['predicted_CT'].values
errors = circular_error(pred_times, true_times)

validation_df = pd.DataFrame({
    'sample_name': Y_test.index,
    'true_CT': true_times,
    'predicted_CT': pred_times,
    'error_minutes': errors,
    'circular_variance': results['circular_variance'].values,
    'osd_id': 'ChronoGauge_test'
})

print(f"  MAE: {np.mean(np.abs(errors)):.1f} min")
print(f"  Median error: {np.median(np.abs(errors)):.1f} min")
print(f"  RMSE: {np.sqrt(np.mean(errors**2)):.1f} min")

fig2_chronogauge_validation(validation_df, FIGURES_DIR)

# Table S2
tableS2_chronogauge_validation(validation_df, TABLES_DIR)

# ============================================================
# Fig S1: PCA on predicted CT values from GeneLab data
# ============================================================
print("\nGenerating Fig S1: PCA/UMAP...")

# Load predictions and metadata
predictions_df = pd.read_csv(os.path.join(DATA_OUT_DIR, 'all_predictions.csv'))
metadata_path = os.path.join(DATA_DIR, 'harmonized_metadata.csv')
if not os.path.exists(metadata_path):
    metadata_path = os.path.join(RESULTS_DIR, 'tables', 'tableS1_sample_metadata.csv')
metadata = pd.read_csv(metadata_path)

# Merge
merged = predictions_df.merge(metadata[['sample_name', 'osd_id', 'condition', 'tissue', 'ecotype']],
                               on=['sample_name', 'osd_id'], how='inner')

# Also load model-level predictions for PCA
model_preds = pd.read_csv(os.path.join(DATA_OUT_DIR, 'all_model_predictions.csv'))
# Select only numeric model prediction columns
model_cols = [c for c in model_preds.columns if c.startswith('0_') or (c.isdigit() if c else False)]
# Actually the columns are model IDs (integers as strings)
model_cols = [c for c in model_preds.columns if c not in ['sample_name', 'osd_id'] and c.replace('.0','').isdigit()]

if len(model_cols) > 5:
    # Use model predictions as features for PCA
    pca_features = model_preds[model_cols].fillna(12)  # Fill missing with noon
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(pca_features)

    # Merge with metadata
    pca_df = pd.DataFrame({
        'sample_name': model_preds['sample_name'],
        'osd_id': model_preds['osd_id'],
        'PC1': pcs[:, 0],
        'PC2': pcs[:, 1]
    })
    pca_df = pca_df.merge(metadata[['sample_name', 'osd_id', 'condition', 'tissue']],
                          on=['sample_name', 'osd_id'], how='inner')

    # Generate figure
    PHYLO_COLORS = {
        'blue': '#0279EE', 'orange': '#FF9400', 'green': '#75A025',
        'pink': '#FD9BED', 'yellow': '#E9ED4C', 'black': '#000000'
    }
    CB_PALETTE = ['#0279EE', '#E66726', '#2CA02C', '#D62728', '#9467BD',
                  '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: PCA colored by condition
    ax = axes[0]
    for cond, color in [('ground_control', PHYLO_COLORS['blue']), ('flight', PHYLO_COLORS['orange'])]:
        mask = pca_df['condition'] == cond
        ax.scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC2'],
                   c=color, label=cond, alpha=0.5, s=20, edgecolors='white', linewidth=0.3)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax.legend()
    ax.set_title('A', loc='left', fontweight='bold', fontsize=14)

    # Panel B: PCA colored by study
    ax = axes[1]
    studies = pca_df['osd_id'].unique()
    for i, study in enumerate(studies):
        mask = pca_df['osd_id'] == study
        ax.scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC2'],
                   c=CB_PALETTE[i % len(CB_PALETTE)], label=study, alpha=0.5, s=20)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax.legend(fontsize=6, ncol=2)
    ax.set_title('B', loc='left', fontweight='bold', fontsize=14)

    plt.tight_layout()
    save_figure(fig, FIGURES_DIR, 'figS1_pca')
else:
    print(f"  Not enough model prediction columns for PCA ({len(model_cols)} found)")

print("\nValidation figures and tables generated successfully.")
print(f"  MAE on held-out test: {np.mean(np.abs(errors)):.1f} min")
print(f"  Figures in: {FIGURES_DIR}")
print(f"  Tables in: {TABLES_DIR}")
