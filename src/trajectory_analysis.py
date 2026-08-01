"""
Circadian trajectory analysis.
Uses the 100 ChronoGauge sub-model predictions per sample as a circadian
fingerprint feature space. Applies PCA and t-SNE to reveal trajectory-like
structure, then quantifies per-study flight-ground centroid separation.

Generates:
  - figS7_circadian_trajectory.png/.svg (4-panel figure)
  - trajectory_analysis.csv (per-study centroid distances + correlation)
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from scipy import stats

rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
rcParams['svg.fonttype'] = 'none'
rcParams['pdf.fonttype'] = 42
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False

PHYLO_COLORS = {
    'blue': '#0279EE', 'orange': '#FF9400', 'green': '#75A025',
    'pink': '#FD9BED', 'yellow': '#E9ED4C', 'black': '#000000',
}
CB_PALETTE = ['#0279EE', '#E66726', '#2CA02C', '#D62728', '#9467BD',
              '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']

src_dir = os.path.dirname(os.path.abspath(__file__))

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


def save_figure(fig, output_dir, name, dpi=300):
    os.makedirs(output_dir, exist_ok=True)
    svg_path = os.path.join(output_dir, f'{name}.svg')
    png_path = os.path.join(output_dir, f'{name}.png')
    fig.savefig(svg_path, format='svg', bbox_inches='tight', dpi=dpi)
    fig.savefig(png_path, format='png', bbox_inches='tight', dpi=dpi)
    plt.close(fig)
    print(f"  Saved {name}.svg and {name}.png")
    return svg_path, png_path


def main():
    print("=" * 60)
    print("CIRCADIAN TRAJECTORY ANALYSIS")
    print("=" * 60)

    # Load data
    model_preds = pd.read_csv(os.path.join(DATA_OUT_DIR, 'all_model_predictions.csv'))
    predictions = pd.read_csv(os.path.join(DATA_OUT_DIR, 'all_predictions.csv'))
    metadata_path = os.path.join(DATA_DIR, 'harmonized_metadata.csv')
    if not os.path.exists(metadata_path):
        metadata_path = os.path.join(RESULTS_DIR, 'tables', 'tableS1_sample_metadata.csv')
    metadata = pd.read_csv(metadata_path)
    per_study = pd.read_csv(os.path.join(DATA_OUT_DIR, 'per_study_results.csv'))

    print(f"Model predictions: {model_preds.shape}")

    # Extract feature columns (model IDs — integers as strings)
    meta_cols = ['sample_name', 'osd_id']
    feature_cols = [c for c in model_preds.columns if c not in meta_cols]
    print(f"Feature columns (sub-models): {len(feature_cols)}")

    # Build feature matrix — fill NaN with 12 (noon)
    X = model_preds[feature_cols].fillna(12).values
    print(f"Feature matrix: {X.shape}")

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # PCA
    print("\nRunning PCA...")
    pca = PCA(n_components=10)
    pcs = pca.fit_transform(X_scaled)
    print(f"  PC1: {pca.explained_variance_ratio_[0]*100:.1f}%")
    print(f"  PC2: {pca.explained_variance_ratio_[1]*100:.1f}%")
    print(f"  PC1-10 cumulative: {sum(pca.explained_variance_ratio_)*100:.1f}%")

    # t-SNE
    print("\nRunning t-SNE...")
    perplexity = min(30, len(X) - 1)
    tsne = TSNE(n_components=2, perplexity=perplexity, max_iter=1000,
                random_state=42, metric='euclidean')
    tsne_results = tsne.fit_transform(X_scaled)
    print(f"  t-SNE shape: {tsne_results.shape}")

    # Merge with metadata
    merged = model_preds[['sample_name', 'osd_id']].copy()
    merged['PC1'] = pcs[:, 0]
    merged['PC2'] = pcs[:, 1]
    merged['tSNE1'] = tsne_results[:, 0]
    merged['tSNE2'] = tsne_results[:, 1]

    # Merge with condition and other metadata
    meta_subset = metadata[['sample_name', 'osd_id', 'condition', 'tissue',
                            'light_regime', 'ecotype']].drop_duplicates(
                            subset=['sample_name', 'osd_id'])
    merged = merged.merge(meta_subset, on=['sample_name', 'osd_id'], how='left')

    # Merge with predicted CT for continuous coloring
    pred_subset = predictions[['sample_name', 'osd_id', 'predicted_CT', 'circular_variance']]
    merged = merged.merge(pred_subset, on=['sample_name', 'osd_id'], how='left')

    print(f"Merged data: {len(merged)} samples")
    print(f"  Flight: {(merged['condition']=='flight').sum()}")
    print(f"  Ground: {(merged['condition']=='ground_control').sum()}")

    # ============================================================
    # Quantitative analysis: per-study centroid distances
    # ============================================================
    print("\nComputing per-study centroid distances...")
    trajectory_results = []

    for osd_id, study_data in merged.groupby('osd_id'):
        flight = study_data[study_data['condition'] == 'flight']
        ground = study_data[study_data['condition'] == 'ground_control']

        if len(flight) < 2 or len(ground) < 2:
            continue

        # Centroids in t-SNE space
        f_centroid = flight[['tSNE1', 'tSNE2']].mean()
        g_centroid = ground[['tSNE1', 'tSNE2']].mean()
        tsne_dist = np.sqrt((f_centroid['tSNE1'] - g_centroid['tSNE1'])**2 +
                            (f_centroid['tSNE2'] - g_centroid['tSNE2'])**2)

        # Centroids in PCA space
        f_centroid_pca = flight[['PC1', 'PC2']].mean()
        g_centroid_pca = ground[['PC1', 'PC2']].mean()
        pca_dist = np.sqrt((f_centroid_pca['PC1'] - g_centroid_pca['PC1'])**2 +
                           (f_centroid_pca['PC2'] - g_centroid_pca['PC2'])**2)

        # Within-group dispersion (mean pairwise distance)
        def mean_pairwise_dist(df, cols):
            from scipy.spatial.distance import pdist
            return np.mean(pdist(df[cols].values)) if len(df) > 1 else 0

        f_disp = mean_pairwise_dist(flight, ['tSNE1', 'tSNE2'])
        g_disp = mean_pairwise_dist(ground, ['tSNE1', 'tSNE2'])
        avg_disp = (f_disp + g_disp) / 2

        # Separation ratio (between-group / within-group)
        sep_ratio = tsne_dist / avg_disp if avg_disp > 0 else np.nan

        # Get phase shift from per_study
        ps_row = per_study[per_study['osd_id'] == osd_id]
        phase_shift = ps_row['phase_shift_hours'].values[0] if len(ps_row) > 0 else np.nan
        p_value = ps_row['p_value'].values[0] if len(ps_row) > 0 else np.nan

        trajectory_results.append({
            'osd_id': osd_id,
            'n_flight': len(flight),
            'n_ground': len(ground),
            'tsne_centroid_distance': tsne_dist,
            'pca_centroid_distance': pca_dist,
            'flight_dispersion': f_disp,
            'ground_dispersion': g_disp,
            'separation_ratio': sep_ratio,
            'phase_shift_hours': phase_shift,
            'p_value': p_value,
        })

    traj_df = pd.DataFrame(trajectory_results)
    traj_df.to_csv(os.path.join(DATA_OUT_DIR, 'trajectory_analysis.csv'), index=False)
    print(f"Saved trajectory_analysis.csv ({len(traj_df)} studies)")

    # Correlation between centroid distance and phase shift
    valid = traj_df.dropna(subset=['tsne_centroid_distance', 'phase_shift_hours'])
    if len(valid) >= 4:
        r_tsne, p_tsne = stats.spearmanr(valid['tsne_centroid_distance'],
                                          valid['phase_shift_hours'].abs())
        r_pca, p_pca = stats.spearmanr(valid['pca_centroid_distance'],
                                        valid['phase_shift_hours'].abs())
        print(f"\nCorrelation (centroid distance vs |phase shift|):")
        print(f"  t-SNE: rho={r_tsne:.3f}, p={p_tsne:.3f}")
        print(f"  PCA:   rho={r_pca:.3f}, p={p_pca:.3f}")

    # ============================================================
    # Figure S7: 4-panel trajectory figure
    # ============================================================
    print("\nGenerating figS7...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Panel A: t-SNE colored by condition
    ax = axes[0, 0]
    for cond, color, label in [('ground_control', PHYLO_COLORS['blue'], 'Ground control'),
                                ('flight', PHYLO_COLORS['orange'], 'Spaceflight')]:
        mask = merged['condition'] == cond
        ax.scatter(merged.loc[mask, 'tSNE1'], merged.loc[mask, 'tSNE2'],
                   c=color, label=label, alpha=0.5, s=25, edgecolors='white', linewidth=0.3)
    ax.set_xlabel('t-SNE 1')
    ax.set_ylabel('t-SNE 2')
    ax.legend(framealpha=0.9)
    ax.set_title('A  t-SNE by condition', loc='left', fontweight='bold', fontsize=13)

    # Panel B: t-SNE colored by study (colorbar)
    ax = axes[0, 1]
    studies = sorted(merged['osd_id'].unique())
    study_to_idx = {s: i for i, s in enumerate(studies)}
    study_indices = merged['osd_id'].map(study_to_idx).values
    scatter = ax.scatter(merged['tSNE1'], merged['tSNE2'], c=study_indices,
                         cmap='tab20', alpha=0.5, s=25, edgecolors='none')
    ax.set_xlabel('t-SNE 1')
    ax.set_ylabel('t-SNE 2')
    ax.set_title('B  t-SNE by study', loc='left', fontweight='bold', fontsize=13)
    cbar = plt.colorbar(scatter, ax=ax, ticks=range(len(studies)), pad=0.02)
    cbar.ax.set_yticklabels(studies, fontsize=5)
    cbar.ax.set_title('Study', fontsize=8)

    # Panel C: PCA colored by condition
    ax = axes[1, 0]
    for cond, color, label in [('ground_control', PHYLO_COLORS['blue'], 'Ground control'),
                                ('flight', PHYLO_COLORS['orange'], 'Spaceflight')]:
        mask = merged['condition'] == cond
        ax.scatter(merged.loc[mask, 'PC1'], merged.loc[mask, 'PC2'],
                   c=color, label=label, alpha=0.5, s=25, edgecolors='white', linewidth=0.3)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax.legend(framealpha=0.9)
    ax.set_title('C  PCA by condition', loc='left', fontweight='bold', fontsize=13)

    # Panel D: t-SNE colored by predicted CT (continuous)
    ax = axes[1, 1]
    scatter = ax.scatter(merged['tSNE1'], merged['tSNE2'],
                         c=merged['predicted_CT'], cmap='twilight_shifted',
                         alpha=0.6, s=25, edgecolors='white', linewidth=0.3,
                         vmin=0, vmax=24)
    ax.set_xlabel('t-SNE 1')
    ax.set_ylabel('t-SNE 2')
    ax.set_title('D  t-SNE by predicted CT', loc='left', fontweight='bold', fontsize=13)
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label('Predicted CT (h)')

    plt.tight_layout()
    save_figure(fig, FIGURES_DIR, 'figS7_circadian_trajectory')

    # ============================================================
    # Supplementary: centroid distance vs phase shift scatter
    # ============================================================
    if len(valid) >= 4:
        fig2, ax = plt.subplots(figsize=(7, 5))
        ax.scatter(valid['tsne_centroid_distance'], valid['phase_shift_hours'].abs(),
                   c=PHYLO_COLORS['blue'], alpha=0.7, s=60, edgecolors='white', linewidth=0.3)
        for _, row in valid.iterrows():
            ax.annotate(row['osd_id'], (row['tsne_centroid_distance'], abs(row['phase_shift_hours'])),
                        fontsize=6, ha='left', va='bottom',
                        xytext=(3, 3), textcoords='offset points')
        ax.set_xlabel('t-SNE centroid distance (flight vs ground)')
        ax.set_ylabel('|Phase shift| (h)')
        ax.set_title(f'Centroid separation vs phase shift (Spearman rho={r_tsne:.2f}, p={p_tsne:.3f})')
        plt.tight_layout()
        save_figure(fig2, FIGURES_DIR, 'figS7b_centroid_vs_shift')

    print("\n" + "=" * 60)
    print("TRAJECTORY ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
