"""
Figure generation module for the spaceflight circadian analysis.
Generates all main and supplementary figures as SVG + PNG.

v2: Fixed readability issues — legends outside plots, labels not clipped,
    figS6 panel A uses actual circular variance data.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Wedge
from matplotlib.lines import Line2D
from pathlib import Path

# Phylo color palette
PHYLO_COLORS = {
    'black': '#000000',
    'cream': '#ECE9E2',
    'white': '#FAF9F3',
    'yellow': '#E9ED4C',
    'orange': '#FF9400',
    'green': '#75A025',
    'pink': '#FD9BED',
    'blue': '#0279EE',
}

# Colorblind-friendly palette
CB_PALETTE = ['#0279EE', '#E66726', '#2CA02C', '#D62728', '#9467BD',
              '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']

# Set global style
rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
rcParams['svg.fonttype'] = 'none'  # Keep SVG text editable
rcParams['pdf.fonttype'] = 42
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False


def save_figure(fig, output_dir, name, dpi=300):
    """Save figure as both SVG and PNG."""
    os.makedirs(output_dir, exist_ok=True)
    svg_path = os.path.join(output_dir, f'{name}.svg')
    png_path = os.path.join(output_dir, f'{name}.png')
    fig.savefig(svg_path, format='svg', bbox_inches='tight', dpi=dpi)
    fig.savefig(png_path, format='png', bbox_inches='tight', dpi=dpi)
    plt.close(fig)
    print(f"  Saved {name}.svg and {name}.png")
    return svg_path, png_path


# ============================================================
# MAIN FIGURES
# ============================================================

def fig1_study_overview(metadata_df, per_study_df, output_dir):
    """Figure 1: Study overview - dataset characteristics and design."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel A: Studies by platform and tissue
    ax = axes[0]
    platform_tissue = metadata_df.groupby(['assay_technology', 'tissue']).size().unstack(fill_value=0)
    platform_tissue.plot(kind='bar', stacked=True, ax=ax, color=CB_PALETTE[:len(platform_tissue.columns)])
    ax.set_title('A', loc='left', fontweight='bold', fontsize=14)
    ax.set_xlabel('Platform')
    ax.set_ylabel('Number of samples')
    ax.legend(title='Tissue', fontsize=7, bbox_to_anchor=(1.02, 1), loc='upper left', framealpha=0.9)

    # Panel B: Samples by condition and hardware
    ax = axes[1]
    cond_hw = metadata_df.groupby(['condition', 'hardware']).size().unstack(fill_value=0)
    cond_hw.plot(kind='bar', stacked=True, ax=ax, color=CB_PALETTE[:len(cond_hw.columns)])
    ax.set_title('B', loc='left', fontweight='bold', fontsize=14)
    ax.set_xlabel('Condition')
    ax.set_ylabel('Number of samples')
    ax.legend(title='Hardware', fontsize=6, bbox_to_anchor=(1.02, 1), loc='upper left', framealpha=0.9)

    # Panel C: Ecotype distribution
    ax = axes[2]
    ecotype_counts = metadata_df['ecotype'].value_counts().head(10)
    ecotype_counts.plot(kind='barh', ax=ax, color=PHYLO_COLORS['blue'])
    ax.set_title('C', loc='left', fontweight='bold', fontsize=14)
    ax.set_xlabel('Number of samples')
    ax.set_ylabel('Ecotype')

    plt.tight_layout()
    return save_figure(fig, output_dir, 'fig1_study_overview')


def fig2_chronogauge_validation(validation_df, output_dir):
    """Figure 2: ChronoGauge validation - predicted vs known ZT."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Predicted vs true CT scatter
    ax = axes[0]
    ax.scatter(validation_df['true_CT'], validation_df['predicted_CT'],
               alpha=0.5, s=30, c=PHYLO_COLORS['blue'], edgecolors='white', linewidth=0.3)
    # Diagonal
    ax.plot([0, 24], [0, 24], 'k--', alpha=0.3, linewidth=1)
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 24)
    ax.set_xticks(range(0, 25, 4))
    ax.set_yticks(range(0, 25, 4))
    ax.set_xlabel('Known circadian time (h)')
    ax.set_ylabel('Predicted circadian time (h)')
    ax.set_title('A', loc='left', fontweight='bold', fontsize=14)

    # Panel B: Error distribution
    ax = axes[1]
    errors = validation_df['error_minutes'].values
    ax.hist(errors, bins=30, color=PHYLO_COLORS['green'], alpha=0.7, edgecolor='white')
    ax.axvline(0, color='k', linestyle='--', linewidth=1)
    ax.axvline(np.median(errors), color=PHYLO_COLORS['orange'], linestyle='-', linewidth=2,
               label=f'Median = {np.median(errors):.0f} min')
    ax.set_xlabel('Circular error (min)')
    ax.set_ylabel('Count')
    ax.legend(framealpha=0.9)
    ax.set_title('B', loc='left', fontweight='bold', fontsize=14)

    plt.tight_layout()
    return save_figure(fig, output_dir, 'fig2_chronogauge_validation')


def fig3_phase_shift_polar(per_study_df, output_dir):
    """Figure 3: Per-study phase shifts on polar plot."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 8), subplot_kw={'projection': 'polar'})

    included = per_study_df[per_study_df['included_in_meta']].copy()

    # Convert phase shifts to radians
    angles = included['phase_shift_hours'].values * 2 * np.pi / 24
    pvals = included['p_value'].values

    # Color by significance
    colors = [PHYLO_COLORS['orange'] if p < 0.05 else PHYLO_COLORS['blue'] for p in pvals]

    # Plot as points on polar
    for i, (angle, color, p) in enumerate(zip(angles, colors, pvals)):
        r = 1 - included.iloc[i].get('p_value', 0.5)  # radius by significance
        ax.scatter(angle, r, c=color, s=80, zorder=5, edgecolors='white', linewidth=0.5)

    # Reference line at 0 (no shift)
    ax.plot([0, 0], [0, 1], 'k--', alpha=0.3, linewidth=1)

    ax.set_theta_zero_location('E')
    ax.set_theta_direction(1)
    ax.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
    ax.set_xticklabels(['0h', '3h', '6h', '9h', '12h', '15h', '18h', '21h'])
    ax.set_title('Per-study circadian phase shifts (flight vs ground)', pad=20)

    # Legend — positioned to avoid clipping
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=PHYLO_COLORS['orange'],
               markersize=8, label='p < 0.05'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=PHYLO_COLORS['blue'],
               markersize=8, label='p >= 0.05'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1.05),
              framealpha=0.9)

    plt.tight_layout()
    return save_figure(fig, output_dir, 'fig3_phase_shift_polar')


def fig4_forest_plot(per_study_df, meta_results, output_dir):
    """Figure 4: Forest plot of phase shifts with meta-analysis summary.
    Uses a two-column layout: text labels on left, plot on right.
    """
    included = per_study_df[per_study_df['included_in_meta']].copy().sort_values('phase_shift_hours')

    n_studies = len(included)

    # Use gridspec for text column + plot column
    fig = plt.figure(figsize=(12, max(6, n_studies * 0.4 + 2)))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 2], wspace=0.05)
    ax_text = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1])

    y_positions = range(n_studies)

    for i, (_, row) in enumerate(included.iterrows()):
        color = PHYLO_COLORS['orange'] if row['p_value'] < 0.05 else PHYLO_COLORS['blue']
        ax.errorbar(row['phase_shift_hours'], i,
                    xerr=[[row['phase_shift_hours'] - row['phase_shift_ci_lower']],
                          [row['phase_shift_ci_upper'] - row['phase_shift_hours']]],
                    fmt='o', color=color, capsize=3, markersize=6, linewidth=1.5)

    # Meta-analysis summary
    if 'overall' in meta_results and 'POOLED_EFFECT' in meta_results['overall']:
        m = meta_results['overall']
        ax.errorbar(m['POOLED_EFFECT'], n_studies + 0.5,
                    xerr=[[m['POOLED_EFFECT'] - m['CI_LOWER']], [m['CI_UPPER'] - m['POOLED_EFFECT']]],
                    fmt='D', color=PHYLO_COLORS['green'], capsize=5, markersize=10, linewidth=2)

    ax.axvline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Phase shift (hours, flight - ground)')
    ax.set_yticks([])
    ax.set_ylim(-0.5, n_studies + 1.5)
    ax.set_title('Circadian phase shift: spaceflight vs ground control', pad=10)

    # Text column — labels aligned with plot rows
    ax_text.set_xlim(0, 1)
    ax_text.set_ylim(-0.5, n_studies + 1.5)
    ax_text.axis('off')

    for i, (_, row) in enumerate(included.iterrows()):
        ax_text.text(0.98, i, f"{row['osd_id']}  (n={int(row['n_flight'])}F/{int(row['n_ground'])}G)",
                     ha='right', va='center', fontsize=8, family='monospace')

    if 'overall' in meta_results and 'POOLED_EFFECT' in meta_results['overall']:
        m = meta_results['overall']
        ax_text.text(0.98, n_studies + 0.5,
                     f"Pooled  (p={m['P_VALUE']:.3f}, I2={m['I2']:.0f}%)",
                     ha='right', va='center', fontsize=9, fontweight='bold', family='monospace')

    return save_figure(fig, output_dir, 'fig4_forest_plot')


def fig5_stratified_analysis(per_study_df, meta_results, output_dir, stratify_by='tissue'):
    """Figure 5: Stratified meta-analysis by tissue/genotype/hardware."""
    included = per_study_df[per_study_df['included_in_meta']].copy()

    strata = included[stratify_by].unique()
    strata = [s for s in strata if s and s != '']

    n_strata = len(strata)
    fig, axes = plt.subplots(1, n_strata, figsize=(5 * n_strata, 5), squeeze=False)

    for idx, stratum in enumerate(strata):
        ax = axes[0, idx]
        stratum_data = included[included[stratify_by] == stratum]

        for i, (_, row) in enumerate(stratum_data.iterrows()):
            color = PHYLO_COLORS['orange'] if row['p_value'] < 0.05 else PHYLO_COLORS['blue']
            ax.errorbar(row['phase_shift_hours'], i,
                        xerr=[[row['phase_shift_hours'] - row['phase_shift_ci_lower']],
                              [row['phase_shift_ci_upper'] - row['phase_shift_hours']]],
                        fmt='o', color=color, capsize=3, markersize=6)

        # Meta-analysis for this stratum
        key = str(stratum)
        if key in meta_results and 'POOLED_EFFECT' in meta_results[key]:
            m = meta_results[key]
            ax.errorbar(m['POOLED_EFFECT'], len(stratum_data) + 0.5,
                        xerr=[[m['POOLED_EFFECT'] - m['CI_LOWER']], [m['CI_UPPER'] - m['POOLED_EFFECT']]],
                        fmt='D', color=PHYLO_COLORS['green'], capsize=5, markersize=8)

        ax.axvline(0, color='k', linestyle='--', alpha=0.3)
        ax.set_xlabel('Phase shift (h)')
        ax.set_title(stratum, fontweight='bold')
        ax.set_yticks([])

    fig.suptitle(f'Phase shift stratified by {stratify_by}', fontsize=14, y=0.98)
    plt.subplots_adjust(top=0.88)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    return save_figure(fig, output_dir, f'fig5_stratified_{stratify_by}')


def fig6_circadian_fingerprint(predictions_df, metadata_df, output_dir):
    """Figure 6: Circadian fingerprint - circular variance comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    merged = predictions_df.merge(metadata_df[['sample_name', 'condition', 'osd_id']], on=['sample_name', 'osd_id'])

    # Panel A: Circular variance by condition
    ax = axes[0]
    flight = merged[merged['condition'] == 'flight']['circular_variance'].dropna()
    ground = merged[merged['condition'] == 'ground_control']['circular_variance'].dropna()

    bp = ax.boxplot([ground, flight], tick_labels=['Ground', 'Flight'],
                    patch_artist=True, widths=0.5)
    bp['boxes'][0].set_facecolor(PHYLO_COLORS['blue'])
    bp['boxes'][1].set_facecolor(PHYLO_COLORS['orange'])
    for box in bp['boxes']:
        box.set_alpha(0.7)
    ax.set_ylabel('Circular variance (1 - R)')
    ax.set_title('A', loc='left', fontweight='bold', fontsize=14)

    # Add individual points
    for i, data in enumerate([ground, flight]):
        x = np.random.normal(i + 1, 0.04, size=len(data))
        ax.scatter(x, data, alpha=0.3, s=10, c='black')

    # Panel B: Circular variance by study
    ax = axes[1]
    study_var = merged.groupby(['osd_id', 'condition'])['circular_variance'].mean().unstack()
    study_var = study_var.dropna()
    if len(study_var) > 0:
        x = np.arange(len(study_var))
        width = 0.35
        ax.bar(x - width/2, study_var['ground_control'], width, label='Ground',
               color=PHYLO_COLORS['blue'], alpha=0.7)
        ax.bar(x + width/2, study_var['flight'], width, label='Flight',
               color=PHYLO_COLORS['orange'], alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(study_var.index, rotation=45, ha='right', fontsize=7)
        ax.legend(framealpha=0.9)
    ax.set_ylabel('Mean circular variance')
    ax.set_title('B', loc='left', fontweight='bold', fontsize=14)

    plt.tight_layout()
    return save_figure(fig, output_dir, 'fig6_circadian_fingerprint')


# ============================================================
# SUPPLEMENTARY FIGURES
# ============================================================

def figS1_pca_umap(expression_df, metadata_df, output_dir):
    """Supplementary Figure S1: PCA of circadian predictions.
    Panel B uses a colorbar instead of a 23-entry legend to avoid overlap.
    """
    from sklearn.decomposition import PCA

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # PCA on predicted CT values
    merged = expression_df.merge(metadata_df[['sample_name', 'condition', 'osd_id']], on=['sample_name', 'osd_id'])

    # Panel A: PCA colored by condition
    ax = axes[0]
    features = merged.select_dtypes(include=[np.number]).drop(columns=['sample_name'], errors='ignore')
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(features.fillna(0))
    for cond, color in [('ground_control', PHYLO_COLORS['blue']), ('flight', PHYLO_COLORS['orange'])]:
        mask = merged['condition'] == cond
        ax.scatter(pcs[mask, 0], pcs[mask, 1], c=color, label=cond, alpha=0.5, s=20)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax.legend(framealpha=0.9)
    ax.set_title('A', loc='left', fontweight='bold', fontsize=14)

    # Panel B: PCA colored by study — use colorbar instead of legend
    ax = axes[1]
    studies = sorted(merged['osd_id'].unique())
    study_to_idx = {s: i for i, s in enumerate(studies)}
    study_indices = merged['osd_id'].map(study_to_idx).values
    scatter = ax.scatter(pcs[:, 0], pcs[:, 1], c=study_indices, cmap='tab20',
                         alpha=0.5, s=20, edgecolors='none')
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax.set_title('B', loc='left', fontweight='bold', fontsize=14)

    # Add colorbar with study labels
    cbar = plt.colorbar(scatter, ax=ax, ticks=range(len(studies)), pad=0.02)
    cbar.ax.set_yticklabels(studies, fontsize=6)
    cbar.ax.set_title('Study', fontsize=8)

    plt.tight_layout()
    return save_figure(fig, output_dir, 'figS1_pca')


def figS2_core_clock_genes(clock_gene_df, metadata_df, output_dir):
    """Supplementary Figure S2: Core clock gene expression."""
    fig, ax = plt.subplots(figsize=(14, 8))

    # Heatmap of clock gene expression, flight vs ground
    merged = clock_gene_df.merge(metadata_df[['sample_name', 'condition', 'osd_id']], on=['sample_name', 'osd_id'])

    # Pivot: genes x samples
    pivot = merged.pivot_table(index='gene', columns='sample_name', values='expression')
    # Sort columns by condition then study
    col_meta = merged[['sample_name', 'condition', 'osd_id']].drop_duplicates().set_index('sample_name')
    pivot = pivot[col_meta.sort_values(['condition', 'osd_id']).index]

    im = ax.imshow(pivot.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_xticks([])
    ax.set_xlabel('Samples (sorted by condition, then study)')
    ax.set_title('Core clock gene expression (z-scored)')

    plt.colorbar(im, ax=ax, label='Z-score')
    plt.tight_layout()
    return save_figure(fig, output_dir, 'figS2_clock_genes')


def figS3_phase_by_hardware(per_study_df, output_dir):
    """Supplementary Figure S3: Phase shift by hardware type."""
    included = per_study_df[per_study_df['included_in_meta']].copy()
    hardware_groups = included.groupby('hardware')

    fig, ax = plt.subplots(figsize=(12, 7))

    positions = []
    labels = []
    data_list = []
    colors_list = []

    for i, (hw, group) in enumerate(hardware_groups):
        if len(group) == 0:
            continue
        positions.append(i)
        # Shorten long hardware names
        short_label = hw
        if len(hw) > 30:
            short_label = hw[:28] + '...'
        labels.append(short_label)
        data_list.append(group['phase_shift_hours'].values)
        colors_list.append(CB_PALETTE[i % len(CB_PALETTE)])

    if data_list:
        bp = ax.boxplot(data_list, positions=positions, patch_artist=True, widths=0.5)
        for box, color in zip(bp['boxes'], colors_list):
            box.set_facecolor(color)
            box.set_alpha(0.7)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('Phase shift (hours)')
    ax.set_title('Circadian phase shift by hardware type')

    plt.subplots_adjust(bottom=0.3)
    plt.tight_layout()
    return save_figure(fig, output_dir, 'figS3_phase_by_hardware')


def figS4_phase_by_light_regime(per_study_df, output_dir):
    """Supplementary Figure S4: Phase shift by light regime."""
    included = per_study_df[per_study_df['included_in_meta']].copy()

    fig, ax = plt.subplots(figsize=(10, 7))
    light_groups = included.groupby('light_regime')

    positions = []
    labels = []
    data_list = []

    for i, (lr, group) in enumerate(light_groups):
        if len(group) == 0:
            continue
        positions.append(i)
        labels.append(lr)
        data_list.append(group['phase_shift_hours'].values)

    if data_list:
        bp = ax.boxplot(data_list, positions=positions, patch_artist=True, widths=0.5)
        for j, box in enumerate(bp['boxes']):
            box.set_facecolor(CB_PALETTE[j % len(CB_PALETTE)])
            box.set_alpha(0.7)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_ylabel('Phase shift (hours)')
    ax.set_title('Circadian phase shift by light regime')

    plt.subplots_adjust(bottom=0.25)
    plt.tight_layout()
    return save_figure(fig, output_dir, 'figS4_phase_by_light_regime')


def figS5_fork_comparison(fork_a_results, fork_b_results, output_dir):
    """Supplementary Figure S5: Fork A vs Fork B comparison.
    Uses two-column layout (text + plot) like fig4.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    for idx, (fork_name, results) in enumerate([('Fork A (Barker 2023)', fork_a_results),
                                                  ('Fork B (All Arabidopsis)', fork_b_results)]):
        ax = axes[idx]
        included = results[results['included_in_meta']].sort_values('phase_shift_hours')
        n = len(included)

        # Use gridspec within each subplot for text + plot
        gs = ax.get_subplotspec().subgridspec(1, 2, width_ratios=[1, 1.5], wspace=0.02)
        ax.remove()
        ax_text = fig.add_subplot(gs[0])
        ax_plot = fig.add_subplot(gs[1])

        for i, (_, row) in enumerate(included.iterrows()):
            color = PHYLO_COLORS['orange'] if row['p_value'] < 0.05 else PHYLO_COLORS['blue']
            ax_plot.errorbar(row['phase_shift_hours'], i,
                        xerr=[[row['phase_shift_hours'] - row['phase_shift_ci_lower']],
                              [row['phase_shift_ci_upper'] - row['phase_shift_hours']]],
                        fmt='o', color=color, capsize=3, markersize=5)

        ax_plot.axvline(0, color='k', linestyle='--', alpha=0.3)
        ax_plot.set_xlabel('Phase shift (h)')
        ax_plot.set_title(fork_name, fontweight='bold')
        ax_plot.set_yticks([])
        ax_plot.set_ylim(-0.5, n - 0.5)

        ax_text.set_xlim(0, 1)
        ax_text.set_ylim(-0.5, n - 0.5)
        ax_text.axis('off')
        for i, (_, row) in enumerate(included.iterrows()):
            ax_text.text(0.98, i, row['osd_id'], ha='right', va='center', fontsize=7, family='monospace')

    fig.suptitle('Fork A vs Fork B: per-study phase shifts', fontsize=14, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    return save_figure(fig, output_dir, 'figS5_fork_comparison')


def figS6_model_uncertainty(per_study_df, predictions_df, output_dir):
    """Supplementary Figure S6: Model uncertainty analysis.
    Panel A now computes mean circular variance per study from predictions_df.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: Circular variance vs phase shift magnitude
    ax = axes[0]
    included = per_study_df[per_study_df['included_in_meta']].copy()

    # Compute mean circular variance per study from predictions_df
    study_var = predictions_df.groupby('osd_id')['circular_variance'].mean()
    included['mean_circular_variance'] = included['osd_id'].map(study_var)

    ax.scatter(included['phase_shift_hours'].abs(), included['mean_circular_variance'],
               c=PHYLO_COLORS['blue'], alpha=0.6, s=50, edgecolors='white', linewidth=0.3)
    ax.set_xlabel('|Phase shift| (h)')
    ax.set_ylabel('Mean circular variance')
    ax.set_title('A', loc='left', fontweight='bold', fontsize=14)

    # Panel B: Number of models used per study
    ax = axes[1]
    if 'n_models_used' in predictions_df.columns:
        model_counts = predictions_df.groupby('osd_id')['n_models_used'].first()
        ax.bar(range(len(model_counts)), model_counts.values, color=PHYLO_COLORS['green'], alpha=0.7)
        ax.set_xticks(range(len(model_counts)))
        ax.set_xticklabels(model_counts.index, rotation=45, ha='right', fontsize=7)
        ax.set_ylabel('Models used (of 100)')
        ax.set_title('B', loc='left', fontweight='bold', fontsize=14)

    plt.tight_layout()
    return save_figure(fig, output_dir, 'figS6_model_uncertainty')


if __name__ == "__main__":
    print("Figure generation module loaded successfully.")
    print("This module is called by the main analysis pipeline.")
    print("Available functions:")
    for name in dir():
        if name.startswith('fig'):
            print(f"  - {name}")
