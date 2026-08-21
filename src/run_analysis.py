"""
Main analysis pipeline.
Runs ChronoGauge on all GeneLab studies, performs statistical analysis,
and generates all figures and tables.

Usage:
    python3 src/run_analysis.py
"""
import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

# Suppress warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
warnings.filterwarnings('ignore')

# Add src to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tables'))

from chronogauge_apply import predict_ct, circular_error, load_ensemble, load_training_scaler
from metadata_curation import curate_all_metadata
from statistical_analysis import per_study_analysis, meta_analysis, circular_mean_hours
from generate_figures import (
    fig1_study_overview, fig2_chronogauge_validation, fig3_phase_shift_polar,
    fig4_forest_plot, fig5_stratified_analysis, fig6_circadian_fingerprint,
    figS1_pca_umap, figS2_core_clock_genes, figS3_phase_by_hardware,
    figS4_phase_by_light_regime, figS5_fork_comparison, figS6_model_uncertainty
)
from generate_tables import (
    table1_dataset_characteristics, table2_per_study_results,
    table3_meta_analysis, tableS1_sample_metadata, tableS2_chronogauge_validation
)

# Paths - detect if running in Docker or locally
if os.path.exists("/workspace/genelab_data"):
    DATA_DIR = "/workspace/genelab_data"
else:
    DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "genelab_data"))

if os.path.exists("/mnt/results"):
    RESULTS_DIR = "/mnt/results"
elif os.path.exists("/results"):
    RESULTS_DIR = "/results"
else:
    RESULTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
TABLES_DIR = os.path.join(RESULTS_DIR, "tables")
DATA_OUT_DIR = os.path.join(RESULTS_DIR, "data")

# Fork definitions
FORK_A_IDS = [7, 17, 37, 38, 44, 46, 120, 121, 136, 147, 205, 208, 218, 251]
FORK_B_FILE = os.path.join(DATA_DIR, "fork_b_osd_ids.json")

# Core clock genes for supplementary analysis
CLOCK_GENES = {
    'morning': ['AT2G46830', 'AT1G01060'],  # CCA1, LHY
    'morning_pr': ['AT2G42830', 'AT5G02810', 'AT5G24470'],  # PRR9, PRR7, PRR5
    'evening': ['AT5G61380', 'AT3G09600'],  # TOC1, GI
    'evening_ec': ['AT2G40080', 'AT3G26740', 'AT1G75530', 'AT3G46640'],  # ELF4, ELF3, LUX, TIC
}
CLOCK_GENE_NAMES = {
    'AT2G46830': 'CCA1', 'AT1G01060': 'LHY',
    'AT2G42830': 'PRR9', 'AT5G02810': 'PRR7', 'AT5G24470': 'PRR5',
    'AT5G61380': 'TOC1', 'AT3G09600': 'GI',
    'AT2G40080': 'ELF4', 'AT3G26740': 'ELF3', 'AT1G75530': 'LUX', 'AT3G46640': 'TIC',
}


def find_expression_file(study_dir, osd_id):
    """Find the normalized expression file for a study."""
    files = os.listdir(study_dir)

    # RNA-seq: Normalized_Counts (prefer rRNArm version as it's cleaner)
    rna_files = [f for f in files if 'Normalized_Counts' in f and f.endswith('.csv')]
    rna_rnarm = [f for f in rna_files if 'rRNArm' in f]
    if rna_rnarm:
        return os.path.join(study_dir, rna_rnarm[0])
    if rna_files:
        return os.path.join(study_dir, rna_files[0])

    # Microarray: normalized_expression_probeset
    array_files = [f for f in files if 'normalized_expression_probeset' in f and f.endswith('.csv')]
    if array_files:
        return os.path.join(study_dir, array_files[0])

    return None


def load_expression_matrix(file_path, platform='rnaseq'):
    """Load expression matrix and return as genes x samples DataFrame."""
    df = pd.read_csv(file_path, index_col=0)

    if platform == 'rnaseq':
        # RNA-seq: genes are already AGI codes in index, samples are columns
        return df
    else:
        # Microarray: may have annotation columns before sample columns
        # Sample columns typically start with 'Atha_' or contain 'FLT'/'GC'
        sample_cols = [c for c in df.columns if 'Atha_' in c or 'FLT' in c or 'GC' in c or '_Rep' in c]
        if not sample_cols:
            # Try to identify sample columns (non-annotation)
            annotation_cols = ['SYMBOL', 'GENENAME', 'REFSEQ', 'ENTREZID', 'STRING_id',
                               'GOSLIM_IDS', 'ProbesetID', 'count_ENSEMBL_mappings']
            sample_cols = [c for c in df.columns if c not in annotation_cols]

        # Index should be AGI codes
        return df[sample_cols]


def determine_platform(study_dir):
    """Determine if a study is RNA-seq or microarray."""
    files = os.listdir(study_dir)
    if any('rna_seq' in f for f in files):
        return 'rnaseq'
    elif any('array' in f or 'microarray' in f for f in files):
        # Determine ATH1 vs AraGene
        # ATH1: Affymetrix ATH1-121501 array (22K probes)
        # AraGene: Affymetrix Arabidopsis Gene 1.0 ST array (28K probes)
        # Check probe count in expression file
        array_file = [f for f in files if 'normalized_expression_probeset' in f]
        if array_file:
            df = pd.read_csv(os.path.join(study_dir, array_file[0]), index_col=0, nrows=1)
            n_rows = sum(1 for _ in open(os.path.join(study_dir, array_file[0]))) - 1
            if n_rows > 25000:
                return 'aragene'
            else:
                return 'ath1'
        return 'ath1'  # Default
    return 'rnaseq'


def run_chronogauge_for_study(study_dir, osd_id, platform, n_models=100):
    """Run ChronoGauge on a single study."""
    expr_file = find_expression_file(study_dir, osd_id)
    if not expr_file:
        print(f"  No expression file found for OSD-{osd_id}")
        return None, None

    print(f"  Loading expression data from {os.path.basename(expr_file)}")
    expr = load_expression_matrix(expr_file, platform)

    # Clean gene IDs (remove version suffixes if any)
    expr.index = expr.index.astype(str).str.replace(r'\.\d+$', '', regex=True)

    # Filter to AGI codes only
    agi_mask = expr.index.str.match(r'^AT[1-5]G\d{5}$', na=False)
    if agi_mask.sum() < 100:
        print(f"  WARNING: Only {agi_mask.sum()} AGI genes found for OSD-{osd_id}")
    expr = expr[agi_mask]

    # Remove duplicate genes (keep first)
    expr = expr[~expr.index.duplicated(keep='first')]

    print(f"  Expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

    try:
        results, all_preds = predict_ct(expr, platform=platform, n_models=n_models)
        results['osd_id'] = f'OSD-{osd_id}'
        results['sample_name'] = results.index
        results = results.reset_index(drop=True)

        all_preds['osd_id'] = f'OSD-{osd_id}'
        all_preds['sample_name'] = all_preds.index
        all_preds = all_preds.reset_index(drop=True)

        print(f"  Predicted CT for {len(results)} samples (mean variance: {results['circular_variance'].mean():.4f})")
        return results, all_preds
    except Exception as e:
        print(f"  ERROR running ChronoGauge for OSD-{osd_id}: {e}")
        return None, None


def run_validation(metadata_df, predictions_df):
    """Validate ChronoGauge on ground control samples where ZT is known."""
    # Most GeneLab studies don't have explicit ZT, so validation is limited
    # We can validate on ChronoGauge's own test data instead
    return None


def extract_clock_gene_expression(study_dir, osd_id, platform):
    """Extract core clock gene expression for a study."""
    expr_file = find_expression_file(study_dir, osd_id)
    if not expr_file:
        return None

    expr = load_expression_matrix(expr_file, platform)
    expr.index = expr.index.astype(str).str.replace(r'\.\d+$', '', regex=True)

    # Filter to clock genes
    all_clock_genes = []
    for genes in CLOCK_GENES.values():
        all_clock_genes.extend(genes)

    clock_expr = expr[expr.index.isin(all_clock_genes)]

    if len(clock_expr) == 0:
        return None

    # Z-score per gene
    clock_expr_z = clock_expr.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

    # Convert to long format - reset_index and use the actual column name
    clock_expr_z = clock_expr_z.reset_index()
    gene_col = clock_expr_z.columns[0]  # First column is the gene index
    clock_long = clock_expr_z.melt(id_vars=[gene_col], var_name='sample_name', value_name='expression')
    clock_long = clock_long.rename(columns={gene_col: 'gene'})
    clock_long['gene_name'] = clock_long['gene'].map(CLOCK_GENE_NAMES)
    clock_long['osd_id'] = f'OSD-{osd_id}'

    return clock_long


def main():
    print("=" * 70)
    print("SPACEFLIGHT CIRCADIAN DECODER - MAIN ANALYSIS PIPELINE")
    print("=" * 70)

    # Create output directories
    for d in [FIGURES_DIR, TABLES_DIR, DATA_OUT_DIR]:
        os.makedirs(d, exist_ok=True)

    # Load Fork B IDs
    with open(FORK_B_FILE) as f:
        fork_b_ids = [int(x.split('-')[1]) for x in json.load(f)]

    all_ids = sorted(set(FORK_A_IDS + fork_b_ids))
    print(f"\nTotal studies to analyze: {len(all_ids)}")
    print(f"Fork A: {len(FORK_A_IDS)} studies")
    print(f"Fork B: {len(fork_b_ids)} studies")

    # Step 1: Curate metadata
    print("\n" + "=" * 70)
    print("STEP 1: METADATA CURATION")
    print("=" * 70)
    
    # Load existing metadata cache if available
    metadata_cache_path = os.path.join(TABLES_DIR, "tableS1_sample_metadata.csv")
    if os.path.exists(metadata_cache_path):
        print(f"Loading cached sample metadata from {os.path.basename(metadata_cache_path)}...")
        cached_metadata = pd.read_csv(metadata_cache_path)
    else:
        cached_metadata = pd.DataFrame()
        
    # Curate new metadata for studies in DATA_DIR
    new_metadata = curate_all_metadata(DATA_DIR, os.path.join(DATA_DIR, "arabidopsis_transcriptomics.json"))
    
    if len(cached_metadata) > 0 and len(new_metadata) > 0:
        new_study_ids = new_metadata['osd_id'].unique().tolist()
        print(f"Merging new metadata for studies {new_study_ids} with cache...")
        # Remove updated studies from cache
        cached_filtered = cached_metadata[~cached_metadata['osd_id'].isin(new_study_ids)]
        metadata = pd.concat([cached_filtered, new_metadata], ignore_index=True)
    elif len(new_metadata) > 0:
        metadata = new_metadata
    else:
        metadata = cached_metadata

    print(f"Total samples: {len(metadata)}")

    # Step 2: Run ChronoGauge on all studies
    print("\n" + "=" * 70)
    print("STEP 2: CHRONOGAUGE APPLICATION")
    print("=" * 70)

    # Check for cached predictions
    cache_path = os.path.join(DATA_OUT_DIR, 'all_predictions.csv')
    cache_model_path = os.path.join(DATA_OUT_DIR, 'all_model_predictions.csv')
    cache_clock_path = os.path.join(DATA_OUT_DIR, 'clock_gene_expression.csv')

    if os.path.exists(cache_path) and os.path.exists(cache_model_path):
        print("Loading cached predictions...")
        predictions_df = pd.read_csv(cache_path)
        model_preds_df = pd.read_csv(cache_model_path)
        if os.path.exists(cache_clock_path):
            clock_expr_df = pd.read_csv(cache_clock_path)
        else:
            clock_expr_df = pd.DataFrame()
        print(f"Loaded {len(predictions_df)} predictions from {predictions_df['osd_id'].nunique()} studies")
    else:
        all_predictions = []
        all_model_preds = []
        all_clock_expr = []

        for osd_id in all_ids:
            study_dir = os.path.join(DATA_DIR, f"OSD-{osd_id}")
            if not os.path.isdir(study_dir):
                print(f"\nOSD-{osd_id}: Directory not found, skipping")
                continue

            # Check if expression file exists
            expr_file = find_expression_file(study_dir, osd_id)
            if not expr_file:
                print(f"\nOSD-{osd_id}: No expression file, skipping")
                continue

            platform = determine_platform(study_dir)
            print(f"\nOSD-{osd_id} (platform={platform}):")

            results, model_preds = run_chronogauge_for_study(study_dir, osd_id, platform, n_models=100)
            if results is not None:
                all_predictions.append(results)
                all_model_preds.append(model_preds)

            # Extract clock gene expression
            clock_expr = extract_clock_gene_expression(study_dir, osd_id, platform)
            if clock_expr is not None:
                all_clock_expr.append(clock_expr)

        # Combine all predictions
        predictions_df = pd.concat(all_predictions, ignore_index=True)
        model_preds_df = pd.concat(all_model_preds, ignore_index=True)
        clock_expr_df = pd.concat(all_clock_expr, ignore_index=True) if all_clock_expr else pd.DataFrame()

        print(f"\nTotal predictions: {len(predictions_df)}")
        print(f"Studies with predictions: {predictions_df['osd_id'].nunique()}")

        # Save predictions
        predictions_df.to_csv(os.path.join(DATA_OUT_DIR, 'all_predictions.csv'), index=False)
        model_preds_df.to_csv(os.path.join(DATA_OUT_DIR, 'all_model_predictions.csv'), index=False)
        if len(clock_expr_df) > 0:
            clock_expr_df.to_csv(os.path.join(DATA_OUT_DIR, 'clock_gene_expression.csv'), index=False)

    # Step 3: Statistical analysis
    print("\n" + "=" * 70)
    print("STEP 3: STATISTICAL ANALYSIS")
    print("=" * 70)

    # Merge predictions with metadata
    merged = predictions_df.merge(metadata, on=['sample_name', 'osd_id'], how='left', suffixes=('', '_meta'))

    # Per-study analysis
    per_study = per_study_analysis(predictions_df, metadata)
    print(f"\nPer-study results: {len(per_study)} studies")
    print(per_study[['osd_id', 'n_flight', 'n_ground', 'phase_shift_hours', 'p_value']].to_string())

    # Meta-analysis - overall
    print("\nOverall meta-analysis:")
    meta_overall = meta_analysis(per_study)
    print(json.dumps(meta_overall, indent=2, default=str))

    # Meta-analysis - stratified by tissue
    print("\nStratified by tissue:")
    meta_tissue = meta_analysis(per_study, stratify_by='tissue')
    print(json.dumps(meta_tissue, indent=2, default=str))

    # Meta-analysis - stratified by hardware
    print("\nStratified by hardware:")
    meta_hardware = meta_analysis(per_study, stratify_by='hardware')
    print(json.dumps(meta_hardware, indent=2, default=str))

    # Meta-analysis - stratified by light regime
    print("\nStratified by light regime:")
    meta_light = meta_analysis(per_study, stratify_by='light_regime')
    print(json.dumps(meta_light, indent=2, default=str))

    # Step 4: Fork A vs Fork B analysis
    print("\n" + "=" * 70)
    print("STEP 4: FORK A vs FORK B ANALYSIS")
    print("=" * 70)

    fork_a_set = {f'OSD-{x}' for x in FORK_A_IDS}
    fork_b_set = {f'OSD-{x}' for x in fork_b_ids}

    fork_a_results = per_study[per_study['osd_id'].isin(fork_a_set)]
    fork_b_results = per_study[per_study['osd_id'].isin(fork_b_set)]

    print(f"Fork A: {len(fork_a_results)} studies, {fork_a_results['included_in_meta'].sum()} in meta-analysis")
    print(f"Fork B: {len(fork_b_results)} studies, {fork_b_results['included_in_meta'].sum()} in meta-analysis")

    fork_a_meta = meta_analysis(fork_a_results)
    fork_b_meta = meta_analysis(fork_b_results)
    print(f"\nFork A meta-analysis: {json.dumps(fork_a_meta, indent=2, default=str)}")
    print(f"\nFork B meta-analysis: {json.dumps(fork_b_meta, indent=2, default=str)}")

    # Step 5: Generate figures
    print("\n" + "=" * 70)
    print("STEP 5: FIGURE GENERATION")
    print("=" * 70)

    print("\nMain figures:")
    fig1_study_overview(metadata, per_study, FIGURES_DIR)
    fig3_phase_shift_polar(per_study, FIGURES_DIR)
    fig4_forest_plot(per_study, meta_overall, FIGURES_DIR)
    fig5_stratified_analysis(per_study, meta_tissue, FIGURES_DIR, stratify_by='tissue')
    fig6_circadian_fingerprint(predictions_df, metadata, FIGURES_DIR)

    print("\nSupplementary figures:")
    figS3_phase_by_hardware(per_study, FIGURES_DIR)
    figS4_phase_by_light_regime(per_study, FIGURES_DIR)
    figS5_fork_comparison(fork_a_results, fork_b_results, FIGURES_DIR)
    figS6_model_uncertainty(per_study, predictions_df, FIGURES_DIR)

    # Clock gene heatmap
    if len(clock_expr_df) > 0:
        figS2_core_clock_genes(clock_expr_df, metadata, FIGURES_DIR)

    # Step 6: Generate tables
    print("\n" + "=" * 70)
    print("STEP 6: TABLE GENERATION")
    print("=" * 70)

    table1_dataset_characteristics(metadata, TABLES_DIR)
    table2_per_study_results(per_study, TABLES_DIR)
    table3_meta_analysis(meta_overall, TABLES_DIR, 'overall')
    table3_meta_analysis(meta_tissue, TABLES_DIR, 'tissue')
    table3_meta_analysis(meta_hardware, TABLES_DIR, 'hardware')
    table3_meta_analysis(meta_light, TABLES_DIR, 'light_regime')
    tableS1_sample_metadata(metadata, TABLES_DIR)

    # Save all results
    per_study.to_csv(os.path.join(DATA_OUT_DIR, 'per_study_results.csv'), index=False)

    # Save meta-analysis results as JSON
    with open(os.path.join(DATA_OUT_DIR, 'meta_analysis_results.json'), 'w') as f:
        json.dump({
            'overall': meta_overall,
            'tissue': meta_tissue,
            'hardware': meta_hardware,
            'light_regime': meta_light,
            'fork_a': fork_a_meta,
            'fork_b': fork_b_meta,
        }, f, indent=2, default=str)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Figures saved to: {FIGURES_DIR}")
    print(f"Tables saved to: {TABLES_DIR}")
    print(f"Data saved to: {DATA_OUT_DIR}")


if __name__ == "__main__":
    main()
