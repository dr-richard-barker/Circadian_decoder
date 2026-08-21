import os
import sys
import pandas as pd

# Add src to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from chronogauge_apply import predict_ct

def main():
    print("Generating ChronoGauge predictions for OSD-522...")
    
    # Load unnormalized counts
    counts_path = 'genelab_data/OSD-522/GLDS-522_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv'
    if not os.path.exists(counts_path):
        print(f"Error: Counts file not found at {counts_path}")
        return

    counts = pd.read_csv(counts_path, index_col=0)

    # CPM normalization
    cpm = counts.div(counts.sum(), axis=1) * 1e6

    # Filter to AGI codes only
    agi_mask = cpm.index.str.match(r'^AT[1-5]G\d{5}$', na=False)
    cpm = cpm[agi_mask]
    cpm = cpm[~cpm.index.duplicated(keep='first')]

    # Run ChronoGauge
    results, model_preds = predict_ct(cpm, platform='rnaseq', n_models=100)
    
    # Format results
    results['osd_id'] = 'OSD-522'
    results['sample_name'] = results.index
    results = results.reset_index(drop=True)

    model_preds['osd_id'] = 'OSD-522'
    model_preds['sample_name'] = model_preds.index
    model_preds = model_preds.reset_index(drop=True)

    # Cache paths
    cache_path = 'data/all_predictions.csv'
    cache_model_path = 'data/all_model_predictions.csv'

    # Append to predictions cache
    if os.path.exists(cache_path):
        preds_df = pd.read_csv(cache_path)
        # Remove any existing OSD-522
        preds_df = preds_df[preds_df['osd_id'] != 'OSD-522']
        preds_df = pd.concat([preds_df, results], ignore_index=True)
    else:
        preds_df = results
    
    # Append to model predictions cache
    if os.path.exists(cache_model_path):
        model_preds_df = pd.read_csv(cache_model_path)
        # Remove any existing OSD-522
        model_preds_df = model_preds_df[model_preds_df['osd_id'] != 'OSD-522']
        model_preds_df = pd.concat([model_preds_df, model_preds], ignore_index=True)
    else:
        model_preds_df = model_preds

    # Save back
    os.makedirs('data', exist_ok=True)
    preds_df.to_csv(cache_path, index=False)
    model_preds_df.to_csv(cache_model_path, index=False)
    print("Successfully integrated OSD-522 predictions into the cache!")

if __name__ == '__main__':
    main()
