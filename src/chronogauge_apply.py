"""
ChronoGauge application module.
Applies the ChronoGauge ensemble to predict circadian time (CT) from
gene expression matrices. Supports RNA-seq, ATH1 microarray, and AraGene microarray.
"""
import os
import warnings
import math
import json
import numpy as np
import pandas as pd
import tensorflow as tf

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
warnings.filterwarnings('ignore')

# ChronoGauge paths
src_dir = os.path.dirname(os.path.abspath(__file__))
if os.path.exists("/workspace/ChronoGauge"):
    CG_ROOT = "/workspace/ChronoGauge"
else:
    CG_ROOT = os.path.abspath(os.path.join(src_dir, "..", "ChronoGauge"))

CG_DATA = os.path.join(CG_ROOT, "data")
CG_MODELS = {
    'rnaseq': os.path.join(CG_ROOT, 'models/hf_rna'),
    'ath1': os.path.join(CG_ROOT, 'models/hf_ath1'),
    'aragene': os.path.join(CG_ROOT, 'models/hf_aragene'),
}

# Training expression matrices for scaler fitting
CG_TRAINING = {
    'rnaseq': os.path.join(CG_DATA, 'expression_matrices/x_training.csv'),
    'ath1': os.path.join(CG_DATA, 'expression_matrices/x_test_ath.csv'),  # ATH1 training data
    'aragene': os.path.join(CG_DATA, 'expression_matrices/x_test_aragene.csv'),
}

# Feature gene sets
CG_FEATURES = os.path.join(CG_DATA, 'model_parameters/gene_features_unadjusted.csv')


def time24(ipreds):
    """Convert circular (cos, sin) predictions to 24-hour CT values."""
    preds = []
    for i in range(ipreds.shape[0]):
        preds.append(math.atan2(ipreds[i, 0], ipreds[i, 1]) / math.pi * 12)
    for i in range(len(preds)):
        if preds[i] < 0:
            preds[i] = preds[i] + 24
    return preds


def cyclic_time(times):
    """Convert 24h times to circular (cos, sin) representation."""
    times = np.asarray(times) % 24
    t_cos = -np.cos((2 * np.pi * times.astype('float64') / 24) + (np.pi / 2))
    t_sin = np.sin((2 * np.pi * times.astype('float64') / 24) + (np.pi / 2))
    return t_cos, t_sin


def circular_mean(predictions_24):
    """Aggregate multiple CT predictions via circular mean."""
    cos_vals = []
    sin_vals = []
    for i in range(predictions_24.shape[1]):
        c, s = cyclic_time(predictions_24.iloc[:, i])
        cos_vals.append(c)
        sin_vals.append(s)
    cos_mean = np.mean(cos_vals, axis=0)
    sin_mean = np.mean(sin_vals, axis=0)
    ct_vals = np.concatenate(
        (np.asarray(cos_mean).reshape(-1, 1), np.asarray(sin_mean).reshape(-1, 1)),
        axis=1
    )
    return time24(ct_vals)


def circular_variance(predictions_24):
    """Compute circular variance across sub-model predictions (clock robustness metric)."""
    cos_vals = []
    sin_vals = []
    for i in range(predictions_24.shape[1]):
        c, s = cyclic_time(predictions_24.iloc[:, i])
        cos_vals.append(c)
        sin_vals.append(s)
    cos_mean = np.mean(cos_vals, axis=0)
    sin_mean = np.mean(sin_vals, axis=0)
    # Resultant length R: 0=max dispersion, 1=max concentration
    R = np.sqrt(cos_mean**2 + sin_mean**2)
    # Circular variance: V = 1 - R
    return 1 - R


def load_training_scaler(platform='rnaseq'):
    """Load training expression matrix and fit StandardScaler."""
    train_path = CG_TRAINING.get(platform)
    if not train_path or not os.path.exists(train_path):
        raise FileNotFoundError(f"Training data not found for platform {platform}: {train_path}")

    X_train = pd.read_csv(train_path, index_col=0)
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    scaler.fit(X_train.T)  # fit on samples (columns), transform genes (rows)
    return scaler, X_train.index


def load_ensemble(platform='rnaseq', n_models=100):
    """Load the ChronoGauge ensemble models and their feature gene sets."""
    model_dir = CG_MODELS.get(platform)
    if not model_dir or not os.path.exists(model_dir):
        raise FileNotFoundError(f"Model directory not found for platform {platform}: {model_dir}")

    features_df = pd.read_csv(CG_FEATURES, index_col=0)

    ensemble = {}
    for i in range(n_models):
        model_path = os.path.join(model_dir, f'model_{i}.h5')
        if not os.path.exists(model_path):
            continue
        model = tf.keras.models.load_model(model_path, compile=False)
        gene_features = features_df.iloc[i].dropna().to_numpy()
        # Truncate features to match model input size (CSV may have extra features)
        input_dim = model.input_shape[-1]
        if input_dim is not None and len(gene_features) > input_dim:
            gene_features = gene_features[:input_dim]
        ensemble[i] = (model, gene_features)

    print(f"Loaded {len(ensemble)} models for platform '{platform}'")
    return ensemble


def predict_ct(expression_matrix, platform='rnaseq', n_models=100):
    """
    Predict circadian time for each sample in an expression matrix.

    Parameters:
    -----------
    expression_matrix : pd.DataFrame
        Gene x Sample matrix (genes as rows, samples as columns).
        Gene IDs must be AGI codes (AT*G*).
    platform : str
        'rnaseq', 'ath1', or 'aragene'
    n_models : int
        Number of ensemble sub-models to use (max 100)

    Returns:
    --------
    pd.DataFrame with columns: predicted_CT, circular_variance, n_models_used
    """
    # Load ensemble and scaler
    ensemble = load_ensemble(platform, n_models)
    scaler, train_genes = load_training_scaler(platform)

    # Align genes with training data
    common_genes = expression_matrix.index.intersection(train_genes)
    if len(common_genes) < len(train_genes) * 0.5:
        print(f"WARNING: Only {len(common_genes)}/{len(train_genes)} training genes found in expression matrix")
    X = expression_matrix.loc[expression_matrix.index.intersection(train_genes)]
    # Reindex to match training gene order
    X = X.reindex(train_genes)
    # Fill missing genes with 0
    X = X.fillna(0)

    # Scale: transform genes (rows) using scaler fit on samples
    X_scaled = pd.DataFrame(
        data=scaler.transform(X.T).T,
        index=X.index,
        columns=X.columns
    )

    # Run each sub-model
    all_preds = {}
    for model_id, (model, gene_features) in ensemble.items():
        # Check that all feature genes are present
        available_features = [g for g in gene_features if g in X_scaled.index]
        if len(available_features) < len(gene_features) * 0.8:
            continue  # Skip model if too many features missing

        # Use reindex to handle missing genes (fill with 0)
        X_sub = X_scaled.reindex(gene_features).fillna(0)
        results = model(X_sub.T)
        preds = time24(np.asarray(results))
        all_preds[model_id] = preds

    if not all_preds:
        raise ValueError("No models could run - insufficient gene overlap")

    # Create predictions DataFrame (samples x models)
    preds_df = pd.DataFrame(all_preds, index=expression_matrix.columns)

    # Aggregate via circular mean
    ct_pred = circular_mean(preds_df)
    ct_var = circular_variance(preds_df)

    results_df = pd.DataFrame({
        'predicted_CT': ct_pred,
        'circular_variance': ct_var,
        'n_models_used': len(all_preds)
    }, index=expression_matrix.columns)

    # Also store individual model predictions for fingerprint analysis
    return results_df, preds_df


def circular_error(pred, true):
    """Compute circular error in minutes."""
    true = np.asarray(true) % 24
    err = np.asarray(pred) - true
    for i in range(len(err)):
        if err[i] > 12:
            err[i] -= 24
        if err[i] < -12:
            err[i] += 24
    return err * 60


if __name__ == "__main__":
    # Quick test on ChronoGauge's own test data
    print("Testing ChronoGauge pipeline on test RNA-seq data...")
    X_test = pd.read_csv(os.path.join(CG_DATA, 'expression_matrices/x_test_rna.csv'), index_col=0)
    Y_test = pd.read_csv(os.path.join(CG_DATA, 'targets/target_test_rna.csv'), index_col=0)

    results, all_preds = predict_ct(X_test, platform='rnaseq', n_models=20)
    true_times = Y_test.iloc[:, 0] % 24
    errors = circular_error(results['predicted_CT'].values, true_times.values)
    mae = np.mean(np.abs(errors))
    print(f"MAE (20 models): {mae:.1f} minutes")
    print(f"Mean circular variance: {results['circular_variance'].mean():.4f}")
    print("\nSample predictions:")
    for idx in range(10):
        print(f"  {Y_test.index[idx]:12s}  true={true_times.iloc[idx]:5.0f}h  pred={results['predicted_CT'].iloc[idx]:5.2f}h  err={errors[idx]:6.1f}min")
