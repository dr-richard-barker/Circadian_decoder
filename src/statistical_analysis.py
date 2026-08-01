"""
Statistical analysis module.
Performs per-study flight-vs-ground circadian phase comparison,
then random-effects meta-analysis across studies.
Stratified by tissue, genotype, and hardware.
"""
import os
import json
import numpy as np
import pandas as pd
import subprocess
import tempfile
from scipy import stats


def circular_distance(a, b):
    """Compute circular distance between two angles (in hours, 0-24)."""
    diff = (np.asarray(a) - np.asarray(b)) % 24
    diff = np.where(diff > 12, 24 - diff, diff)
    return diff


def circular_mean_hours(angles):
    """Compute circular mean of angles in hours (0-24)."""
    angles = np.asarray(angles, dtype=float) % 24
    cos_vals = np.cos(2 * np.pi * angles / 24)
    sin_vals = np.sin(2 * np.pi * angles / 24)
    mean_angle = np.arctan2(np.mean(sin_vals), np.mean(cos_vals))
    mean_hours = (mean_angle * 12 / np.pi) % 24
    return mean_hours


def circular_std_hours(angles):
    """Compute circular standard deviation in hours."""
    angles = np.asarray(angles, dtype=float) % 24
    cos_vals = np.cos(2 * np.pi * angles / 24)
    sin_vals = np.sin(2 * np.pi * angles / 24)
    R = np.sqrt(np.mean(cos_vals)**2 + np.mean(sin_vals)**2)
    if R == 0:
        return 12.0
    circ_std = np.sqrt(-2 * np.log(R)) * 24 / (2 * np.pi)
    return circ_std


def watson_williams_test(group1, group2):
    """
    Watson-Williams test for equality of circular means.
    Returns (F_statistic, p_value).
    """
    g1 = np.asarray(group1, dtype=float) % 24
    g2 = np.asarray(group2, dtype=float) % 24
    n1, n2 = len(g1), len(g2)
    n = n1 + n2

    if n1 < 2 or n2 < 2:
        return np.nan, np.nan

    # Resultant lengths
    R1 = np.sqrt(np.sum(np.cos(2*np.pi*g1/24))**2 + np.sum(np.sin(2*np.pi*g1/24))**2) / n1
    R2 = np.sqrt(np.sum(np.cos(2*np.pi*g2/24))**2 + np.sum(np.sin(2*np.pi*g2/24))**2) / n2

    # Combined
    all_angles = np.concatenate([g1, g2])
    R_all = np.sqrt(np.sum(np.cos(2*np.pi*all_angles/24))**2 + np.sum(np.sin(2*np.pi*all_angles/24))**2) / n

    # Watson-Williams F statistic with correction factor
    if R_all >= 0.95:
        k = 1 + 3/(8 * 2 * R_all)
    else:
        k = 1

    if (1 - R_all) == 0:
        return np.nan, np.nan

    F = k * (n - 2) * (n1 * R1 + n2 * R2 - n * R_all) / ((n - 2) * (1 - R_all) + 0.001)
    p_value = 1 - stats.f.cdf(F, 1, n - 2)

    return F, p_value


def per_study_analysis(predictions_df, metadata_df, min_replicates=3):
    """
    Perform per-study flight-vs-ground comparison of predicted CT.
    """
    # Merge predictions with metadata on both sample_name and osd_id
    # to avoid column conflicts
    merge_cols = [c for c in ['sample_name', 'osd_id'] if c in predictions_df.columns and c in metadata_df.columns]
    merged = predictions_df.merge(metadata_df, on=merge_cols, how='inner')

    results = []

    for osd_id, study_data in merged.groupby('osd_id'):
        flight = study_data[study_data['condition'] == 'flight']['predicted_CT'].values
        ground = study_data[study_data['condition'] == 'ground_control']['predicted_CT'].values

        # Get study-level metadata (mode of each field)
        def safe_mode(s):
            s = s.dropna()
            s = s[s != '']
            return s.mode().iloc[0] if len(s) > 0 else ''

        tissue = safe_mode(study_data['tissue']) if 'tissue' in study_data.columns else ''
        ecotype = safe_mode(study_data['ecotype']) if 'ecotype' in study_data.columns else ''
        hardware = safe_mode(study_data['hardware']) if 'hardware' in study_data.columns else ''
        light_regime = safe_mode(study_data['light_regime']) if 'light_regime' in study_data.columns else ''

        if len(flight) < min_replicates or len(ground) < min_replicates:
            results.append({
                'osd_id': osd_id,
                'n_flight': len(flight),
                'n_ground': len(ground),
                'flight_mean_CT': circular_mean_hours(flight) if len(flight) > 0 else np.nan,
                'ground_mean_CT': circular_mean_hours(ground) if len(ground) > 0 else np.nan,
                'phase_shift_hours': np.nan,
                'phase_shift_ci_lower': np.nan,
                'phase_shift_ci_upper': np.nan,
                'p_value': np.nan,
                'F_statistic': np.nan,
                'tissue': tissue,
                'ecotype': ecotype,
                'hardware': hardware,
                'light_regime': light_regime,
                'included_in_meta': False,
                'reason': f'Insufficient replicates (flight={len(flight)}, ground={len(ground)})'
            })
            continue

        # Circular means
        flight_mean = circular_mean_hours(flight)
        ground_mean = circular_mean_hours(ground)

        # Phase shift (flight - ground), circular
        shift = (flight_mean - ground_mean) % 24
        if shift > 12:
            shift -= 24

        # Watson-Williams test
        F_stat, p_val = watson_williams_test(flight, ground)

        # Bootstrap CI for phase shift
        n_boot = 1000
        boot_shifts = []
        for _ in range(n_boot):
            f_boot = np.random.choice(flight, len(flight), replace=True)
            g_boot = np.random.choice(ground, len(ground), replace=True)
            f_mean = circular_mean_hours(f_boot)
            g_mean = circular_mean_hours(g_boot)
            s = (f_mean - g_mean) % 24
            if s > 12:
                s -= 24
            boot_shifts.append(s)

        ci_lower = np.percentile(boot_shifts, 2.5)
        ci_upper = np.percentile(boot_shifts, 97.5)

        results.append({
            'osd_id': osd_id,
            'n_flight': len(flight),
            'n_ground': len(ground),
            'flight_mean_CT': flight_mean,
            'ground_mean_CT': ground_mean,
            'phase_shift_hours': shift,
            'phase_shift_ci_lower': ci_lower,
            'phase_shift_ci_upper': ci_upper,
            'p_value': p_val,
            'F_statistic': F_stat,
            'tissue': tissue,
            'ecotype': ecotype,
            'hardware': hardware,
            'light_regime': light_regime,
            'included_in_meta': True,
            'reason': 'OK'
        })

    return pd.DataFrame(results)


def _run_metafor(data):
    """Run metafor random-effects meta-analysis via R subprocess."""
    csv_path = tempfile.mktemp(suffix='.csv')
    data[['osd_id', 'phase_shift_hours', 'se']].to_csv(csv_path, index=False)

    r_script = f"""
library(metafor)
data <- read.csv("{csv_path}")
res <- rma(yi=phase_shift_hours, sei=se, data=data, method="REML")
cat("POOLED_EFFECT:", res$b, "\\n")
cat("SE:", res$se, "\\n")
cat("CI_LOWER:", res$ci.lb, "\\n")
cat("CI_UPPER:", res$ci.ub, "\\n")
cat("Z:", res$zval, "\\n")
cat("P_VALUE:", res$pval, "\\n")
cat("Q:", res$QE, "\\n")
cat("Q_P:", res$QEp, "\\n")
cat("I2:", res$I2, "\\n")
cat("TAU2:", res$tau2, "\\n")
cat("N_STUDIES:", res$k, "\\n")
"""

    r_path = tempfile.mktemp(suffix='.R')
    with open(r_path, 'w') as f:
        f.write(r_script)

    try:
        result = subprocess.run(['Rscript', r_path], capture_output=True, text=True, timeout=60)
        output = result.stdout

        parsed = {}
        for line in output.split('\n'):
            if ':' in line and any(k in line for k in ['POOLED_EFFECT', 'SE:', 'CI_', 'Z:', 'P_VALUE', 'Q:', 'Q_P', 'I2', 'TAU2', 'N_STUDIES']):
                key, val = line.split(':', 1)
                try:
                    parsed[key.strip()] = float(val.strip())
                except ValueError:
                    parsed[key.strip()] = val.strip()

        os.unlink(csv_path)
        os.unlink(r_path)
        return parsed
    except Exception as e:
        if os.path.exists(csv_path):
            os.unlink(csv_path)
        if os.path.exists(r_path):
            os.unlink(r_path)
        return {'error': str(e)}


def meta_analysis(per_study_results, stratify_by=None):
    """
    Random-effects meta-analysis of per-study phase shifts.
    """
    included = per_study_results[per_study_results['included_in_meta']].copy()

    if len(included) < 2:
        return {'error': 'Too few studies for meta-analysis', 'n_studies': len(included)}

    # Compute standard errors from CIs
    included['se'] = (included['phase_shift_ci_upper'] - included['phase_shift_ci_lower']) / (2 * 1.96)

    results = {}

    if stratify_by is None:
        results['overall'] = _run_metafor(included)
    else:
        for stratum, stratum_data in included.groupby(stratify_by):
            if len(stratum_data) >= 2:
                results[str(stratum)] = _run_metafor(stratum_data)
            else:
                results[str(stratum)] = {
                    'n_studies': len(stratum_data),
                    'pooled_effect': stratum_data['phase_shift_hours'].mean() if len(stratum_data) == 1 else np.nan,
                    'error': 'Too few studies in stratum'
                }

    return results


if __name__ == "__main__":
    # Test with synthetic data
    np.random.seed(42)

    all_preds = []
    all_meta = []
    for study_idx in range(5):
        shift = np.random.uniform(-3, 3)
        for cond in ['flight', 'ground_control']:
            n = np.random.randint(3, 8)
            base_ct = np.random.uniform(0, 24)
            if cond == 'flight':
                ct = (base_ct + shift + np.random.normal(0, 2, n)) % 24
            else:
                ct = (base_ct + np.random.normal(0, 2, n)) % 24

            for i, c in enumerate(ct):
                all_preds.append({
                    'sample_name': f'sim_{study_idx}_{cond}_{i}',
                    'predicted_CT': c,
                    'circular_variance': np.random.uniform(0.1, 0.3)
                })
                all_meta.append({
                    'sample_name': f'sim_{study_idx}_{cond}_{i}',
                    'osd_id': f'OSD-{100+study_idx}',
                    'circular_variance': np.random.uniform(0.1, 0.3),
                    'condition': cond,
                    'tissue': np.random.choice(['root', 'leaf']),
                    'ecotype': 'Col-0',
                    'hardware': 'BRIC',
                    'light_regime': 'continuous_light'
                })

    preds_df = pd.DataFrame(all_preds)
    meta_df = pd.DataFrame(all_meta)

    print("Running per-study analysis...")
    per_study = per_study_analysis(preds_df, meta_df)
    print(per_study[['osd_id', 'n_flight', 'n_ground', 'phase_shift_hours', 'p_value']].to_string())

    print("\nRunning overall meta-analysis...")
    meta_results = meta_analysis(per_study)
    print(json.dumps(meta_results, indent=2, default=str))

    print("\nRunning stratified meta-analysis by tissue...")
    strat_results = meta_analysis(per_study, stratify_by='tissue')
    print(json.dumps(strat_results, indent=2, default=str))
