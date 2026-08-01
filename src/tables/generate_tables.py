"""
Table generation module for the spaceflight circadian analysis.
Generates all main and supplementary tables as CSV and formatted text.
"""
import os
import numpy as np
import pandas as pd
from pathlib import Path


def table1_dataset_characteristics(metadata_df, output_dir):
    """Table 1: Dataset characteristics for all studies."""
    rows = []

    for osd_id, study in metadata_df.groupby('osd_id'):
        platform = study['assay_technology'].iloc[0] if 'assay_technology' in study.columns else ''
        n_flight = (study['condition'] == 'flight').sum()
        n_ground = (study['condition'] == 'ground_control').sum()
        tissue = study['tissue'].mode().iloc[0] if len(study['tissue'].dropna()) > 0 else ''
        ecotype = study['ecotype'].mode().iloc[0] if len(study['ecotype'].dropna()) > 0 else ''
        hardware = study['hardware'].mode().iloc[0] if len(study['hardware'].dropna()) > 0 else ''
        light_regime = study['light_regime'].mode().iloc[0] if len(study['light_regime'].dropna()) > 0 else ''

        rows.append({
            'OSD_ID': osd_id,
            'Platform': platform,
            'Tissue': tissue,
            'Ecotype': ecotype,
            'Hardware': hardware,
            'Light_Regime': light_regime,
            'N_Flight': n_flight,
            'N_Ground': n_ground,
            'N_Total': len(study),
        })

    table = pd.DataFrame(rows).sort_values('OSD_ID')
    os.makedirs(output_dir, exist_ok=True)
    table.to_csv(os.path.join(output_dir, 'table1_dataset_characteristics.csv'), index=False)
    print(f"  Saved table1_dataset_characteristics.csv ({len(table)} studies)")
    return table


def table2_per_study_results(per_study_df, output_dir):
    """Table 2: Per-study circadian phase shift results."""
    table = per_study_df.copy()
    # Format columns
    table['phase_shift_hours'] = table['phase_shift_hours'].round(2)
    table['phase_shift_ci_lower'] = table['phase_shift_ci_lower'].round(2)
    table['phase_shift_ci_upper'] = table['phase_shift_ci_upper'].round(2)
    table['p_value'] = table['p_value'].apply(lambda x: f'{x:.4f}' if pd.notna(x) else 'NA')
    table['flight_mean_CT'] = table['flight_mean_CT'].round(2)
    table['ground_mean_CT'] = table['ground_mean_CT'].round(2)

    # Rename for publication
    table = table.rename(columns={
        'osd_id': 'OSD_ID',
        'n_flight': 'N_Flight',
        'n_ground': 'N_Ground',
        'flight_mean_CT': 'Flight_Mean_CT',
        'ground_mean_CT': 'Ground_Mean_CT',
        'phase_shift_hours': 'Phase_Shift_h',
        'phase_shift_ci_lower': 'CI_95_Lower',
        'phase_shift_ci_upper': 'CI_95_Upper',
        'p_value': 'P_Value',
        'tissue': 'Tissue',
        'ecotype': 'Ecotype',
        'hardware': 'Hardware',
        'light_regime': 'Light_Regime',
        'included_in_meta': 'Included_in_Meta',
        'reason': 'Notes',
    })

    os.makedirs(output_dir, exist_ok=True)
    table.to_csv(os.path.join(output_dir, 'table2_per_study_results.csv'), index=False)
    print(f"  Saved table2_per_study_results.csv ({len(table)} studies)")
    return table


def table3_meta_analysis(meta_results, output_dir, stratify_label='overall'):
    """Table 3: Meta-analysis results by stratum."""
    rows = []

    for stratum, result in meta_results.items():
        if 'error' in result:
            rows.append({
                'Stratum': stratum,
                'N_Studies': result.get('n_studies', np.nan),
                'Pooled_Effect_h': np.nan,
                'CI_95_Lower': np.nan,
                'CI_95_Upper': np.nan,
                'P_Value': np.nan,
                'I2_pct': np.nan,
                'Tau2': np.nan,
                'Notes': result['error']
            })
        else:
            rows.append({
                'Stratum': stratum,
                'N_Studies': result.get('N_STUDIES', np.nan),
                'Pooled_Effect_h': result.get('POOLED_EFFECT', np.nan),
                'CI_95_Lower': result.get('CI_LOWER', np.nan),
                'CI_95_Upper': result.get('CI_UPPER', np.nan),
                'P_Value': result.get('P_VALUE', np.nan),
                'I2_pct': result.get('I2', np.nan),
                'Tau2': result.get('TAU2', np.nan),
                'Notes': ''
            })

    table = pd.DataFrame(rows)

    # Round numeric columns
    for col in ['Pooled_Effect_h', 'CI_95_Lower', 'CI_95_Upper', 'I2_pct', 'Tau2']:
        table[col] = table[col].round(3)
    table['P_Value'] = table['P_Value'].apply(lambda x: f'{x:.4f}' if pd.notna(x) else 'NA')

    os.makedirs(output_dir, exist_ok=True)
    fname = f'table3_meta_analysis_{stratify_label}.csv'
    table.to_csv(os.path.join(output_dir, fname), index=False)
    print(f"  Saved {fname}")
    return table


def tableS1_sample_metadata(metadata_df, output_dir):
    """Table S1: Full sample-level metadata."""
    table = metadata_df.copy()
    os.makedirs(output_dir, exist_ok=True)
    table.to_csv(os.path.join(output_dir, 'tableS1_sample_metadata.csv'), index=False)
    print(f"  Saved tableS1_sample_metadata.csv ({len(table)} samples)")
    return table


def tableS2_chronogauge_validation(validation_df, output_dir):
    """Table S2: ChronoGauge validation results."""
    # Per-study validation metrics
    rows = []
    for osd_id, study in validation_df.groupby('osd_id'):
        errors = study['error_minutes'].values
        rows.append({
            'OSD_ID': osd_id,
            'N_Samples': len(study),
            'MAE_minutes': np.mean(np.abs(errors)),
            'Median_Error_minutes': np.median(np.abs(errors)),
            'RMSE_minutes': np.sqrt(np.mean(errors**2)),
            'R_value': study['true_CT'].corr(study['predicted_CT']) if len(study) > 2 else np.nan,
        })

    table = pd.DataFrame(rows)
    for col in ['MAE_minutes', 'Median_Error_minutes', 'RMSE_minutes', 'R_value']:
        table[col] = table[col].round(2)

    os.makedirs(output_dir, exist_ok=True)
    table.to_csv(os.path.join(output_dir, 'tableS2_chronogauge_validation.csv'), index=False)
    print(f"  Saved tableS2_chronogauge_validation.csv")
    return table


if __name__ == "__main__":
    print("Table generation module loaded successfully.")
    print("Available functions:")
    for name in dir():
        if name.startswith('table'):
            print(f"  - {name}")
