"""
Metadata curation module.
Parses ISA-Tab metadata from GeneLab studies and harmonizes into a
single sample-level table with tissue, genotype, hardware, light regime, ZT, etc.

Handles the actual ISA-Tab column structure observed in GeneLab studies:
- Factor Value[Spaceflight]: 'Space Flight' / 'Ground Control'
- Factor Value[Ecotype] or Characteristics[ecotype] or Characteristics[Strain]
- Parameter Value[Hardware]
- Parameter Value[light cycle] or Parameter Value[light regimen]
- Characteristics[Material Type] or Characteristics[organism part] or Factor Value[Tissue]
- Parameter Value[Age at sample harvest] or Parameter Value[Age at sample collection]
"""
import os
import re
import json
import zipfile
import numpy as np
import pandas as pd
from pathlib import Path


def find_isa_files(study_dir):
    """Find ISA-Tab s_*.txt sample files in a study directory."""
    isa_files = []
    # Check metadata subdirectory (from extracted zip)
    meta_dir = os.path.join(study_dir, 'metadata')
    if os.path.isdir(meta_dir):
        for root, dirs, files in os.walk(meta_dir):
            for f in files:
                if f.startswith('s_') and f.endswith('.txt'):
                    isa_files.append(os.path.join(root, f))
    # Check direct files
    if not isa_files:
        for f in os.listdir(study_dir):
            if f.startswith('s_') and f.endswith('.txt'):
                isa_files.append(os.path.join(study_dir, f))
    return isa_files


def parse_isa_tab(isa_path):
    """Parse an ISA-Tab sample file into a DataFrame."""
    try:
        df = pd.read_csv(isa_path, sep='\t', encoding='utf-8', dtype=str)
        return df
    except Exception as e:
        print(f"  Error parsing {isa_path}: {e}")
        return None


def find_column(df, patterns):
    """Find the first matching column name from a list of patterns (case-insensitive)."""
    for p in patterns:
        for c in df.columns:
            if p.lower() == c.lower():
                return c
    # Partial match
    for p in patterns:
        for c in df.columns:
            if p.lower() in c.lower():
                return c
    return None


def normalize_condition(value):
    """Normalize spaceflight condition values."""
    v = str(value).lower().strip()
    if 'space flight' in v or 'spaceflight' in v or 'flight' in v or 'flt' in v:
        return 'flight'
    elif 'ground control' in v or 'ground' in v or 'gc' in v or '1g' in v:
        return 'ground_control'
    elif 'space' in v:
        return 'flight'
    elif v in ('', 'nan', 'none'):
        return ''
    return v


def normalize_tissue(value, developmental_stage=''):
    """Normalize tissue/material type values."""
    v = str(value).lower().strip()
    if v in ('', 'nan', 'none'):
        v = str(developmental_stage).lower().strip()

    if 'root' in v:
        if 'zone' in v:
            return 'root_zone'
        return 'root'
    elif 'leaf' in v or 'leaves' in v:
        return 'leaf'
    elif 'shoot' in v:
        return 'shoot'
    elif 'hypocotyl' in v:
        return 'hypocotyl'
    elif 'seedling' in v or 'whole' in v or 'plant' in v:
        return 'whole_seedling'
    elif 'callus' in v or 'culture' in v:
        return 'cell_culture'
    elif 'cotyledon' in v:
        return 'cotyledon'
    elif 'stem' in v:
        return 'stem'
    elif 'flower' in v or 'inflorescence' in v:
        return 'flower'
    elif v and v not in ('', 'nan', 'none'):
        return v
    return ''


def normalize_ecotype(value, strain='', genotype=''):
    """Normalize ecotype/accession values."""
    v = str(value).strip()
    if v in ('', 'nan', 'None'):
        v = str(strain).strip()
    if v in ('', 'nan', 'None'):
        v = str(genotype).strip()

    # Map full names to standard abbreviations
    v_lower = v.lower()
    if 'wassilewskija' in v_lower or 'wasselewskija' in v_lower or v_lower == 'ws-0' or v_lower == 'ws':
        return 'Ws'
    elif 'col-0' in v_lower or 'col0' in v_lower or 'columbia' in v_lower:
        if 'phyd' in v_lower:
            return 'Col-0 phyD'
        return 'Col-0'
    elif 'ler' in v_lower:
        return 'Ler-0'
    elif 'cvi' in v_lower:
        return 'Cvi-0'
    elif 'ws-2' in v_lower:
        return 'Ws-2'
    elif 'wild type' in v_lower or 'wildtype' in v_lower:
        return 'wild_type'
    elif v and v not in ('', 'nan', 'None'):
        return v
    return ''


def normalize_light_regime(light_cycle, light_regimen=''):
    """Normalize light regime values."""
    v = str(light_cycle).lower().strip()
    if v in ('', 'nan', 'none', 'not available'):
        v = str(light_regimen).lower().strip()

    if 'continuous' in v or '24 light' in v or '24light' in v or v == '24' or 'll' in v or '24 hr light' in v or '24 white light' in v or '24white' in v:
        return 'continuous_light'
    elif 'complete darkness' in v or ('dark' in v and 'light' not in v) or v == 'dark':
        return 'dark'
    elif '16' in v and '8' in v:
        return 'LD_16:8'
    elif '12' in v and '12' in v:
        return 'LD_12:12'
    elif '4:2' in v or '4:2 light:dark' in v:
        return 'LD_4:2'
    elif 'photoperiod' in v or 'light dark' in v or 'light:dark' in v or 'ld' in v:
        return 'photoperiod'
    elif v and v not in ('', 'nan', 'none', 'not available'):
        return v
    return 'unknown'


def extract_harvest_age(df):
    """Extract harvest age from various possible columns."""
    age_col = find_column(df, [
        'Parameter Value[Age at sample harvest]',
        'Parameter Value[Age at sample collection]',
        'Parameter Value[Age at Sample Harvest]',
        'Characteristics[age]',
        'Factor Value[Age]',
    ])
    if age_col:
        return df[age_col].apply(lambda x: _parse_age(x))
    return pd.Series([np.nan] * len(df))


def _parse_age(value):
    """Parse age value to days (numeric)."""
    v = str(value).strip()
    if v in ('', 'nan', 'None'):
        return np.nan
    # Extract number
    m = re.search(r'(\d+)', v)
    if m:
        return float(m.group(1))
    return np.nan


def extract_sample_metadata(df, osd_id):
    """Extract harmonized metadata from an ISA-Tab DataFrame."""
    samples = []

    # Find columns using flexible matching
    sample_name_col = find_column(df, ['Sample Name'])
    organism_col = find_column(df, ['Characteristics[organism]', 'Characteristics[Organism]'])

    # Condition (Spaceflight)
    condition_col = find_column(df, [
        'Factor Value[Spaceflight]',
        'Factor Value[Space Flight]',
        'Factor Value[Condition]',
        'Factor Value[Treatment]',
        'Factor Value[Gravity]',
    ])

    # Ecotype/Genotype
    ecotype_col = find_column(df, [
        'Factor Value[Ecotype]',
        'Characteristics[ecotype]',
        'Characteristics[Ecotype]',
        'Characteristics[Cultivar]',
        'Factor Value[Cultivar]',
        'Characteristics[Strain]',
        'Characteristics[Accession]',
    ])
    genotype_col = find_column(df, [
        'Characteristics[genotype]',
        'Characteristics[Genotype]',
        'Factor Value[Genotype]',
    ])

    # Tissue
    tissue_col = find_column(df, [
        'Characteristics[Material Type]',
        'Factor Value[Material Type]',
        'Characteristics[organism part]',
        'Factor Value[Tissue]',
        'Characteristics[Tissue]',
        'Factor Value[Organism Part]',
    ])
    dev_stage_col = find_column(df, [
        'Characteristics[developmental stage]',
        'Characteristics[Developmental Stage]',
    ])
    factor_tissue_col = find_column(df, ['Factor Value[Tissue]'])
    if not tissue_col and factor_tissue_col:
        tissue_col = factor_tissue_col

    # Hardware
    hardware_col = find_column(df, [
        'Parameter Value[Hardware]',
        'Parameter Value[hardware]',
    ])

    # Light regime
    light_cycle_col = find_column(df, [
        'Parameter Value[light cycle]',
        'Parameter Value[Light Cycle]',
        'Parameter Value[Growth Light Cycle]',
        'Parameter Value[growth light cycle]',
        'Parameter Value[light]',
        'Parameter Value[Light]',
    ])
    light_regimen_col = find_column(df, [
        'Parameter Value[light regimen]',
        'Parameter Value[Light Regimen]',
        'Parameter Value[Growth Light Regimen]',
    ])

    # Growth environment (helps identify ground-only studies)
    growth_env_col = find_column(df, [
        'Parameter Value[growth environment]',
        'Parameter Value[Growth Environment]',
    ])

    for _, row in df.iterrows():
        sample = {'osd_id': f'OSD-{osd_id}'}

        # Sample name
        sample['sample_name'] = str(row[sample_name_col]).strip() if sample_name_col else ''

        # Organism
        sample['organism'] = str(row[organism_col]).strip() if organism_col else ''

        # Condition
        cond_val = str(row[condition_col]).strip() if condition_col else ''
        growth_env = str(row[growth_env_col]).strip() if growth_env_col else ''
        sample['condition'] = normalize_condition(cond_val)
        # If no Spaceflight factor but growth environment says Earth surface, it's ground-only
        if not sample['condition'] and 'earth' in growth_env.lower():
            sample['condition'] = 'ground_control'

        # Ecotype
        eco_val = str(row[ecotype_col]).strip() if ecotype_col else ''
        geno_val = str(row[genotype_col]).strip() if genotype_col else ''
        sample['ecotype'] = normalize_ecotype(eco_val, strain='', genotype=geno_val)

        # Also keep raw genotype
        sample['genotype'] = geno_val if geno_val not in ('', 'nan', 'None') else ''

        # Tissue
        tissue_val = str(row[tissue_col]).strip() if tissue_col else ''
        dev_val = str(row[dev_stage_col]).strip() if dev_stage_col else ''
        sample['tissue'] = normalize_tissue(tissue_val, dev_val)

        # Hardware
        sample['hardware'] = str(row[hardware_col]).strip() if hardware_col else ''
        if osd_id == 522:
            sample['hardware'] = 'BRIC-LED'
        elif sample['hardware'] in ('', 'nan', 'None'):
            sample['hardware'] = ''

        # Light regime
        lc_val = str(row[light_cycle_col]).strip() if light_cycle_col else ''
        lr_val = str(row[light_regimen_col]).strip() if light_regimen_col else ''
        sample['light_regime'] = normalize_light_regime(lc_val, lr_val)

        # Harvest age
        sample['harvest_age_days'] = _parse_age(row.get(age_col, '')) if (age_col := find_column(df, [
            'Parameter Value[Age at sample harvest]',
            'Parameter Value[Age at sample collection]',
            'Parameter Value[Age at Sample Harvest]',
            'Characteristics[age]',
            'Factor Value[Age]',
        ])) else np.nan

        # ZT: most GeneLab studies don't have explicit ZT
        sample['zt_known'] = False
        sample['zt_hours'] = np.nan

        samples.append(sample)

    return samples


def curate_all_metadata(data_dir, study_info_path=None):
    """Curate metadata for all downloaded studies."""
    all_samples = []

    # Load study info if available
    study_info = {}
    if study_info_path and os.path.exists(study_info_path):
        try:
            with open(study_info_path) as f:
                study_info_list = json.load(f)
                study_info = {s['osd_id']: s for s in study_info_list}
        except Exception as e:
            print(f"  Error loading study info JSON: {e}")

    if 'OSD-522' not in study_info:
        study_info['OSD-522'] = {
            'title': 'Integrative Transcriptomics and Proteomics Profiling of Arabidopsis thaliana Elucidates Novel Mechanisms Underlying Spaceflight Adaptation Study',
            'assay_technology': 'RNA Sequencing (RNA-Seq)',
            'platform': 'Illumina NovaSeq 6000',
            'study_type': 'Transcription Profiling'
        }

    # Find all study directories
    study_dirs = [d for d in os.listdir(data_dir) if d.startswith('OSD-') and os.path.isdir(os.path.join(data_dir, d))]

    for study_dir_name in sorted(study_dirs, key=lambda x: int(x.split('-')[1])):
        osd_id = int(study_dir_name.split('-')[1])
        study_dir = os.path.join(data_dir, study_dir_name)

        isa_files = find_isa_files(study_dir)
        if not isa_files:
            print(f"  No ISA files for {study_dir_name}")
            continue

        for isa_path in isa_files:
            df = parse_isa_tab(isa_path)
            if df is not None and len(df) > 0:
                samples = extract_sample_metadata(df, osd_id)
                all_samples.extend(samples)
                print(f"  {study_dir_name}: {len(samples)} samples from {os.path.basename(isa_path)}")

    metadata_df = pd.DataFrame(all_samples)

    # Add study-level info
    if study_info:
        metadata_df['study_title'] = metadata_df['osd_id'].map(
            lambda x: study_info.get(x, {}).get('title', '')
        )
        metadata_df['assay_technology'] = metadata_df['osd_id'].map(
            lambda x: study_info.get(x, {}).get('assay_technology', '')
        )
        metadata_df['platform'] = metadata_df['osd_id'].map(
            lambda x: study_info.get(x, {}).get('platform', '')
        )
        metadata_df['study_type'] = metadata_df['osd_id'].map(
            lambda x: study_info.get(x, {}).get('study_type', '')
        )

    return metadata_df


if __name__ == "__main__":
    if os.path.exists("/workspace/genelab_data"):
        DATA_DIR = "/workspace/genelab_data"
    else:
        DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "genelab_data"))
    STUDY_INFO = os.path.join(DATA_DIR, "arabidopsis_transcriptomics.json")

    print("Curating metadata for all studies...")
    metadata = curate_all_metadata(DATA_DIR, STUDY_INFO)

    print(f"\nTotal samples: {len(metadata)}")
    print(f"Studies: {metadata['osd_id'].nunique()}")
    print(f"\nCondition distribution:\n{metadata['condition'].value_counts()}")
    print(f"\nTissue distribution:\n{metadata['tissue'].value_counts()}")
    print(f"\nEcotype distribution:\n{metadata['ecotype'].value_counts().head(10)}")
    print(f"\nHardware distribution:\n{metadata['hardware'].value_counts()}")
    print(f"\nLight regime distribution:\n{metadata['light_regime'].value_counts()}")

    output_path = os.path.join(DATA_DIR, "harmonized_metadata.csv")
    metadata.to_csv(output_path, index=False)
    print(f"\nSaved to {output_path}")
