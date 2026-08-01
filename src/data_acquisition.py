"""
Data acquisition from NASA GeneLab OSDR - corrected version.
Downloads processed RNA-seq/microarray normalized data and ISA-Tab metadata
for Arabidopsis spaceflight studies.
"""
import requests
import json
import os
import time
import zipfile
import io
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

FILE_API = "https://osdr.nasa.gov/osdr/data/osd/files"
BDATA_API = "https://visualization.osdr.nasa.gov/biodata/api"

# Fork A: 15 accessions from Barker et al. 2023
FORK_A_GLDS = [7, 17, 37, 38, 44, 46, 120, 121, 136, 147, 205, 208, 213, 218, 251]
# GLDS-213 not found in OSD; map the rest directly
FORK_A_OSD = [7, 17, 37, 38, 44, 46, 120, 121, 136, 147, 205, 208, 218, 251]  # 14 (213 missing)


def get_study_files(osd_id):
    """Get file listing for a study via the file API."""
    url = f"{FILE_API}/{osd_id}?all_files=true"
    try:
        r = requests.get(url, timeout=120)
        if r.status_code == 200:
            data = r.json()
            study_key = f"OSD-{osd_id}"
            if study_key in data.get("studies", {}):
                return data["studies"][study_key]
    except Exception as e:
        print(f"  Error fetching files for OSD-{osd_id}: {e}")
    return None


def download_file(url, dest_path, max_retries=3):
    """Download a file with retries. Handles relative URLs."""
    if url.startswith('/'):
        url = 'https://osdr.nasa.gov' + url
    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=600, stream=True)
            if r.status_code == 200:
                with open(dest_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=65536):
                        f.write(chunk)
                return True
            else:
                print(f"    HTTP {r.status_code} for {url[-80:]}")
        except Exception as e:
            print(f"    Attempt {attempt+1} error: {e}")
            time.sleep(5)
    return False


def categorize_files(study_files):
    """Categorize files by type."""
    result = {
        'rnaseq_normalized': [],
        'rnaseq_deg': [],
        'microarray_normalized': [],
        'microarray_deg': [],
        'metadata_zip': [],
        'isa_txt': [],
        'other_processed': []
    }

    for f in study_files:
        fname = f.get('file_name', '')
        furl = f.get('remote_url', '')
        fcat = f.get('category', '')

        if 'Normalized_Counts' in fname and fname.endswith('.csv'):
            result['rnaseq_normalized'].append((fname, furl))
        elif 'differential_expression' in fname and fname.endswith('.csv'):
            result['rnaseq_deg'].append((fname, furl))
        elif 'normalized_expression' in fname and fname.endswith('.csv'):
            result['microarray_normalized'].append((fname, furl))
        elif 'normalized_intensities_probe' in fname and fname.endswith('.csv'):
            pass  # Skip probe-level intensities (too large, we use probeset-level)
        elif '_metadata_' in fname and fname.endswith('.zip'):
            result['metadata_zip'].append((fname, furl))
        elif fname.startswith('s_') and fname.endswith('.txt'):
            result['isa_txt'].append((fname, furl))
        elif 'GeneLab Processed' in fcat and fname.endswith('.csv') and 'raw' not in fname.lower() and 'differential_expression' not in fname.lower():
            result['other_processed'].append((fname, furl))

    return result


def download_study(osd_id, base_dir):
    """Download all relevant files for a study."""
    study_dir = Path(base_dir) / f"OSD-{osd_id}"
    study_dir.mkdir(parents=True, exist_ok=True)

    study_data = get_study_files(osd_id)
    if not study_data:
        print(f"  No file data for OSD-{osd_id}")
        return None

    files = study_data.get('study_files', [])
    categorized = categorize_files(files)

    result = {
        'osd_id': osd_id,
        'file_count': study_data.get('file_count', 0),
        'categories': {k: len(v) for k, v in categorized.items()},
        'files': {}
    }

    for cat_name, file_list in categorized.items():
        if not file_list:
            continue
        result['files'][cat_name] = []
        for fname, furl in file_list:
            dest = study_dir / fname
            if not dest.exists() or dest.stat().st_size == 0:
                print(f"  Downloading {fname}...")
                if download_file(furl, dest):
                    result['files'][cat_name].append(str(dest))
            else:
                result['files'][cat_name].append(str(dest))

    return result


def extract_metadata_zip(zip_path, study_dir):
    """Extract ISA-Tab metadata from zip file."""
    try:
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(study_dir / 'metadata')
        return True
    except Exception as e:
        print(f"  Error extracting {zip_path}: {e}")
        return False


if __name__ == "__main__":
    if os.path.exists("/workspace/genelab_data"):
        DATA_DIR = "/workspace/genelab_data"
    else:
        DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "genelab_data"))
    Path(DATA_DIR).mkdir(parents=True, exist_ok=True)

    # Load Fork B IDs (all Arabidopsis spaceflight transcriptomics)
    with open(f"{DATA_DIR}/fork_b_osd_ids.json") as f:
        fork_b_ids = json.load(f)

    # Fork A is subset of Fork B plus some non-spaceflight studies
    # Add Fork A IDs that aren't in Fork B
    fork_a_set = set(FORK_A_OSD)
    fork_b_set = set(int(x.split('-')[1]) for x in fork_b_ids)
    all_ids = sorted(fork_a_set | fork_b_set)

    print(f"Fork A: {len(fork_a_set)} studies")
    print(f"Fork B: {len(fork_b_set)} studies")
    print(f"Combined unique: {len(all_ids)} studies")
    print()

    # Download all studies
    all_results = {}
    for osd_id in all_ids:
        fork = "A+B" if osd_id in fork_a_set and osd_id in fork_b_set else ("A" if osd_id in fork_a_set else "B")
        print(f"\nProcessing OSD-{osd_id} (Fork {fork})...")
        result = download_study(osd_id, DATA_DIR)
        if result:
            all_results[osd_id] = result
            cats = result['categories']
            print(f"  RNA-seq norm: {cats['rnaseq_normalized']}, microarray norm: {cats['microarray_normalized']}, "
                  f"metadata: {cats['metadata_zip']}, DEG: {cats['rnaseq_deg']}")

            # Extract metadata zips
            for mp in result['files'].get('metadata_zip', []):
                extract_metadata_zip(mp, Path(DATA_DIR) / f"OSD-{osd_id}")

    # Save download log
    with open(f"{DATA_DIR}/download_log.json", "w") as f:
        json.dump(all_results, f, indent=2)

    # Summary
    print(f"\n{'='*60}")
    print(f"Download complete: {len(all_results)} studies")
    rnaseq_studies = [k for k, v in all_results.items() if v['categories']['rnaseq_normalized'] > 0]
    microarray_studies = [k for k, v in all_results.items() if v['categories']['microarray_normalized'] > 0]
    print(f"  RNA-seq studies: {len(rnaseq_studies)} -> {rnaseq_studies}")
    print(f"  Microarray studies: {len(microarray_studies)} -> {microarray_studies}")
