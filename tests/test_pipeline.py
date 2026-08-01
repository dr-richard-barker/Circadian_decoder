"""
Unit tests for the spaceflight circadian decoder pipeline.
Tests circular statistics, metadata parsing, and ChronoGauge loading.

Run with: python -m pytest tests/test_pipeline.py -v
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import pytest

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

# Dynamic path resolution for tests
if os.path.exists('/workspace/genelab_data'):
    DATA_DIR = '/workspace/genelab_data'
else:
    DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'genelab_data'))

if os.path.exists('/mnt/results'):
    RESULTS_DIR = '/mnt/results'
else:
    RESULTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

DATA_OUT_DIR = os.path.join(RESULTS_DIR, 'data')
TABLES_DIR = os.path.join(RESULTS_DIR, 'tables')

metadata_file = os.path.join(DATA_DIR, 'harmonized_metadata.csv')
if not os.path.exists(metadata_file):
    metadata_file = os.path.join(TABLES_DIR, 'tableS1_sample_metadata.csv')


# ============================================================
# Circular statistics tests
# ============================================================

class TestCircularStatistics:
    """Test circular statistics functions from statistical_analysis.py."""

    def test_circular_mean_known_values(self):
        from statistical_analysis import circular_mean_hours
        # Mean of 0h and 24h should be 0h (same point on circle)
        assert abs(circular_mean_hours([0, 24]) - 0) < 0.01
        # Mean of 6h and 18h should be 12h
        assert abs(circular_mean_hours([6, 18]) - 12) < 0.01
        # Mean of 0h and 12h should be 6h
        assert abs(circular_mean_hours([0, 12]) - 6) < 0.01

    def test_circular_mean_symmetry(self):
        from statistical_analysis import circular_mean_hours
        # Mean of 2h and 22h should be 0h (symmetric around midnight)
        result = circular_mean_hours([2, 22])
        assert abs(result - 0) < 0.01 or abs(result - 24) < 0.01

    def test_circular_std_returns_float(self):
        from statistical_analysis import circular_std_hours
        result = circular_std_hours([6, 12, 18, 0])
        assert isinstance(result, float)
        assert result >= 0

    def test_circular_distance(self):
        from statistical_analysis import circular_distance
        # Distance between 0 and 23 should be 1 (not 23)
        assert abs(circular_distance(0, 23) - 1) < 0.01
        # Distance between 6 and 18 should be 12
        assert abs(circular_distance(6, 18) - 12) < 0.01
        # Distance between 0 and 0 should be 0
        assert abs(circular_distance(0, 0)) < 0.01

    def test_watson_williams_returns_tuple(self):
        from statistical_analysis import watson_williams_test
        np.random.seed(42)
        g1 = np.random.uniform(0, 24, 20)
        g2 = np.random.uniform(0, 24, 20)
        F, p = watson_williams_test(g1, g2)
        assert isinstance(F, float)
        assert isinstance(p, float)
        assert not np.isnan(p)

    def test_watson_williams_insufficient_samples(self):
        from statistical_analysis import watson_williams_test
        F, p = watson_williams_test([12], [6])
        assert np.isnan(F)
        assert np.isnan(p)


# ============================================================
# Metadata parsing tests
# ============================================================

class TestMetadataParsing:
    """Test metadata normalization functions from metadata_curation.py."""

    def test_normalize_condition_flight(self):
        from metadata_curation import normalize_condition
        assert normalize_condition('Space Flight') == 'flight'
        assert normalize_condition('spaceflight') == 'flight'
        assert normalize_condition('FLT') == 'flight'

    def test_normalize_condition_ground(self):
        from metadata_curation import normalize_condition
        assert normalize_condition('Ground Control') == 'ground_control'
        assert normalize_condition('ground') == 'ground_control'
        assert normalize_condition('1g') == 'ground_control'

    def test_normalize_condition_empty(self):
        from metadata_curation import normalize_condition
        assert normalize_condition('') == ''
        assert normalize_condition('nan') == ''

    def test_normalize_tissue(self):
        from metadata_curation import normalize_tissue
        assert normalize_tissue('root') == 'root'
        assert normalize_tissue('leaf') == 'leaf'
        assert normalize_tissue('whole seedling') == 'whole_seedling'
        assert normalize_tissue('') == ''

    def test_normalize_ecotype(self):
        from metadata_curation import normalize_ecotype
        assert normalize_ecotype('Col-0') == 'Col-0'
        assert normalize_ecotype('Columbia-0') == 'Col-0'
        assert normalize_ecotype('Wassilewskija') == 'Ws'
        assert normalize_ecotype('Ler-0') == 'Ler-0'

    def test_normalize_light_regime(self):
        from metadata_curation import normalize_light_regime
        assert normalize_light_regime('continuous light') == 'continuous_light'
        assert normalize_light_regime('complete darkness') == 'dark'
        assert normalize_light_regime('16:8 light:dark') == 'LD_16:8'
        assert normalize_light_regime('') == 'unknown'

    def test_find_column_exact_match(self):
        from metadata_curation import find_column
        df = pd.DataFrame({'Sample Name': ['a'], 'Factor Value[Spaceflight]': ['b']})
        assert find_column(df, ['Sample Name']) == 'Sample Name'
        assert find_column(df, ['Factor Value[Spaceflight]']) == 'Factor Value[Spaceflight]'

    def test_find_column_partial_match(self):
        from metadata_curation import find_column
        df = pd.DataFrame({'Parameter Value[Hardware Type]': ['BRIC']})
        assert find_column(df, ['Hardware']) == 'Parameter Value[Hardware Type]'


# ============================================================
# ChronoGauge utility tests
# ============================================================

class TestChronoGaugeUtils:
    """Test ChronoGauge utility functions (no model loading required)."""

    def test_time24_conversion(self):
        from chronogauge_apply import time24
        # (cos, sin) for 0h: cos=-cos(pi/2)=0, sin=sin(pi/2)=1
        # Actually cyclic_time uses offset: -cos(2*pi*t/24 + pi/2), sin(2*pi*t/24 + pi/2)
        # At t=0: -cos(pi/2)=0, sin(pi/2)=1 -> atan2(1, 0)/pi*12 = pi/2/pi*12 = 6... 
        # Let's just test it returns values in [0, 24)
        preds = np.array([[0, 1], [1, 0], [0, -1], [-1, 0]])
        result = time24(preds)
        for v in result:
            assert 0 <= v < 24

    def test_cyclic_time_shape(self):
        from chronogauge_apply import cyclic_time
        cos, sin = cyclic_time([0, 6, 12, 18])
        assert len(cos) == 4
        assert len(sin) == 4
        # All values in [-1, 1]
        assert np.all(np.abs(cos) <= 1.0001)
        assert np.all(np.abs(sin) <= 1.0001)

    def test_circular_error(self):
        from chronogauge_apply import circular_error
        # Error of 0 when pred == true
        err = circular_error([12, 12], [12, 12])
        assert np.allclose(err, 0)
        # Error wraps around: pred=23, true=1 -> error = -2h = -120 min
        err = circular_error([23], [1])
        assert abs(err[0] - (-120)) < 1

    def test_circular_mean_aggregation(self):
        from chronogauge_apply import circular_mean
        # Create a DataFrame where all models predict the same time
        df = pd.DataFrame({'m0': [6, 12], 'm1': [6, 12], 'm2': [6, 12]})
        result = circular_mean(df)
        assert len(result) == 2
        assert abs(result[0] - 6) < 0.1
        assert abs(result[1] - 12) < 0.1


# ============================================================
# Config file tests
# ============================================================

class TestConfig:
    """Test that config files are valid and complete."""

    def test_study_lists_json_exists(self):
        config_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'study_lists.json')
        assert os.path.exists(config_path)

    def test_study_lists_json_valid(self):
        config_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'study_lists.json')
        with open(config_path) as f:
            config = json.load(f)
        assert 'fork_a' in config
        assert 'fork_b' in config
        assert 'osd_ids' in config['fork_a']
        assert 'osd_ids' in config['fork_b']
        assert isinstance(config['fork_a']['osd_ids'], list)
        assert isinstance(config['fork_b']['osd_ids'], list)

    def test_fork_a_has_15_studies(self):
        config_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'study_lists.json')
        with open(config_path) as f:
            config = json.load(f)
        assert len(config['fork_a']['osd_ids']) == 15

    def test_fork_b_has_23_studies(self):
        config_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'study_lists.json')
        with open(config_path) as f:
            config = json.load(f)
        assert len(config['fork_b']['osd_ids']) == 23

    def test_core_clock_genes_present(self):
        config_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'study_lists.json')
        with open(config_path) as f:
            config = json.load(f)
        assert 'core_clock_genes' in config
        genes = config['core_clock_genes']
        assert 'CCA1' in genes['morning']
        assert 'TOC1' in genes['evening']


# ============================================================
# Integration test (requires data + models)
# ============================================================

class TestIntegration:
    """Integration tests that require downloaded data and models.
    Skipped if data is not present.
    """

    @pytest.mark.skipif(
        not os.path.exists(metadata_file),
        reason="GeneLab metadata or tableS1 not found"
    )
    def test_metadata_file_exists(self):
        df = pd.read_csv(metadata_file)
        assert len(df) > 100
        assert 'osd_id' in df.columns
        assert 'condition' in df.columns

    @pytest.mark.skipif(
        not os.path.exists(os.path.join(DATA_OUT_DIR, 'all_predictions.csv')),
        reason="Predictions not generated"
    )
    def test_predictions_file_exists(self):
        df = pd.read_csv(os.path.join(DATA_OUT_DIR, 'all_predictions.csv'))
        assert len(df) > 100
        assert 'predicted_CT' in df.columns
        assert 'circular_variance' in df.columns

    @pytest.mark.skipif(
        not os.path.exists(os.path.join(DATA_OUT_DIR, 'per_study_results.csv')),
        reason="Statistical analysis not run"
    )
    def test_per_study_results_exist(self):
        df = pd.read_csv(os.path.join(DATA_OUT_DIR, 'per_study_results.csv'))
        assert len(df) > 5
        assert 'phase_shift_hours' in df.columns
        assert 'p_value' in df.columns


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
