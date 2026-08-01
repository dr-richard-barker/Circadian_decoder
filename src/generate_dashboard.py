"""
Generate a self-contained interactive Plotly dashboard for GitHub Pages.
Produces static HTML + data.js + local plotly.min.js that can be deployed
to GitHub Pages (docs/ folder convention) with no server-side processing.

Dashboard views:
  1. Circadian phase overview (predicted CT vs circular variance, colored by condition)
  2. Phase shift forest plot (per-study with CIs and pooled estimate)
  3. Circadian fingerprint (circular variance by condition + study)
  4. Clock gene expression heatmap (11 core clock genes)
  5. t-SNE trajectory (from trajectory analysis)
"""
import os
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import requests

src_dir = os.path.dirname(os.path.abspath(__file__))

if os.path.exists("/mnt/results"):
    RESULTS_DIR = "/mnt/results"
elif os.path.exists("/results"):
    RESULTS_DIR = "/results"
else:
    RESULTS_DIR = os.path.abspath(os.path.join(src_dir, ".."))

DATA_OUT_DIR = os.path.join(RESULTS_DIR, "data")

if os.path.exists("/workspace/spaceflight-circadian-decoder/docs"):
    DOCS_DIR = "/workspace/spaceflight-circadian-decoder/docs"
else:
    DOCS_DIR = os.path.abspath(os.path.join(src_dir, "..", "dashboard"))

if os.path.exists("/workspace/genelab_data"):
    DATA_DIR = "/workspace/genelab_data"
else:
    DATA_DIR = os.path.abspath(os.path.join(src_dir, "..", "genelab_data"))

os.makedirs(DOCS_DIR, exist_ok=True)

# Phylo color palette
PHYLO_COLORS = {
    'blue': '#0279EE', 'orange': '#FF9400', 'green': '#75A025',
    'pink': '#FD9BED', 'yellow': '#E9ED4C', 'black': '#000000',
    'cream': '#ECE9E2', 'white': '#FAF9F3',
}


def download_plotly_js():
    """Download plotly.min.js for local hosting (no CDN dependency)."""
    js_path = os.path.join(DOCS_DIR, 'plotly.min.js')
    if os.path.exists(js_path) and os.path.getsize(js_path) > 1000000:
        print("  plotly.min.js already exists")
        return
    print("  Downloading plotly.min.js...")
    try:
        r = requests.get('https://cdn.plot.ly/plotly-2.35.2.min.js', timeout=60)
        if r.status_code == 200:
            with open(js_path, 'w') as f:
                f.write(r.text)
            print(f"  Downloaded plotly.min.js ({len(r.text):,} chars)")
        else:
            print(f"  HTTP {r.status_code}, will use CDN fallback")
    except Exception as e:
        print(f"  Download failed: {e}, will use CDN fallback")


def prepare_data():
    """Load and prepare all data for the dashboard."""
    print("Loading data...")

    # Predictions
    predictions = pd.read_csv(os.path.join(DATA_OUT_DIR, 'all_predictions.csv'))
    metadata_path = os.path.join(DATA_DIR, 'harmonized_metadata.csv')
    if not os.path.exists(metadata_path):
        metadata_path = os.path.join(RESULTS_DIR, 'tables', 'tableS1_sample_metadata.csv')
    metadata = pd.read_csv(metadata_path)

    # Merge predictions with metadata
    merged = predictions.merge(
        metadata[['sample_name', 'osd_id', 'condition', 'tissue', 'light_regime',
                  'ecotype', 'hardware']],
        on=['sample_name', 'osd_id'], how='left'
    )

    # Per-study results
    per_study = pd.read_csv(os.path.join(DATA_OUT_DIR, 'per_study_results.csv'))

    # Meta-analysis results
    with open(os.path.join(DATA_OUT_DIR, 'meta_analysis_results.json')) as f:
        meta_results = json.load(f)

    # Clock gene expression
    clock_expr = pd.read_csv(os.path.join(DATA_OUT_DIR, 'clock_gene_expression.csv'))

    # Trajectory analysis
    trajectory = None
    traj_path = os.path.join(DATA_OUT_DIR, 'trajectory_analysis.csv')
    if os.path.exists(traj_path):
        trajectory = pd.read_csv(traj_path)

    # t-SNE results (recompute from model predictions for the dashboard)
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler
    model_preds = pd.read_csv(os.path.join(DATA_OUT_DIR, 'all_model_predictions.csv'))
    feature_cols = [c for c in model_preds.columns if c not in ['sample_name', 'osd_id']]
    X = model_preds[feature_cols].fillna(12).values
    X_scaled = StandardScaler().fit_transform(X)
    tsne = TSNE(n_components=2, perplexity=30, max_iter=1000, random_state=42)
    tsne_results = tsne.fit_transform(X_scaled)

    tsne_df = model_preds[['sample_name', 'osd_id']].copy()
    tsne_df['tSNE1'] = tsne_results[:, 0]
    tsne_df['tSNE2'] = tsne_results[:, 1]
    tsne_df = tsne_df.merge(
        metadata[['sample_name', 'osd_id', 'condition', 'tissue']],
        on=['sample_name', 'osd_id'], how='left'
    )
    tsne_df = tsne_df.merge(predictions[['sample_name', 'osd_id', 'predicted_CT']],
                            on=['sample_name', 'osd_id'], how='left')

    return merged, per_study, meta_results, clock_expr, trajectory, tsne_df


def create_dashboard_data_js(merged, per_study, meta_results, clock_expr, trajectory, tsne_df):
    """Create data.js with all data embedded as JSON."""
    print("Creating data.js...")

    # Prepare scatter data for phase overview
    scatter_data = []
    for _, row in merged.iterrows():
        scatter_data.append({
            'sample_name': str(row['sample_name']),
            'osd_id': str(row['osd_id']),
            'condition': str(row.get('condition', '')),
            'tissue': str(row.get('tissue', '')),
            'light_regime': str(row.get('light_regime', '')),
            'ecotype': str(row.get('ecotype', '')),
            'predicted_CT': float(row['predicted_CT']) if pd.notna(row['predicted_CT']) else None,
            'circular_variance': float(row['circular_variance']) if pd.notna(row['circular_variance']) else None,
        })

    # Prepare forest plot data
    forest_data = []
    for _, row in per_study.iterrows():
        forest_data.append({
            'osd_id': str(row['osd_id']),
            'n_flight': int(row['n_flight']) if pd.notna(row['n_flight']) else 0,
            'n_ground': int(row['n_ground']) if pd.notna(row['n_ground']) else 0,
            'phase_shift': float(row['phase_shift_hours']) if pd.notna(row['phase_shift_hours']) else None,
            'ci_lower': float(row['phase_shift_ci_lower']) if pd.notna(row['phase_shift_ci_lower']) else None,
            'ci_upper': float(row['phase_shift_ci_upper']) if pd.notna(row['phase_shift_ci_upper']) else None,
            'p_value': float(row['p_value']) if pd.notna(row['p_value']) else None,
            'tissue': str(row.get('tissue', '')),
            'light_regime': str(row.get('light_regime', '')),
            'included': bool(row.get('included_in_meta', False)),
        })

    # Pooled effect
    pooled = meta_results.get('overall', {}).get('overall', {})

    # Prepare clock gene heatmap data (subset for performance)
    # Get unique genes and a subset of samples
    clock_genes = sorted(clock_expr['gene_name'].dropna().unique())
    # Sample subset: take up to 200 samples, balanced flight/ground
    clock_samples = clock_expr[['sample_name', 'osd_id']].drop_duplicates()
    # Merge with condition
    meta_cond = merged[['sample_name', 'osd_id', 'condition']].drop_duplicates()
    clock_samples = clock_samples.merge(meta_cond, on=['sample_name', 'osd_id'], how='left')
    flight_samples = clock_samples[clock_samples['condition'] == 'flight']['sample_name'].tolist()
    ground_samples = clock_samples[clock_samples['condition'] == 'ground_control']['sample_name'].tolist()
    # Take up to 100 each
    selected_samples = flight_samples[:100] + ground_samples[:100]

    clock_pivot = clock_expr[clock_expr['sample_name'].isin(selected_samples)].pivot_table(
        index='gene_name', columns='sample_name', values='expression'
    )
    # Sort: ground first, then flight
    sample_conditions = clock_samples.set_index('sample_name')['condition'].to_dict()
    sorted_samples = sorted(selected_samples, key=lambda s: (sample_conditions.get(s, ''), s))
    clock_pivot = clock_pivot[sorted_samples]

    heatmap_z = clock_pivot.values.tolist()
    heatmap_x = sorted_samples
    heatmap_y = list(clock_pivot.index)

    # Prepare t-SNE data
    tsne_data = []
    for _, row in tsne_df.iterrows():
        tsne_data.append({
            'sample_name': str(row['sample_name']),
            'osd_id': str(row['osd_id']),
            'condition': str(row.get('condition', '')),
            'tissue': str(row.get('tissue', '')),
            'tSNE1': float(row['tSNE1']),
            'tSNE2': float(row['tSNE2']),
            'predicted_CT': float(row['predicted_CT']) if pd.notna(row['predicted_CT']) else None,
        })

    # Trajectory analysis data
    traj_data = []
    if trajectory is not None:
        for _, row in trajectory.iterrows():
            traj_data.append({
                'osd_id': str(row['osd_id']),
                'tsne_centroid_distance': float(row['tsne_centroid_distance']),
                'phase_shift': float(row['phase_shift_hours']) if pd.notna(row['phase_shift_hours']) else None,
                'p_value': float(row['p_value']) if pd.notna(row['p_value']) else None,
            })

    # Studies list for dropdowns
    studies = sorted(merged['osd_id'].unique().tolist())
    tissues = sorted([t for t in merged['tissue'].dropna().unique().tolist() if t])
    light_regimes = sorted([l for l in merged['light_regime'].dropna().unique().tolist() if l])

    data_obj = {
        'scatter': scatter_data,
        'forest': forest_data,
        'pooled': {
            'effect': float(pooled.get('POOLED_EFFECT', 0)),
            'ci_lower': float(pooled.get('CI_LOWER', 0)),
            'ci_upper': float(pooled.get('CI_UPPER', 0)),
            'p_value': float(pooled.get('P_VALUE', 1)),
            'i2': float(pooled.get('I2', 0)),
            'n_studies': int(pooled.get('N_STUDIES', 0)),
        },
        'heatmap': {
            'z': heatmap_z,
            'x': heatmap_x,
            'y': heatmap_y,
            'sample_conditions': {s: sample_conditions.get(s, '') for s in heatmap_x},
        },
        'tsne': tsne_data,
        'trajectory': traj_data,
        'studies': studies,
        'tissues': tissues,
        'light_regimes': light_regimes,
    }

    # Write data.js
    js_content = "// Auto-generated dashboard data for Spaceflight Circadian Decoder\n"
    js_content += "// Generated from NASA GeneLab Arabidopsis transcriptomics analysis\n\n"
    js_content += "const DASHBOARD_DATA = " + json.dumps(data_obj, separators=(',', ':')) + ";\n"

    data_path = os.path.join(DOCS_DIR, 'data.js')
    with open(data_path, 'w') as f:
        f.write(js_content)
    print(f"  Saved data.js ({len(js_content):,} chars)")


def create_index_html():
    """Create the main dashboard HTML file."""
    print("Creating index.html...")

    html = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Spaceflight Circadian Decoder - Interactive Dashboard</title>
    <script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
    <script>window.Plotly || document.write('<script src="plotly.min.js"><\/script>')</script>
    <script src="data.js"></script>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: #FAF9F3;
            color: #000000;
            padding: 20px;
        }
        h1 {
            color: #0279EE;
            font-size: 24px;
            margin-bottom: 5px;
        }
        .subtitle {
            color: #666;
            font-size: 14px;
            margin-bottom: 20px;
        }
        .controls {
            background: #ECE9E2;
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 20px;
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
            align-items: center;
        }
        .controls label {
            font-size: 13px;
            font-weight: 600;
            color: #333;
        }
        .controls select {
            padding: 6px 12px;
            border: 1px solid #ccc;
            border-radius: 4px;
            font-size: 13px;
            background: white;
        }
        .tab-bar {
            display: flex;
            gap: 5px;
            margin-bottom: 15px;
            flex-wrap: wrap;
        }
        .tab {
            padding: 10px 20px;
            background: #ECE9E2;
            border: none;
            border-radius: 8px 8px 0 0;
            cursor: pointer;
            font-size: 14px;
            font-weight: 600;
            color: #666;
            transition: all 0.2s;
        }
        .tab.active {
            background: #0279EE;
            color: white;
        }
        .tab:hover:not(.active) {
            background: #ddd;
        }
        .panel {
            display: none;
            background: white;
            border-radius: 0 8px 8px 8px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        .panel.active { display: block; }
        .plot-container { width: 100%; }
        .stats-box {
            background: #ECE9E2;
            padding: 15px;
            border-radius: 8px;
            margin-top: 15px;
            font-size: 13px;
        }
        .stats-box strong { color: #0279EE; }
        .footer {
            margin-top: 30px;
            padding-top: 15px;
            border-top: 1px solid #ddd;
            font-size: 12px;
            color: #999;
        }
        #plot-phase, #plot-forest, #plot-fingerprint, #plot-clock, #plot-tsne {
            width: 100%;
        }
    </style>
</head>
<body>
    <h1>Spaceflight Circadian Decoder</h1>
    <p class="subtitle">Interactive exploration of circadian clock phase disruption in spaceflight-grown <em>Arabidopsis thaliana</em></p>

    <div class="tab-bar">
        <button class="tab active" onclick="showTab('phase', this)">Phase Overview</button>
        <button class="tab" onclick="showTab('forest', this)">Forest Plot</button>
        <button class="tab" onclick="showTab('fingerprint', this)">Circadian Fingerprint</button>
        <button class="tab" onclick="showTab('clock', this)">Clock Genes</button>
        <button class="tab" onclick="showTab('tsne', this)">t-SNE Trajectory</button>
    </div>

    <!-- Phase Overview Panel -->
    <div id="panel-phase" class="panel active">
        <div class="controls">
            <label>Study:</label>
            <select id="filter-study" onchange="updatePhasePlot()">
                <option value="">All studies</option>
            </select>
            <label>Tissue:</label>
            <select id="filter-tissue" onchange="updatePhasePlot()">
                <option value="">All tissues</option>
            </select>
            <label>Light regime:</label>
            <select id="filter-light" onchange="updatePhasePlot()">
                <option value="">All regimes</option>
            </select>
        </div>
        <div id="plot-phase" class="plot-container"></div>
        <div class="stats-box" id="phase-stats"></div>
    </div>

    <!-- Forest Plot Panel -->
    <div id="panel-forest" class="panel">
        <div id="plot-forest" class="plot-container"></div>
        <div class="stats-box" id="forest-stats"></div>
    </div>

    <!-- Fingerprint Panel -->
    <div id="panel-fingerprint" class="panel">
        <div id="plot-fingerprint" class="plot-container"></div>
    </div>

    <!-- Clock Genes Panel -->
    <div id="panel-clock" class="panel">
        <div id="plot-clock" class="plot-container"></div>
    </div>

    <!-- t-SNE Panel -->
    <div id="panel-tsne" class="panel">
        <div class="controls">
            <label>Color by:</label>
            <select id="tsne-colorby" onchange="updateTsnePlot()">
                <option value="condition">Condition</option>
                <option value="study">Study</option>
                <option value="ct">Predicted CT</option>
            </select>
        </div>
        <div id="plot-tsne" class="plot-container"></div>
    </div>

    <div class="footer">
        Data: NASA GeneLab (23 studies, 603 samples) | Model: ChronoGauge (100-model ensemble) |
        Analysis: Watson-Williams test + random-effects meta-analysis (metafor)
    </div>

    <script>
    // ============================================================
    // Tab switching
    // ============================================================
    function showTab(tabName, btnEl) {
        document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
        document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
        btnEl.classList.add('active');
        document.getElementById('panel-' + tabName).classList.add('active');

        // Trigger plot resize
        setTimeout(() => {
            if (tabName === 'phase') updatePhasePlot();
            if (tabName === 'forest') drawForestPlot();
            if (tabName === 'fingerprint') drawFingerprint();
            if (tabName === 'clock') drawClockHeatmap();
            if (tabName === 'tsne') updateTsnePlot();
        }, 50);
    }

    // ============================================================
    // Populate dropdowns
    // ============================================================
    (function initDropdowns() {
        const studySel = document.getElementById('filter-study');
        DASHBOARD_DATA.studies.forEach(s => {
            const opt = document.createElement('option');
            opt.value = s; opt.textContent = s;
            studySel.appendChild(opt);
        });
        const tissueSel = document.getElementById('filter-tissue');
        DASHBOARD_DATA.tissues.forEach(t => {
            const opt = document.createElement('option');
            opt.value = t; opt.textContent = t;
            tissueSel.appendChild(opt);
        });
        const lightSel = document.getElementById('filter-light');
        DASHBOARD_DATA.light_regimes.forEach(l => {
            const opt = document.createElement('option');
            opt.value = l; opt.textContent = l;
            lightSel.appendChild(opt);
        });
    })();

    // ============================================================
    // Panel 1: Phase Overview
    // ============================================================
    function updatePhasePlot() {
        const studyFilter = document.getElementById('filter-study').value;
        const tissueFilter = document.getElementById('filter-tissue').value;
        const lightFilter = document.getElementById('filter-light').value;

        let data = DASHBOARD_DATA.scatter.filter(d => {
            if (studyFilter && d.osd_id !== studyFilter) return false;
            if (tissueFilter && d.tissue !== tissueFilter) return false;
            if (lightFilter && d.light_regime !== lightFilter) return false;
            return d.predicted_CT !== null && d.circular_variance !== null;
        });

        const flight = data.filter(d => d.condition === 'flight');
        const ground = data.filter(d => d.condition === 'ground_control');

        const traces = [
            {
                x: ground.map(d => d.predicted_CT),
                y: ground.map(d => d.circular_variance),
                mode: 'markers', type: 'scatter',
                name: 'Ground control',
                marker: { color: '#0279EE', size: 8, opacity: 0.6,
                          line: { width: 0.5, color: 'white' } },
                text: ground.map(d => `${d.sample_name}<br>${d.osd_id}<br>${d.tissue}`),
                hoverinfo: 'text+x+y',
            },
            {
                x: flight.map(d => d.predicted_CT),
                y: flight.map(d => d.circular_variance),
                mode: 'markers', type: 'scatter',
                name: 'Spaceflight',
                marker: { color: '#FF9400', size: 8, opacity: 0.6,
                          line: { width: 0.5, color: 'white' } },
                text: flight.map(d => `${d.sample_name}<br>${d.osd_id}<br>${d.tissue}`),
                hoverinfo: 'text+x+y',
            }
        ];

        const layout = {
            xaxis: { title: 'Predicted circadian time (h)', range: [0, 24], dtick: 4 },
            yaxis: { title: 'Circular variance (1 - R)' },
            margin: { l: 60, r: 20, t: 30, b: 50 },
            legend: { x: 0.02, y: 0.98 },
            plot_bgcolor: '#FAF9F3',
            paper_bgcolor: 'white',
            height: 500,
        };

        Plotly.newPlot('plot-phase', traces, layout, { responsive: true });

        // Stats
        const nFlight = flight.length, nGround = ground.length;
        const meanVarF = nFlight > 0 ? (flight.reduce((s, d) => s + d.circular_variance, 0) / nFlight).toFixed(3) : 'N/A';
        const meanVarG = nGround > 0 ? (ground.reduce((s, d) => s + d.circular_variance, 0) / nGround).toFixed(3) : 'N/A';
        document.getElementById('phase-stats').innerHTML =
            `<strong>Samples:</strong> ${nGround} ground, ${nFlight} flight | ` +
            `<strong>Mean circular variance:</strong> Ground = ${meanVarG}, Flight = ${meanVarF}`;
    }

    // ============================================================
    // Panel 2: Forest Plot
    // ============================================================
    function drawForestPlot() {
        const included = DASHBOARD_DATA.forest.filter(d => d.included && d.phase_shift !== null);
        included.sort((a, b) => a.phase_shift - b.phase_shift);

        const pooled = DASHBOARD_DATA.pooled;
        const nStudies = included.length;

        const traces = [{
            x: included.map(d => d.phase_shift),
            y: included.map((d, i) => i),
            mode: 'markers', type: 'scatter',
            name: 'Phase shift',
            marker: {
                color: included.map(d => d.p_value < 0.05 ? '#FF9400' : '#0279EE'),
                size: 10, line: { width: 1, color: 'white' }
            },
            error_x: {
                type: 'data',
                symmetric: false,
                array: included.map(d => d.ci_upper - d.phase_shift),
                arrayminus: included.map(d => d.phase_shift - d.ci_lower),
                thickness: 1.5, width: 3,
            },
            text: included.map(d => `${d.osd_id}<br>n=${d.n_flight}F / ${d.n_ground}G<br>p=${d.p_value.toExponential(2)}<br>${d.tissue}, ${d.light_regime}`),
            hoverinfo: 'text',
        }];

        // Pooled estimate
        traces.push({
            x: [pooled.effect],
            y: [nStudies + 0.5],
            mode: 'markers', type: 'scatter',
            name: 'Pooled (REML)',
            marker: { color: '#75A025', size: 14, symbol: 'diamond' },
            error_x: {
                type: 'data', symmetric: false,
                array: [pooled.ci_upper - pooled.effect],
                arrayminus: [pooled.effect - pooled.ci_lower],
                thickness: 2, width: 5,
            },
            text: [`Pooled: ${pooled.effect.toFixed(3)}h (p=${pooled.p_value.toFixed(3)}, I2=${pooled.i2.toFixed(1)}%)`],
            hoverinfo: 'text',
        });

        // Zero line
        traces.push({
            x: [0, 0], y: [-0.5, nStudies + 1.5],
            mode: 'lines', type: 'scatter',
            name: 'No effect',
            line: { color: 'gray', dash: 'dash', width: 1 },
            showlegend: false,
        });

        const layout = {
            xaxis: { title: 'Phase shift (hours, flight - ground)' },
            yaxis: {
                tickvals: included.map((d, i) => i).concat([nStudies + 0.5]),
                ticktext: included.map(d => `${d.osd_id}`).concat(['Pooled']),
                autorange: 'reversed',
            },
            margin: { l: 120, r: 20, t: 30, b: 50 },
            legend: { x: 0.02, y: 0.98 },
            plot_bgcolor: '#FAF9F3',
            paper_bgcolor: 'white',
            height: 600,
        };

        Plotly.newPlot('plot-forest', traces, layout, { responsive: true });

        document.getElementById('forest-stats').innerHTML =
            `<strong>Pooled phase shift:</strong> ${pooled.effect.toFixed(3)}h ` +
            `(95% CI: ${pooled.ci_lower.toFixed(3)} to ${pooled.ci_upper.toFixed(3)}, ` +
            `p = ${pooled.p_value.toFixed(3)}, I2 = ${pooled.i2.toFixed(1)}%, k = ${pooled.n_studies} studies)`;
    }

    // ============================================================
    // Panel 3: Circadian Fingerprint
    // ============================================================
    function drawFingerprint() {
        const data = DASHBOARD_DATA.scatter.filter(d => d.circular_variance !== null);

        // Box plot: circular variance by condition
        const flight = data.filter(d => d.condition === 'flight').map(d => d.circular_variance);
        const ground = data.filter(d => d.condition === 'ground_control').map(d => d.circular_variance);

        const traces = [
            { y: ground, name: 'Ground control', type: 'box',
              marker: { color: '#0279EE' }, boxpoints: 'all', jitter: 0.3, pointpos: 0 },
            { y: flight, name: 'Spaceflight', type: 'box',
              marker: { color: '#FF9400' }, boxpoints: 'all', jitter: 0.3, pointpos: 0 },
        ];

        const layout = {
            yaxis: { title: 'Circular variance (1 - R)' },
            margin: { l: 60, r: 20, t: 30, b: 50 },
            plot_bgcolor: '#FAF9F3',
            paper_bgcolor: 'white',
            height: 500,
        };

        Plotly.newPlot('plot-fingerprint', traces, layout, { responsive: true });
    }

    // ============================================================
    // Panel 4: Clock Gene Heatmap
    // ============================================================
    function drawClockHeatmap() {
        const hm = DASHBOARD_DATA.heatmap;

        const traces = [{
            z: hm.z,
            x: hm.x.map(s => s.substring(0, 20)),
            y: hm.y,
            type: 'heatmap',
            colorscale: 'RdBu_r',
            zmin: -2, zmax: 2,
            colorbar: { title: 'Z-score' },
        }];

        const layout = {
            xaxis: { title: 'Samples', showticklabels: false },
            yaxis: { title: 'Clock gene', autorange: 'reversed' },
            margin: { l: 80, r: 60, t: 30, b: 50 },
            plot_bgcolor: '#FAF9F3',
            paper_bgcolor: 'white',
            height: 400,
        };

        Plotly.newPlot('plot-clock', traces, layout, { responsive: true });
    }

    // ============================================================
    // Panel 5: t-SNE Trajectory
    // ============================================================
    function updateTsnePlot() {
        const colorBy = document.getElementById('tsne-colorby').value;
        const data = DASHBOARD_DATA.tsne;

        let traces = [];

        if (colorBy === 'condition') {
            const flight = data.filter(d => d.condition === 'flight');
            const ground = data.filter(d => d.condition === 'ground_control');
            traces = [
                { x: ground.map(d => d.tSNE1), y: ground.map(d => d.tSNE2),
                  mode: 'markers', type: 'scatter', name: 'Ground control',
                  marker: { color: '#0279EE', size: 7, opacity: 0.6 },
                  text: ground.map(d => `${d.sample_name}<br>${d.osd_id}`),
                  hoverinfo: 'text' },
                { x: flight.map(d => d.tSNE1), y: flight.map(d => d.tSNE2),
                  mode: 'markers', type: 'scatter', name: 'Spaceflight',
                  marker: { color: '#FF9400', size: 7, opacity: 0.6 },
                  text: flight.map(d => `${d.sample_name}<br>${d.osd_id}`),
                  hoverinfo: 'text' },
            ];
        } else if (colorBy === 'study') {
            const studyColors = ['#0279EE','#FF9400','#75A025','#FD9BED','#E9ED4C','#000000',
                '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f',
                '#bcbd22','#17becf','#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#c49c94',
                '#f7b6d2','#dbdb8d','#9edae5','#393b79','#637939','#8c6d31','#843c39','#7b4173'];
            const studies = [...new Set(data.map(d => d.osd_id))].sort();
            studies.forEach((s, i) => {
                const subset = data.filter(d => d.osd_id === s);
                traces.push({
                    x: subset.map(d => d.tSNE1), y: subset.map(d => d.tSNE2),
                    mode: 'markers', type: 'scatter', name: s,
                    marker: { size: 7, opacity: 0.6,
                              color: studyColors[i % studyColors.length] },
                    text: subset.map(d => `${d.sample_name}<br>${d.osd_id}<br>${d.condition}`),
                    hoverinfo: 'text',
                });
            });
        } else { // ct
            traces = [{
                x: data.map(d => d.tSNE1), y: data.map(d => d.tSNE2),
                mode: 'markers', type: 'scatter', name: 'Predicted CT',
                marker: {
                    color: data.map(d => d.predicted_CT),
                    colorscale: 'Twilight', cmin: 0, cmax: 24,
                    size: 7, opacity: 0.6,
                    colorbar: { title: 'CT (h)' },
                },
                text: data.map(d => `${d.sample_name}<br>${d.osd_id}<br>CT=${d.predicted_CT ? d.predicted_CT.toFixed(2) : 'N/A'}h`),
                hoverinfo: 'text',
            }];
        }

        const layout = {
            xaxis: { title: 't-SNE 1' },
            yaxis: { title: 't-SNE 2' },
            margin: { l: 60, r: 20, t: 30, b: 50 },
            legend: { x: 1.02, y: 0.98, font: { size: 8 } },
            plot_bgcolor: '#FAF9F3',
            paper_bgcolor: 'white',
            height: 550,
        };

        Plotly.newPlot('plot-tsne', traces, layout, { responsive: true });
    }

    // ============================================================
    // Initialize all plots on load
    // ============================================================
    updatePhasePlot();
    </script>
</body>
</html>'''

    html_path = os.path.join(DOCS_DIR, 'index.html')
    with open(html_path, 'w') as f:
        f.write(html)
    print(f"  Saved index.html ({len(html):,} chars)")


def create_docs_readme():
    """Create deployment instructions for GitHub Pages."""
    readme = '''# Interactive Dashboard - GitHub Pages Deployment

## Overview

This folder contains a self-contained interactive dashboard for exploring the
spaceflight circadian decoder results. It requires no server-side processing —
just static HTML, JavaScript, and data files.

## Files

- `index.html` — Main dashboard with 5 interactive views
- `data.js` — All analysis data embedded as JSON
- `plotly.min.js` — Plotly.js library (local copy, no CDN dependency)

## Dashboard views

1. **Phase Overview** — Predicted circadian time vs circular variance, with
   dropdown filters for study, tissue, and light regime
2. **Forest Plot** — Per-study phase shifts with confidence intervals and
   pooled meta-analysis estimate
3. **Circadian Fingerprint** — Box plot of circular variance by condition
4. **Clock Genes** — Heatmap of 11 core clock genes across samples
5. **t-SNE Trajectory** — Non-linear embedding of circadian fingerprints,
   colorable by condition, study, or predicted CT

## Deploy to GitHub Pages

### Option 1: Using the dashboard/ folder

1. Push this repository to GitHub
2. Go to **Settings** > **Pages**
3. Under **Source**, select **Deploy from a branch**
4. Select **main** branch and **/dashboard** folder (or root and access via `/dashboard/index.html`)
5. Click **Save**
6. Your dashboard will be available at:
   `https://<username>.github.io/Circadian_decoder/dashboard/`

### Option 2: Using GitHub Actions (automatic deployment)

Add this workflow to `.github/workflows/deploy-dashboard.yml`:

```yaml
name: Deploy Dashboard
on:
  push:
    branches: [main]
    paths: ['dashboard/**']
jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: peaceiris/actions-gh-pages@v3
        with:
          github_token: ${{ secrets.GITHUB_TOKEN }}
          publish_dir: ./dashboard
```

## Local preview

```bash
cd dashboard
python3 -m http.server 8000
# Open http://localhost:8000 in your browser
```

## Data sources

- NASA GeneLab (23 Arabidopsis spaceflight transcriptomics studies)
- ChronoGauge (100-model ensemble circadian time predictor)
- Watson-Williams test + random-effects meta-analysis (metafor)
'''

    readme_path = os.path.join(DOCS_DIR, 'README.md')
    with open(readme_path, 'w') as f:
        f.write(readme)
    print(f"  Saved dashboard/README.md")


def main():
    print("=" * 60)
    print("INTERACTIVE DASHBOARD GENERATION")
    print("=" * 60)

    # Download plotly.js
    download_plotly_js()

    # Prepare data
    merged, per_study, meta_results, clock_expr, trajectory, tsne_df = prepare_data()

    # Create data.js
    create_dashboard_data_js(merged, per_study, meta_results, clock_expr, trajectory, tsne_df)

    # Create index.html
    create_index_html()

    # Create deployment README
    create_docs_readme()

    print("\n" + "=" * 60)
    print("DASHBOARD GENERATION COMPLETE")
    print("=" * 60)
    print(f"Output: {DOCS_DIR}")
    for f in sorted(os.listdir(DOCS_DIR)):
        size = os.path.getsize(os.path.join(DOCS_DIR, f))
        print(f"  {f}: {size:,} bytes")


if __name__ == "__main__":
    main()
