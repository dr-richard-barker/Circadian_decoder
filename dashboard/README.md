# Interactive Dashboard - GitHub Pages Deployment

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
