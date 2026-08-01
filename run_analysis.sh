#!/bin/bash
# Run the Spaceflight Circadian Decoder analysis pipeline

echo "============================================================"
echo "Starting Spaceflight Circadian Decoder analysis pipeline..."
echo "============================================================"

# Run Python analysis pipeline
python3 src/run_analysis.py

# Run R-based differential expression & GSEA pathway analysis
if command -v Rscript >/dev/null 2>&1; then
  echo ""
  echo "Running Gene Set Enrichment Analysis (R)..."
  Rscript src/enrichment_analysis.R
else
  echo ""
  echo "WARNING: Rscript not found. Skipping Gene Set Enrichment Analysis."
  echo "Pre-generated GSEA results are preserved in tables/tableS3_enrichment_results.csv."
fi

# Run interactive dashboard generator
python3 src/generate_dashboard.py

echo ""
echo "============================================================"
echo "Analysis pipeline completed successfully!"
echo "Figures saved to: ./figures/"
echo "Tables saved to: ./tables/"
echo "Data saved to: ./data/"
echo "Interactive dashboard updated in: ./dashboard/"
echo "============================================================"
