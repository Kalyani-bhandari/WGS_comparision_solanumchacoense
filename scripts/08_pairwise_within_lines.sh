#!/bin/bash
set -euo pipefail

# ================================
# Environment
# ================================
source ~/.bashrc
conda activate wgs_env

# ================================
# Directories
# ================================
EMSDIR="results/EMS_vcfs"
OUTDIR="results/compare_within_lines"

mkdir -p "$OUTDIR"

echo "Running pairwise comparisons within EMS lines..."

# ================================
# Extract unique line IDs (e.g., 4, 5)
# ================================
LINES=$(ls $EMSDIR/*.EMS.vcf.gz \
  | sed 's/.*\///' \
  | sed 's/-Rest.EMS.vcf.gz//' \
  | sed 's/-Susp.EMS.vcf.gz//' \
  | sort -u)

for LINE in $LINES; do
  REST="$EMSDIR/${LINE}-Rest.EMS.vcf.gz"
  SUSP="$EMSDIR/${LINE}-Susp.EMS.vcf.gz"

  echo "--------------------------------------------"
  echo "Line: $LINE"
  echo "Resistant: $REST"
  echo "Susceptible: $SUSP"

  if [[ ! -f "$REST" || ! -f "$SUSP" ]]; then
    echo "Skipping line $LINE (missing files)"
    continue
  fi

  bcftools isec \
    -p "$OUTDIR/line_${LINE}" \
    "$REST" "$SUSP"

  # Compress and index outputs
  for f in "$OUTDIR/line_${LINE}"/*.vcf; do
    bgzip "$f"
    bcftools index "$f.gz"
  done

done

echo "Within-line comparisons completed."

