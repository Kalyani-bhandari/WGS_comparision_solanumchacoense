#!/bin/bash
set -euo pipefail

source ~/.bashrc
conda activate wgs_env
source config/paths.sh

echo "Building SnpEff database for Solanum chacoense M6..."

mkdir -p "$SNPEFF_DATA"

# Copy reference files
cp "$REF" "$SNPEFF_DATA/sequences.fa"
cp "$GFF" "$SNPEFF_DATA/genes.gff"

# Register genome in config
CONFIG="$SNPEFF_DIR/snpEff.config"
grep -q "solanum_chacoense_m6.genome" "$CONFIG" || \
  echo "solanum_chacoense_m6.genome : Solanum_chacoense_M6" >> "$CONFIG"

# Build database
cd "$SNPEFF_DIR"
java -Xmx8g -jar snpEff.jar build -gff3 -v solanum_chacoense_m6

echo "SnpEff database build complete."

