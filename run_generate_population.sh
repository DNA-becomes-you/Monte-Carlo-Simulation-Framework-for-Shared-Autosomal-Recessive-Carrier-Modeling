#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default parameters
USERS=${1:-2000}
GENES=${2:-500}
AVG_CARRIERS=${3:-3}
VARIANCE=${4:-1.5}
STRUCTURE=${5:-0.02}
SEED=${6:-42}

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

OUTPUT="$SCRIPT_DIR/output/population_${USERS}u_${GENES}g_${AVG_CARRIERS}c_struct${STRUCTURE}_${TIMESTAMP}.json"

mkdir -p "$SCRIPT_DIR/output"

echo "Running DNABY population generator..."
echo "Users: $USERS"
echo "Genes: $GENES"
echo "Avg carriers: $AVG_CARRIERS"
echo "Structure: $STRUCTURE"
echo "Seed: $SEED"
echo "Output: $OUTPUT"

python3 "$SCRIPT_DIR/generate_population.py" \
  --users "$USERS" \
  --genes "$GENES" \
  --avg-carriers "$AVG_CARRIERS" \
  --variance "$VARIANCE" \
  --population-structure "$STRUCTURE" \
  --seed "$SEED" \
  --output "$OUTPUT"

echo "Done."
