#!/bin/bash
# =============================================================================
# Recover synthetic mock community sequences filtered out by ITSx
#
# Purpose: ITSx may remove synthetic mock community sequences as it cannot
# identify them as natural ITS regions. This script recovers them from
# ITSx's "no_detections" output and appends them back to the main ITSx output.
#
# Usage: recover_synmock_from_itsx.sh <synmock_fasta> <itsx_output_fasta> \
#                                      <itsx_no_detections> <output_fasta>
#
# Args:
#   synmock_fasta:      Path to synmock.fasta (contains synmock_1, synmock_2, ...)
#   itsx_output_fasta:  Path to main ITSx output (e.g., Centroid.ITSx.ITS1.fasta)
#   itsx_no_detections: Path to ITSx no_detections output (sequences not matched)
#   output_fasta:       Path to write combined output (itsx_output + recovered synmock)
# =============================================================================

set -euo pipefail

SYNMOCK_FASTA="${1:-}"
ITSX_OUTPUT="${2:-}"
NO_DETECTIONS="${3:-}"
OUTPUT_FASTA="${4:-}"

[ -n "$SYNMOCK_FASTA" ] || { echo "ERROR: synmock_fasta argument required"; exit 1; }
[ -n "$ITSX_OUTPUT" ] || { echo "ERROR: itsx_output_fasta argument required"; exit 1; }
[ -n "$NO_DETECTIONS" ] || { echo "ERROR: itsx_no_detections argument required"; exit 1; }
[ -n "$OUTPUT_FASTA" ] || { echo "ERROR: output_fasta argument required"; exit 1; }

[ -f "$SYNMOCK_FASTA" ] || { echo "ERROR: synmock_fasta not found: $SYNMOCK_FASTA"; exit 1; }
[ -f "$ITSX_OUTPUT" ] || { echo "ERROR: itsx_output_fasta not found: $ITSX_OUTPUT"; exit 1; }
[ -f "$NO_DETECTIONS" ] || { echo "ERROR: no_detections file not found: $NO_DETECTIONS"; exit 1; }

# Extract synmock sequence IDs from synmock.fasta
SYNMOCK_IDS=$(grep "^>" "$SYNMOCK_FASTA" | sed 's/>//g' | sort)

echo "Checking for synmock sequences in ITSx output..."

# Create output directory if needed
mkdir -p "$(dirname "$OUTPUT_FASTA")"

# Start with main ITSx output
cp "$ITSX_OUTPUT" "$OUTPUT_FASTA"

RECOVERED_COUNT=0

# Check if each synmock ID is in the ITSx output
for id in $SYNMOCK_IDS; do
    if grep -q "^>$id" "$ITSX_OUTPUT"; then
        # Already in ITSx output, skip
        :
    elif grep -q "^>$id" "$NO_DETECTIONS"; then
        # Found in no_detections, recover it
        echo "  Recovering: $id"
        # Extract the sequence from no_detections and append to output
        awk -v target_id="$id" '
          $0 ~ "^>" target_id { found=1; print; next }
          found && /^>/ { found=0 }
          found { print }
        ' "$NO_DETECTIONS" >> "$OUTPUT_FASTA"
        RECOVERED_COUNT=$((RECOVERED_COUNT + 1))
    fi
done

if [ "$RECOVERED_COUNT" -gt 0 ]; then
    echo "Recovered $RECOVERED_COUNT synmock sequences from ITSx no_detections"
else
    echo "No synmock sequences needed recovery (all found in main ITSx output)"
fi

echo "Mock community recovery complete. Output: $OUTPUT_FASTA"
