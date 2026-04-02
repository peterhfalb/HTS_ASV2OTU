#!/bin/bash
# =============================================================================
# Pipeline submission wrapper
# Reads email from config.sh and submits the SLURM job with --mail-user set.
# Symlinked to ~/bin/run_asv2otu by setup.sh — run from anywhere.
#
# USAGE:
#   run_asv2otu <project_dir> <asv_table_path> <proj_name> <primer_set> [--skip-itsx] [--db <database>]
# =============================================================================

SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
CONFIG="$SCRIPT_DIR/config.sh"

[ -f "$CONFIG" ] || { echo "ERROR: config.sh not found. Run setup.sh first."; exit 1; }
source "$CONFIG"
[ -n "${SLURM_EMAIL:-}" ] || { echo "ERROR: SLURM_EMAIL not set in config.sh"; exit 1; }

PROJECT_DIR="${1:-}"
[ -n "$PROJECT_DIR" ] || { echo "ERROR: project_dir argument is required."; echo "Usage: run_asv2otu <project_dir> <asv_table> <proj_name> <primer_set> [--skip-itsx] [--db <database>]"; exit 1; }
[ -d "$PROJECT_DIR" ] || { echo "ERROR: project directory not found: $PROJECT_DIR"; exit 1; }

sbatch --mail-user="$SLURM_EMAIL" \
    --output="$PROJECT_DIR/pipeline_%j.out" \
    --error="$PROJECT_DIR/pipeline_%j.err" \
    --export=ALL,PIPELINE_DIR="$SCRIPT_DIR" \
    "$SCRIPT_DIR/ASVtoOTU_msiSLURM.sh" "$@"
