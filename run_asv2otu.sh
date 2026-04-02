#!/bin/bash
# =============================================================================
# Pipeline submission wrapper
# Reads email from config.sh and submits the SLURM job with --mail-user set.
# Symlinked to ~/bin/run_asv2otu by setup.sh — run from anywhere.
#
# USAGE:
#   run_asv2otu <project_dir> <asv_table_path> <proj_name> <primer_set> [--skip-itsx] [--db <database>]
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/config.sh"

[ -f "$CONFIG" ] || { echo "ERROR: config.sh not found. Run setup.sh first."; exit 1; }
source "$CONFIG"
[ -n "${SLURM_EMAIL:-}" ] || { echo "ERROR: SLURM_EMAIL not set in config.sh"; exit 1; }

sbatch --mail-user="$SLURM_EMAIL" "$SCRIPT_DIR/ASVtoOTU_msiSLURM.sh" "$@"
