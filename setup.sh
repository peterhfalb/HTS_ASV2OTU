#!/bin/bash
# =============================================================================
# ASV to OTU Pipeline Setup Script
# Run once after cloning to configure the pipeline for your account
# =============================================================================

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo ""
echo "============================================================"
echo "  Amplicon Pipeline Setup"
echo "============================================================"
echo "  Repo directory: $REPO_DIR"
echo ""

# ------------------------------------------------------------------------------
# 1. Set PIPELINE_DIR in the SLURM script
# ------------------------------------------------------------------------------

echo "--- Configuring pipeline directory ---"
sed -i "s|PIPELINE_DIR=.*|PIPELINE_DIR=\"$REPO_DIR\"|" "$REPO_DIR/ASVtoOTU_msiSLURM.sh"
sed -i "s|PIPELINE_DIR=.*|PIPELINE_DIR=\"$REPO_DIR\"|" "$REPO_DIR/ASVtoOTU_msiSLURM.sh" 2>/dev/null || true
echo "  PIPELINE_DIR set to: $REPO_DIR"

# ------------------------------------------------------------------------------
# 2. Set email in SLURM script
# ------------------------------------------------------------------------------

echo ""
echo "--- Configuring SLURM email ---"
read -p "  Enter your email address for SLURM notifications: " USER_EMAIL
sed -i "s|#SBATCH --mail-user=.*|#SBATCH --mail-user=$USER_EMAIL|" "$REPO_DIR/ASVtoOTU_msiSLURM.sh"
echo "  Email set to: $USER_EMAIL"

# ------------------------------------------------------------------------------
# 3. Create conda environment
# ------------------------------------------------------------------------------

echo ""
echo "--- Setting up conda environment ---"

# Source conda
CONDA_BASE=$(conda info --base 2>/dev/null) || { echo "ERROR: conda not found. Please load/install conda first."; exit 1; }
source "$CONDA_BASE/etc/profile.d/conda.sh"

if conda env list | grep -q "itsx_env"; then
    echo "  itsx_env already exists, skipping creation."
else
    echo "  Creating itsx_env (this may take several minutes)..."
    conda create -n itsx_env -c bioconda -c conda-forge \
        itsx vsearch blast bioconductor-dada2 \
        r-dplyr r-tidyr r-readr r-tibble r-purrr r-stringr -y
    echo "  itsx_env created successfully."
fi

# Activate and configure
conda activate itsx_env
export PATH="$HOME/.conda/envs/itsx_env/bin:$PATH"
unset R_HOME 2>/dev/null || true

# ------------------------------------------------------------------------------
# 4. Validate R packages
# ------------------------------------------------------------------------------

echo ""
echo "--- Validating R packages ---"

CONDA_R="$HOME/.conda/envs/itsx_env/bin/Rscript"

if [ ! -f "$CONDA_R" ]; then
    echo "ERROR: conda Rscript not found at $CONDA_R"
    exit 1
fi

PACKAGES="dada2 dplyr tidyr readr tibble purrr stringr"
ALL_OK=true

for pkg in $PACKAGES; do
    if $CONDA_R -e "library($pkg)" 2>/dev/null; then
        echo "  $pkg: OK"
    else
        echo "  $pkg: MISSING"
        ALL_OK=false
    fi
done

if [ "$ALL_OK" = false ]; then
    echo ""
    echo "WARNING: Some packages are missing. Try running:"
    echo "  conda install -n itsx_env -c conda-forge -c bioconda <package_name>"
    exit 1
fi

# ------------------------------------------------------------------------------
# 5. Compile mumu
# ------------------------------------------------------------------------------

echo ""
echo "--- Setting up mumu ---"

MUMU_BIN="$HOME/packages/mumu/mumu"

if [ -f "$MUMU_BIN" ]; then
    echo "  mumu binary already exists, skipping compilation."
else
    echo "  Compiling mumu..."
    # Check if gcc/13.1.0 is available via module
    if module load gcc/13.1.0-5z64cho 2>/dev/null; then
        echo "  Loaded gcc/13.1.0"
    else
        echo "  WARNING: Could not load gcc/13.1.0-5z64cho — will try system gcc"
    fi

    mkdir -p "$HOME/packages"
    if [ ! -d "$HOME/packages/mumu" ]; then
        git clone https://github.com/frederic-mahe/mumu.git "$HOME/packages/mumu"
    fi
    (cd "$HOME/packages/mumu" && make && make check)
    chmod +x "$MUMU_BIN"
    echo "  mumu compiled successfully."
fi

# Validate mumu
if "$MUMU_BIN" --version 2>/dev/null || "$MUMU_BIN" --help 2>/dev/null | head -1; then
    echo "  mumu: OK"
else
    echo "  WARNING: mumu binary exists but may not be functional. Test manually:"
    echo "  $MUMU_BIN --help"
fi

# ------------------------------------------------------------------------------
# 6. Check databases
# ------------------------------------------------------------------------------

echo ""
echo "--- Checking databases ---"

SHARED_DB="/projects/standard/kennedyp/shared/taxonomy"

declare -A DB_FILES
DB_FILES["ITS1/ITS2"]="sh_general_release_dynamic_all_19.02.2025.fasta"
DB_FILES["16S-V4"]="silva_nr99_v138.2_toGenus_trainset.fa"
DB_FILES["18S-V4"]="pr2_version_5.1.1_SSU_dada2.fasta"
DB_FILES["18S-AMF"]="maarjam_dada2.txt"

ALL_DB_OK=true
for primer in "${!DB_FILES[@]}"; do
    db="${DB_FILES[$primer]}"
    if [ -f "$SHARED_DB/$db" ]; then
        echo "  $primer: OK"
    else
        echo "  $primer: NOT FOUND at $SHARED_DB/$db"
        ALL_DB_OK=false
    fi
done

if [ "$ALL_DB_OK" = false ]; then
    echo ""
    echo "WARNING: Some databases not found. Contact Peter F (falb0011@umn.edu) or check shared directory:"
    echo "  $SHARED_DB"
fi

# ------------------------------------------------------------------------------
# Done
# ------------------------------------------------------------------------------

echo ""
echo "============================================================"
echo "  Setup complete!"
echo ""
echo "  Run the pipeline with:"
echo "  sbatch ASVtoOTU_msiSLURM.sh <project_dir> <asv_table_path> <proj_name> <primer_set> [--skip-itsx]"
echo ""
echo "  Primer Set Options:"
echo "  ITS1    - fungal ITS1 region (UNITE database, ITSx optional)"
echo "  ITS2    - fungal ITS2 region (UNITE database, ITSx optional)"
echo "  16S-V4  - bacterial 16S V4 region (SILVA database, no ITSx)"
echo "  18S-V4  - microeukaryote 18S V4 region (PR2 database, no ITSx)"
echo "  18S-AMF - arbuscular mycorrhizal fungi 18S (MaarjAM database, no ITSx)"
echo ""
echo "  Example:"
echo "  sbatch ASVtoOTU_msiSLURM.sh /path/to/project /path/to/ASVtable.tsv FAB2 ITS2"
echo ""
echo "  NOTE: You need to be in the HTS_ASV2OTU/ directory (which contains the ASVtoOTU_msiSLURM.sh script) to run the sbatch command."
echo " "
echo "============================================================"
echo ""