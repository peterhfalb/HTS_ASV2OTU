#!/bin/bash

# ==============================================================================
# ASV to OTU Pipeline - SLURM Cluster Version (UMN MSI)
# ==============================================================================
#
# USAGE: Navigate to the HTS_ASV2OTU directory then run the command
#   run_asv2otu <project_dir> <asv_table_path> <proj_name> <primer_set> [--skip-itsx] [--db <database>]
#
# PRIMER SET OPTIONS:
#   ITS1       — fungal ITS1 region (UNITE database, ITSx optional)
#   ITS2       — fungal ITS2 region (UNITE database, ITSx optional)
#   16S-V4     — bacterial 16S V4 region (SILVA database, no ITSx)
#   18S-V4     — microeukaryote 18S V4 region (PR2 database, no ITSx)
#   18S-AMF    — arbuscular mycorrhizal fungi 18S (MaarjAM database, no ITSx)
#
# DATABASE OPTIONS (if manually chosen, using flag --db):
#   SILVA        - bacteria SSU 16S rRNA sequences
#   UNITE        - fungal ITS1 and ITS2 regions
#   PR2          - full-length SSU 18S sequences, covers whole eukaryote tree but with focus on protists
#   Maarjam      - AMF 18S SSU sequences
#   EukaryomeITS - ITS1 and ITS2 sequences with good coverage across the eukaryote tree
#   EukaryomeSSU - 18S SSU sequences with good coverage across the eukaryote tree (especially for AMF?)
#
# EXAMPLES:
#   run_asv2otu /path/to/project /path/to/table.tsv FAB2 ITS2
#   run_asv2otu /path/to/project /path/to/table.tsv FAB2 ITS2 --skip-itsx
#   run_asv2otu /path/to/project /path/to/table.tsv FAB2 16S-V4
#
# NOTE: Use run_asv2otu to submit — it reads your email from config.sh
#       and passes it to sbatch automatically.
#
# ==============================================================================

#SBATCH --job-name=ASVtoOTUpipeline
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=5G
#SBATCH -p msismall
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=pipeline_%j.out
#SBATCH --error=pipeline_%j.err

set -euo pipefail

# ------------------------------------------------------------------------------
# Load user configuration (PIPELINE_DIR, SLURM_EMAIL set by setup.sh)
# ------------------------------------------------------------------------------

[ -n "${PIPELINE_DIR:-}" ] || { echo "ERROR: PIPELINE_DIR not set. Submit via run_asv2otu, not sbatch directly."; exit 1; }
CONFIG="$PIPELINE_DIR/config.sh"
[ -f "$CONFIG" ] || { echo "ERROR: config.sh not found at $PIPELINE_DIR. Run setup.sh first."; exit 1; }
source "$CONFIG"

# ------------------------------------------------------------------------------
# Load modules
# ------------------------------------------------------------------------------

module purge
module load conda
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh

module load R/4.2.2-gcc-8.2.0-vp7tyde
module load blast-plus/2.13.0-gcc-8.2.0-vo4mr4d

# Add local packages (mumu binary) to PATH
export PATH=$HOME/packages:$PATH

# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------

if [ "$#" -lt 4 ]; then
    echo "Usage: run_asv2otu <project_dir> <asv_table> <proj_name> <primer_set> [--run-itsx] [--db <database>]"
    echo ""
    echo "  primer_set options: ITS1, ITS2, 16S-V4, 18S-V4, 18S-AMF"
    echo "  --db options:       UNITE, SILVA, PR2, Maarjam, EukaryomeITS, EukaryomeSSU"
    echo "  --run-itsx:         enable ITSx extraction (ITS1/ITS2 only; disabled by default)"
    exit 1
fi

PROJECT_DIR="$1"
ASV_TABLE="$2"
PROJ="$3"
PRIMER_SET="$4"

# Parse optional flags
RUN_ITSX_FLAG=false
DB_OVERRIDE=""

for arg in "${@:5}"; do
    case "$arg" in
        --run-itsx)
            RUN_ITSX_FLAG=true
            ;;
        --db)
            # handled by next iteration — see below
            ;;
        UNITE|SILVA|PR2|Maarjam|EukaryomeITS|EukaryomeSSU)
            DB_OVERRIDE="$arg"
            ;;
        *)
            echo "ERROR: Unrecognized argument: '$arg'"
            echo "       Valid optional arguments: --run-itsx, --db <UNITE|SILVA|PR2|Maarjam|EukaryomeITS|EukaryomeSSU>"
            exit 1
            ;;
    esac
done

[ -d "$PROJECT_DIR" ] || { echo "ERROR: project directory not found: $PROJECT_DIR"; exit 1; }
[ -f "$ASV_TABLE" ]   || { echo "ERROR: ASV table not found: $ASV_TABLE"; exit 1; }

cd "$PROJECT_DIR"

# ------------------------------------------------------------------------------
# Output subdirectories
# ------------------------------------------------------------------------------

DIR_INPUT="$PROJECT_DIR/01_input"
DIR_ITSX="$PROJECT_DIR/02_itsx"
DIR_VSEARCH="$PROJECT_DIR/03_vsearch"
DIR_MUMU="$PROJECT_DIR/04_mumu"
DIR_TAXONOMY="$PROJECT_DIR/05_taxonomy"

mkdir -p "$DIR_INPUT" "$DIR_VSEARCH" "$DIR_MUMU" "$DIR_TAXONOMY"
# DIR_ITSX created only if ITSx runs

# ------------------------------------------------------------------------------
# Set database — use override if provided, otherwise use primer set default
# ------------------------------------------------------------------------------

DB_DIR="/projects/standard/kennedyp/shared/taxonomy"

# Default databases per primer set
case "$PRIMER_SET" in
    ITS1|ITS2)
        DEFAULT_DB="UNITE"
        RUN_ITSX=false  # ITSx is disabled by default; use --run-itsx to enable
        ITSX_REGION="$PRIMER_SET"
        ;;
    16S-V4)
        DEFAULT_DB="SILVA"
        RUN_ITSX=false
        ;;
    18S-V4)
        DEFAULT_DB="PR2"
        RUN_ITSX=false
        ;;
    18S-AMF)
        DEFAULT_DB="Maarjam"
        RUN_ITSX=false
        ;;
    *)
        echo "ERROR: Unrecognized primer set: '$PRIMER_SET'"
        echo "       Valid options: ITS1, ITS2, 16S-V4, 18S-V4, 18S-AMF"
        exit 1
        ;;
esac

# Apply override if specified
DB_NAME="${DB_OVERRIDE:-$DEFAULT_DB}"

# Resolve database name to file path
case "$DB_NAME" in
    UNITE)
        DB="$DB_DIR/sh_general_release_dynamic_all_19.02.2025_wKennedySynmock.fasta"
        ;;
    SILVA)
        DB="$DB_DIR/silva_nr99_v138.1_train_set.fa"
        ;;
    PR2)
        DB="$DB_DIR/pr2_version_5.1.1_SSU_dada2.fasta"
        ;;
    Maarjam)
        DB="$DB_DIR/maarjam_dada2_SSU_2021.fasta"
        ;;
    EukaryomeITS)
        DB="$DB_DIR/DADA2_EUK_ITS_v2.0_wKennedySynmock.fasta"
        ;;
    EukaryomeSSU)
        DB="$DB_DIR/DADA2_EUK_SSU_v2.0.fasta"
        ;;
esac

[ -f "$DB" ] || { echo "ERROR: Database file not found: $DB"; exit 1; }

if [ "$RUN_ITSX_FLAG" = true ]; then
    RUN_ITSX=true
fi

# ------------------------------------------------------------------------------
# Logging setup
# ------------------------------------------------------------------------------

LOG="$PROJECT_DIR/pipeline_run.log"
> "$LOG"

plog() {
    echo "$1" | tee -a "$LOG"
}

plog_section() {
    plog ""
    plog "============================================================"
    plog "  $1"
    plog "============================================================"
}

plog_file() {
    local filepath="$1"
    local description="$2"
    if [ -f "$filepath" ]; then
        local size
        size=$(du -sh "$filepath" | cut -f1)
        plog "  [FILE] $filepath ($size)"
        plog "         $description"
    else
        plog "  [MISSING] $filepath — $description"
    fi
}

plog "ASV to OTU Pipeline Log — UMN MSI Cluster Version"
plog "Run started:   $(date)"
plog "SLURM job ID:  ${SLURM_JOB_ID:-N/A}"
plog "Project dir:   $PROJECT_DIR"
plog "ASV table:     $ASV_TABLE"
plog "Project name:  $PROJ"
plog "Primer set:    $PRIMER_SET"
plog "Database:      $DB ($DB_NAME)"
if [ "$RUN_ITSX" = true ]; then
    plog "ITSx:          enabled (region: $ITSX_REGION)"
else
    if [[ "$PRIMER_SET" == ITS* ]] && [ "$SKIP_ITSX" = true ]; then
        plog "ITSx:          skipped (--skip-itsx flag)"
    else
        plog "ITSx:          not applicable for $PRIMER_SET"
    fi
fi

# ------------------------------------------------------------------------------
# Detect available threads from SLURM allocation
# ------------------------------------------------------------------------------

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
    THREADS=$SLURM_CPUS_PER_TASK
elif command -v nproc &>/dev/null; then
    THREADS=$(nproc)
else
    THREADS=8
    plog "  WARNING: Could not detect CPU count, defaulting to $THREADS threads"
fi

plog "  Using $THREADS threads (SLURM allocation)"

# ------------------------------------------------------------------------------
# Pre-flight checks
# ------------------------------------------------------------------------------

[ -f "$DB" ]                                                      || { echo "ERROR: database not found at $DB"; exit 1; }
[ -f "$PIPELINE_DIR/back_ground_scripts/01_prepare_input.R" ]       || { echo "ERROR: 01_prepare_input.R not found"; exit 1; }
[ -f "$PIPELINE_DIR/back_ground_scripts/03_build_otu_table.R" ]     || { echo "ERROR: 03_build_otu_table.R not found"; exit 1; }
[ -f "$PIPELINE_DIR/back_ground_scripts/05_assign_taxonomy_rdp.R" ] || { echo "ERROR: 05_assign_taxonomy_rdp.R not found"; exit 1; }
[ -f "$PIPELINE_DIR/back_ground_scripts/06_combine_otu_taxonomy.R" ] || { echo "ERROR: 06_combine_otu_taxonomy.R not found"; exit 1; }

# ------------------------------------------------------------------------------
# Step 0: Create conda virtual environment (skip if already exists)
# ------------------------------------------------------------------------------

plog_section "Step 0: Conda Environment"

if ! conda env list | grep -q "itsx_env"; then
    plog "  Creating itsx_env..."
    conda create -n itsx_env -c bioconda -c conda-forge \
        itsx vsearch blast bioconductor-dada2 \
        r-dplyr r-tidyr r-readr r-tibble r-purrr r-stringr -y
else
    plog "  itsx_env already exists, skipping."
fi

# Activate conda env and ensure it takes priority over system R
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate itsx_env
export PATH="$HOME/.conda/envs/itsx_env/bin:$PATH"
unset R_HOME 2>/dev/null || true

# ------------------------------------------------------------------------------
# Check R packages
# ------------------------------------------------------------------------------

plog_section "Checking R package dependencies"
CONDA_R="$HOME/.conda/envs/itsx_env/bin/Rscript"
if [ -f "$CONDA_R" ]; then
    plog "  Using conda R from itsx_env"
    RSCRIPT="$CONDA_R"
    $RSCRIPT -e "library(dada2)" 2>/dev/null && plog "  dada2: OK" || \
        { plog "ERROR: dada2 not found in itsx_env conda R"; exit 1; }
else
    plog "  Using system R, running setup..."
    RSCRIPT="Rscript"
    Rscript "$PIPELINE_DIR/back_ground_scripts/setup_r_packages.R" 2>&1 | tee -a "$LOG"
fi

# ------------------------------------------------------------------------------
# Step 1: Prepare input
# ------------------------------------------------------------------------------

plog_section "Step 1: Prepare Input FASTA and Abundance Table"

$RSCRIPT "$PIPELINE_DIR/back_ground_scripts/01_prepare_input.R" \
    "$ASV_TABLE" \
    "$DIR_INPUT" 2>&1 | tee -a "$LOG"

cat > "$DIR_INPUT/README.txt" << 'README'
Input Preparation Files
=======================
Generated by: 01_prepare_input.R

Centroid.fasta
  All input ASVs as FASTA with >SeqXXXXXX;size=<total_reads> headers.
  The size= annotation is required by VSEARCH for abundance-weighted clustering.

asv_abundance_table.txt
  ASV abundance matrix (SeqXXXXXX rows x sample columns).
  Used by 03_build_otu_table.R to reconstruct OTU-level abundances from VSEARCH results.

Map_file.csv
  Maps bioinformatics sample IDs (S001, S002...) back to original sample names.
  Used at the final step to restore original column names in the output table.
README

N_ASVS=$(grep -c "^>" "$DIR_INPUT/Centroid.fasta" || true)
plog "  ASVs in input FASTA: $N_ASVS"
plog_file "$DIR_INPUT/Centroid.fasta"          "All input ASVs with size= annotations"
plog_file "$DIR_INPUT/asv_abundance_table.txt" "ASV abundance matrix (SeqXXXXXX rows x sample columns)"
plog_file "$DIR_INPUT/Map_file.csv"            "Mapping of bioinformatics sample IDs to original sample names"

# ------------------------------------------------------------------------------
# Step 2: ITSx - Extract ITS region (ITS1/ITS2 only, if not skipped)
# ------------------------------------------------------------------------------

VSEARCH=$(which vsearch)

if [ "$RUN_ITSX" = true ]; then

    plog_section "Step 2: ITSx ${ITSX_REGION} Extraction"

    mkdir -p "$DIR_ITSX"

    ITSx -i "$DIR_INPUT/Centroid.fasta" \
         --complement F \
         -t F \
         --preserve T \
         --cpu ${THREADS} \
         -o "$DIR_ITSX/Centroid.ITSx"

    # Fix ITSx header format — ensure trailing semicolon after size= tag
    sed -i 's/;size=\([0-9]*\)/;size=\1;/g' "$DIR_ITSX/Centroid.ITSx.${ITSX_REGION}.fasta"

    # Remove zero-abundance sequences
    awk '/^>/{keep=1; if(/;size=0;/) keep=0} keep{print}' \
        "$DIR_ITSX/Centroid.ITSx.${ITSX_REGION}.fasta" > "$DIR_ITSX/Centroid.ITSx.${ITSX_REGION}.filtered.fasta"

    TMP_INPUT="$DIR_ITSX/Centroid.ITSx.${ITSX_REGION}.filtered.fasta"

    cat > "$DIR_ITSX/README.txt" << README
ITSx Extraction Files
=====================
Generated by: ITSx (${ITSX_REGION} region)

Centroid.ITSx.${ITSX_REGION}.fasta
  Raw ITSx output. ITSx identifies and extracts the ${ITSX_REGION} subregion,
  removing conserved flanking 5.8S sequence that inflates clustering similarity.

Centroid.ITSx.${ITSX_REGION}.filtered.fasta
  Final filtered FASTA used as input to VSEARCH clustering.
  Zero-abundance sequences removed; size= annotations restored.

Other Centroid.ITSx.* files
  Additional ITSx output files (other regions, summaries). Not used downstream.
README

    SEQS_AFTER_ITSX=$(grep -c "^>" "${TMP_INPUT}" || true)
    plog "  ASVs before ITSx:                    $N_ASVS"
    plog "  ASVs after ${ITSX_REGION} extraction:         $SEQS_AFTER_ITSX"
    plog "  ASVs removed (non-${ITSX_REGION}):            $((N_ASVS - SEQS_AFTER_ITSX))"

else

    plog_section "Step 2: ITSx — Skipped"
    if [[ "$PRIMER_SET" == ITS* ]]; then
        plog "  ITSx skipped by user request — using full ASV sequences"
    else
        plog "  ITSx not applicable for $PRIMER_SET — using full ASV sequences"
    fi
    TMP_INPUT="$DIR_INPUT/Centroid.fasta"
    SEQS_AFTER_ITSX=$N_ASVS

fi

# ------------------------------------------------------------------------------
# Step 3: VSEARCH Clustering and Cleaning
# ------------------------------------------------------------------------------

plog_section "Step 3: VSEARCH Clustering and Cleaning"

CLUSTER=0.97
TMP_SORTED=$(mktemp)
TMP_CLUSTERED=$(mktemp)
TMP_NOCHIM=$(mktemp)

# Sort by abundance
"${VSEARCH}" --fasta_width 0 \
    --sortbysize "${TMP_INPUT}" \
    --minseqlength 10 \
    --sizein --sizeout \
    --threads ${THREADS} \
    --output "${TMP_SORTED}" 2>/dev/null

# Cluster
"${VSEARCH}" --cluster_size "${TMP_SORTED}" \
    --id ${CLUSTER} \
    --sizein --sizeout \
    --minseqlength 10 \
    --qmask none \
    --threads ${THREADS} \
    --centroids "${TMP_CLUSTERED}" 2>&1 | tee -a "$LOG" > /dev/null

SEQS_AFTER_CLUSTER=$(grep -c "^>" "${TMP_CLUSTERED}" || true)
plog "  OTUs after clustering at ${CLUSTER} identity: $SEQS_AFTER_CLUSTER"
plog "  Sequences collapsed by clustering:   $((SEQS_AFTER_ITSX - SEQS_AFTER_CLUSTER))"

# Chimera checking
"${VSEARCH}" --uchime_denovo "${TMP_CLUSTERED}" \
    --sizein --sizeout \
    --qmask none \
    --threads ${THREADS} \
    --nonchimeras "${TMP_NOCHIM}" 2>&1 | tee -a "$LOG" > /dev/null

SEQS_AFTER_CHIMERA=$(grep -c "^>" "${TMP_NOCHIM}" || true)
CHIMERAS_REMOVED=$((SEQS_AFTER_CLUSTER - SEQS_AFTER_CHIMERA))
plog "  OTUs after chimera removal:          $SEQS_AFTER_CHIMERA"
plog "  Chimeras removed:                    $CHIMERAS_REMOVED ($(echo "scale=1; $CHIMERAS_REMOVED * 100 / $SEQS_AFTER_CLUSTER" | bc)%)"

# Sort & write final centroids
"${VSEARCH}" --fasta_width 0 \
    --minseqlength 10 \
    --sortbysize "${TMP_NOCHIM}" \
    --threads ${THREADS} \
    --sizein --sizeout \
    --output "$DIR_VSEARCH/${PROJ}.centroids" 2>/dev/null

# Map original ASVs to OTU centroids
"${VSEARCH}" --usearch_global "${TMP_INPUT}" \
    --db "$DIR_VSEARCH/${PROJ}.centroids" \
    --minseqlength 10 \
    --strand plus \
    --id ${CLUSTER} \
    --maxaccepts 0 \
    --qmask none \
    --threads ${THREADS} \
    --uc "$DIR_VSEARCH/${PROJ}.uc" 2>&1 | tee -a "$LOG" > /dev/null

rm -f "${TMP_SORTED}" "${TMP_CLUSTERED}" "${TMP_NOCHIM}"

# Build OTU abundance table
$RSCRIPT "$PIPELINE_DIR/back_ground_scripts/03_build_otu_table.R" \
    "$DIR_VSEARCH/${PROJ}.uc" \
    "$DIR_INPUT/asv_abundance_table.txt" \
    "$DIR_VSEARCH/${PROJ}.otutable" 2>&1 | tee -a "$LOG"

rm -f "$DIR_VSEARCH/${PROJ}.uc"

OTUS_IN_TABLE=$(awk 'NR>1' "$DIR_VSEARCH/${PROJ}.otutable" | wc -l | tr -d ' ')
plog "  OTUs in table before mumu:           $OTUS_IN_TABLE"

cat > "$DIR_VSEARCH/README.txt" << README
VSEARCH Clustering Files
========================
Generated by: VSEARCH + 03_build_otu_table.R

${PROJ}.centroids
  OTU centroid sequences after clustering at 97% identity and de novo chimera removal.
  Each sequence represents one OTU. Size annotations reflect total cluster abundance.

${PROJ}.otutable
  OTU abundance table (OTUs x samples) before mumu curation.
  Built by mapping all original ASVs back to OTU centroids and summing abundances.
README

plog_file "$DIR_VSEARCH/${PROJ}.centroids" "OTU centroid sequences after clustering and chimera removal"
plog_file "$DIR_VSEARCH/${PROJ}.otutable"  "OTU abundance table (OTUs x samples) before mumu curation"


# ------------------------------------------------------------------------------
# Step 4: Mumu Curation
# ------------------------------------------------------------------------------

plog_section "Step 4: Mumu Curation"

# Deactivate conda and load modules needed for mumu
conda deactivate
module purge
module load gcc/13.1.0-5z64cho
module load blast-plus/2.13.0-gcc-8.2.0-vo4mr4d

MUMU_BIN="$HOME/packages/mumu/mumu"

# Compile mumu if not already present
if [ ! -f "$MUMU_BIN" ]; then
    plog "  Compiling mumu..."
    mkdir -p "$HOME/packages"
    if [ ! -d "$HOME/packages/mumu" ]; then
        git clone https://github.com/frederic-mahe/mumu.git "$HOME/packages/mumu"
    fi
    (cd "$HOME/packages/mumu" && make && make check)
    chmod +x "$MUMU_BIN"
    plog "  mumu compiled successfully."
else
    plog "  mumu binary found, skipping compilation."
fi

[ -f "$MUMU_BIN" ] || { plog "ERROR: mumu binary not found after compilation."; exit 1; }

# Run mumu curation
cp "$DIR_VSEARCH/${PROJ}.centroids" "$DIR_MUMU/OTU_centroids"
sed -i 's/;.*//' "$DIR_MUMU/OTU_centroids"

makeblastdb -in "$DIR_MUMU/OTU_centroids" -parse_seqids -dbtype nucl > /dev/null 2>&1

# Set BLAST identity threshold for match list based on primer set
case "$PRIMER_SET" in
    16S-V4)
        MUMU_BLAST_ID=94
        ;;
    *)
        MUMU_BLAST_ID=84
        ;;
esac

plog "  mumu BLAST identity threshold: $MUMU_BLAST_ID%"

blastn -db "$DIR_MUMU/OTU_centroids" \
    -outfmt '6 qseqid sseqid pident' \
    -out "$DIR_MUMU/match_list.txt" \
    -qcov_hsp_perc 80 \
    -perc_identity ${MUMU_BLAST_ID} \
    -num_threads ${THREADS} \
    -query "$DIR_MUMU/OTU_centroids" 2>/dev/null

awk 'BEGIN{FS=OFS="\t"} {sub(/;.*/, "", $1); print}' "$DIR_VSEARCH/${PROJ}.otutable" > "$DIR_MUMU/mumu_table.txt"

# Set mumu minimum_ratio based on primer set
# Bacteria (16S) requires higher minimum_ratio to prevent merging of distinct taxa
case "$PRIMER_SET" in
    16S-V4)
        MUMU_MIN_RATIO=100
        ;;
    *)
        MUMU_MIN_RATIO=1
        ;;
esac

plog "  mumu minimum_ratio: $MUMU_MIN_RATIO"

"$MUMU_BIN" \
    --otu_table "$DIR_MUMU/mumu_table.txt" \
    --match_list "$DIR_MUMU/match_list.txt" \
    --log "$DIR_MUMU/mumu.log" \
    --new_otu_table "$DIR_MUMU/${PROJ}_mumu_curated.txt" \
    --minimum_ratio $MUMU_MIN_RATIO 2>&1 | tee -a "$LOG"

OTUS_AFTER_MUMU=$(awk 'NR>1' "$DIR_MUMU/${PROJ}_mumu_curated.txt" | wc -l | tr -d ' ')
MUMU_REMOVED=$((OTUS_IN_TABLE - OTUS_AFTER_MUMU))
plog "  OTUs after mumu curation: $OTUS_AFTER_MUMU"
plog "  OTUs removed by mumu: $MUMU_REMOVED"

grep -A 1 -Ff <(awk 'NR>1 {print $1}' "$DIR_MUMU/${PROJ}_mumu_curated.txt") "$DIR_MUMU/OTU_centroids" \
    | sed '/--/d' > "$DIR_MUMU/Centroid_mumu_curated.fas"

cat > "$DIR_MUMU/README.txt" << README
Mumu Curation Files
===================
Generated by: BLAST + mumu

OTU_centroids
  Centroid sequences with size= annotations stripped. Used as the BLAST database.

match_list.txt
  All-vs-all BLAST pairwise match list (query ID, subject ID, percent identity).
  Input to mumu for identifying potential parent-daughter OTU pairs.

mumu_table.txt
  OTU table reformatted for mumu input (OTU ID column stripped of size annotations).

mumu.log
  Detailed mumu curation log. Shows which OTUs were identified as daughters of
  parent OTUs and the abundance pattern evidence used for each decision.

${PROJ}_mumu_curated.txt
  OTU abundance table after mumu curation. Rare artefactual OTUs merged into parents.

Centroid_mumu_curated.fas
  Centroid sequences for the mumu-curated OTU set. Used as input for taxonomy assignment.
README

plog_file "$DIR_MUMU/${PROJ}_mumu_curated.txt"   "Mumu-curated OTU table with artefactual OTUs removed"
plog_file "$DIR_MUMU/mumu.log"                    "Detailed mumu curation log showing which OTUs were merged and why"
plog_file "$DIR_MUMU/match_list.txt"              "BLAST pairwise match list used as mumu input"
plog_file "$DIR_MUMU/Centroid_mumu_curated.fas"   "Centroid sequences for mumu-curated OTUs — used for taxonomy assignment"

# Reactivate conda for R steps
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate itsx_env
export PATH="$HOME/.conda/envs/itsx_env/bin:$PATH"

# ------------------------------------------------------------------------------
# Step 5: RDP Taxonomy Assignment via DADA2
# ------------------------------------------------------------------------------

plog_section "Step 5: RDP Taxonomy Assignment via DADA2"

START=$(date +%s)
$RSCRIPT "$PIPELINE_DIR/back_ground_scripts/05_assign_taxonomy_rdp.R" \
    "$DIR_MUMU/Centroid_mumu_curated.fas" \
    "$DB" \
    "$DIR_TAXONOMY/Taxonomy_rdp_${PRIMER_SET}.txt" \
    "${PRIMER_SET}" \
    "${DB_NAME}" 2>&1 | tee -a "$LOG"
END=$(date +%s)
plog "  Taxonomy assignment runtime: $(( (END - START) / 60 )) min $(( (END - START) % 60 )) sec"

CLASSIFIED=$(awk -F'\t' 'NR>1 && $3 != "NA"' "$DIR_TAXONOMY/Taxonomy_rdp_${PRIMER_SET}.txt" | wc -l | tr -d ' ')
TOTAL_OTUS=$(awk 'NR>1' "$DIR_TAXONOMY/Taxonomy_rdp_${PRIMER_SET}.txt" | wc -l | tr -d ' ')
plog "  OTUs with Kingdom assigned:          $CLASSIFIED of $TOTAL_OTUS"

cat > "$DIR_TAXONOMY/README.txt" << README
Taxonomy Assignment Files
=========================
Generated by: DADA2 assignTaxonomy (RDP Naive Bayesian classifier)
Database: ${DB_NAME}

Taxonomy_rdp_${PRIMER_SET}.txt
  Taxonomy assignments for each OTU centroid. Columns are standard taxonomic ranks
  (Kingdom through Species, or Domain through Species for PR2).

Taxonomy_rdp_${PRIMER_SET}_bootstraps.txt
  Bootstrap confidence values (0-100) for each taxonomic rank assignment.
  Higher values indicate more confident classification at that rank.

Taxonomy_rdp_${PRIMER_SET}_combined.txt
  Taxonomy and bootstrap values combined into a single table.
  Used as input to the final combination step (Step 6).
README

plog_file "$DIR_TAXONOMY/Taxonomy_rdp_${PRIMER_SET}.txt"            "RDP taxonomy assignments"
plog_file "$DIR_TAXONOMY/Taxonomy_rdp_${PRIMER_SET}_bootstraps.txt" "Per-rank bootstrap confidence values"
plog_file "$DIR_TAXONOMY/Taxonomy_rdp_${PRIMER_SET}_combined.txt"   "Taxonomy + bootstraps combined table"

# Deactivate conda environment
conda deactivate

# ------------------------------------------------------------------------------
# Step 6: Combine OTU Table with Taxonomy
# ------------------------------------------------------------------------------

plog_section "Step 6: Combine OTU Table with Taxonomy"

# For 18S-AMF datasets, output file will be named accordingly (unfiltered version created by Step 6)
FINAL_OTU_TAX="$PROJECT_DIR/${PROJ}_OTU_with_taxonomy_${PRIMER_SET}_${DB_NAME}.txt"

$RSCRIPT "$PIPELINE_DIR/back_ground_scripts/06_combine_otu_taxonomy.R" \
    "${PROJ}" \
    "${PRIMER_SET}" \
    "$DIR_MUMU/${PROJ}_mumu_curated.txt" \
    "$DIR_TAXONOMY/Taxonomy_rdp_${PRIMER_SET}_combined.txt" \
    "$DIR_INPUT/Map_file.csv" \
    "$FINAL_OTU_TAX" 2>&1 | tee -a "$LOG"

plog_file "$FINAL_OTU_TAX" \
    "Final output — OTU abundance table with original sample names and taxonomy"

# ------------------------------------------------------------------------------
# Pipeline Summary
# ------------------------------------------------------------------------------

plog_section "Pipeline Summary — Quality Control"
plog "  Primer set:                                      $PRIMER_SET"
plog "  SLURM job ID:                                    ${SLURM_JOB_ID:-N/A}"
plog "  Input ASVs:                                      $N_ASVS"
if [ "$RUN_ITSX" = true ]; then
plog "  After ITSx ${ITSX_REGION} extraction:                   $SEQS_AFTER_ITSX"
fi
plog "  After clustering (${CLUSTER} identity):                $SEQS_AFTER_CLUSTER"
plog "  Chimeras removed:                                $CHIMERAS_REMOVED ($(echo "scale=1; $CHIMERAS_REMOVED * 100 / $SEQS_AFTER_CLUSTER" | bc)%)"
plog "  OTUs after chimera removal:                      $SEQS_AFTER_CHIMERA"
plog "  OTUs removed by mumu:                            $MUMU_REMOVED ($(echo "scale=1; $MUMU_REMOVED * 100 / $OTUS_IN_TABLE" | bc)%)"
plog "  Final OTU count (post-mumu):                     $OTUS_AFTER_MUMU"
plog "  OTUs with taxonomy assigned:                     $CLASSIFIED of $TOTAL_OTUS"

plog_section "Output Files"
plog "  01_input/"
plog_file "$DIR_INPUT/Centroid.fasta"                               "All input ASVs with size annotations"
plog_file "$DIR_INPUT/asv_abundance_table.txt"                      "ASV abundance matrix used for OTU reconstruction"
plog_file "$DIR_INPUT/Map_file.csv"                                 "Sample ID mapping (S001... to original names)"
if [ "$RUN_ITSX" = true ]; then
plog "  02_itsx/"
plog_file "$DIR_ITSX/Centroid.ITSx.${ITSX_REGION}.filtered.fasta"  "Filtered ITSx output used for clustering"
fi
plog "  03_vsearch/"
plog_file "$DIR_VSEARCH/${PROJ}.centroids"                          "OTU centroids post-clustering and chimera removal"
plog_file "$DIR_VSEARCH/${PROJ}.otutable"                           "OTU abundance table before mumu"
plog "  04_mumu/"
plog_file "$DIR_MUMU/match_list.txt"                                "BLAST pairwise match list for mumu"
plog_file "$DIR_MUMU/mumu.log"                                      "Detailed mumu curation decisions"
plog_file "$DIR_MUMU/${PROJ}_mumu_curated.txt"                      "OTU table after mumu curation"
plog_file "$DIR_MUMU/Centroid_mumu_curated.fas"                     "Centroid sequences for mumu-curated OTUs"
plog "  05_taxonomy/"
plog_file "$DIR_TAXONOMY/Taxonomy_rdp_${PRIMER_SET}.txt"            "RDP taxonomy assignments"
plog_file "$DIR_TAXONOMY/Taxonomy_rdp_${PRIMER_SET}_bootstraps.txt" "Bootstrap confidence values per rank"
plog_file "$DIR_TAXONOMY/Taxonomy_rdp_${PRIMER_SET}_combined.txt"   "Taxonomy + bootstraps combined table"
plog "  (top level)"
plog_file "$PROJECT_DIR/pipeline_run.log"                           "This log file"

if [ "$PRIMER_SET" = "18S-AMF" ] && [ "$SKIP_AMF_FILTER" = false ]; then
    plog_file "$FINAL_OTU_TAX" \
        "FINAL OUTPUT (UNFILTERED) — OTU table with MaarjAM taxonomy (all OTUs)"
    plog_file "$FILTERED_OTU_TAX" \
        "FINAL OUTPUT (FILTERED) — OTU table with MaarjAM taxonomy (Mucoromycota only, SILVA-validated)"
    plog_file "$FILTER_SUMMARY" \
        "Filtering summary — details on which OTUs were retained/removed and why"
else
    plog_file "$FINAL_OTU_TAX" \
        "FINAL OUTPUT — OTU table with original sample names and taxonomy"
fi

plog ""
plog "Run completed: $(date)"

echo ""
echo "Pipeline complete. Log saved to: $LOG"