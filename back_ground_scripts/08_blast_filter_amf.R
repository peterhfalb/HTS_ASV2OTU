# AMF Dataset Filtering: BLAST-based Quality Control
#
# Purpose: For 18S-AMF datasets, this script performs BLAST-based filtering
# against a curated 18S fungal sequence database to validate MaarjAM assignments.
# Only sequences with a top BLAST hit identified as Mucoromycota are retained.
#
# This approach:
# - Uses MaarjAM for detailed AMF identification (species-level resolution)
# - Uses BLAST against 18S fungal reference for Mucoromycota validation
# - Keeps MaarjAM taxonomy in final output (BLAST used only for filtering)
#
# Output:
#   - Two OTU+taxonomy files (unfiltered MaarjAM, filtered Mucoromycota)
#   - BLAST summary showing all hits and filtering decisions
#
# Usage: Rscript 08_blast_filter_amf.R <project_name> <maarjam_combined_tax> \
#          <fasta_path> <otu_abundance_table> <map_file> <blast_output> \
#          <output_unfiltered> <output_filtered> <summary_report>

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ------------------------------------------------------------------------------
# Library paths — must be set before any library() calls
# ------------------------------------------------------------------------------

conda_lib <- file.path(Sys.getenv("HOME"), ".conda/envs/itsx_env/lib/R/library")
personal_lib <- Sys.getenv("R_LIBS_USER")
if (personal_lib == "") {
  personal_lib <- file.path(path.expand("~"), "R",
                            paste0(R.version$platform, "-library"),
                            paste0(R.version$major, ".",
                                   substr(R.version$minor, 1, 1)))
}
dir.create(personal_lib, recursive = TRUE, showWarnings = FALSE)

if (file.exists(conda_lib)) {
  .libPaths(c(conda_lib, personal_lib))
} else {
  .libPaths(c(personal_lib, .libPaths()))
}

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: Rscript 08_blast_filter_amf.R <proj_name> <maarjam_combined_tax> ",
       "<fasta_path> <otu_abundance_table> <map_file> <blast_output> ",
       "<output_unfiltered> <output_filtered> <summary_report>")
}

PROJ             <- args[1]
MAARJAM_TAX_PATH <- args[2]
FASTA_PATH       <- args[3]
OTU_TABLE_PATH   <- args[4]
MAP_FILE_PATH    <- args[5]
BLAST_OUTPUT     <- args[6]
OUT_UNFILTERED   <- args[7]
OUT_FILTERED     <- args[8]
SUMMARY_PATH     <- args[9]

cat("Project name:             ", PROJ, "\n")
cat("MaarjAM taxonomy:         ", MAARJAM_TAX_PATH, "\n")
cat("Input FASTA:              ", FASTA_PATH, "\n")
cat("OTU abundance table:      ", OTU_TABLE_PATH, "\n")
cat("Map file:                 ", MAP_FILE_PATH, "\n")
cat("BLAST output:             ", BLAST_OUTPUT, "\n")
cat("Output (unfiltered):      ", OUT_UNFILTERED, "\n")
cat("Output (filtered):        ", OUT_FILTERED, "\n")
cat("Summary report:           ", SUMMARY_PATH, "\n\n")

# Derive taxonomy folder from summary_path
TAXONOMY_DIR <- dirname(SUMMARY_PATH)

# ------------------------------------------------------------------------------
# Read input files
# ------------------------------------------------------------------------------

cat("Reading input files...\n")

# Read MaarjAM taxonomy assignments (combined with bootstraps)
maarjam_tax <- read_delim(MAARJAM_TAX_PATH, delim = "\t", col_types = cols(.default = "c"))

# Read OTU abundance table
otu_table <- read_delim(OTU_TABLE_PATH, delim = "\t", col_types = cols(.default = "c"))

# Read map file (sample mapping)
map_file <- read_csv(MAP_FILE_PATH, col_types = cols(.default = "c"))

# Read BLAST results (format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle)
cat("Reading BLAST results...\n")
blast_results <- read_delim(BLAST_OUTPUT,
                            delim = "\t",
                            col_names = c("OTUId", "subject_id", "pident", "length", "mismatch",
                                         "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                                         "bitscore", "stitle"),
                            col_types = cols(.default = "c"),
                            skip = 0)

cat("BLAST results read: ", nrow(blast_results), " total hits\n")

# Get top BLAST hit per OTU (first hit is top hit from BLAST output)
top_hits <- blast_results %>%
  group_by(OTUId) %>%
  slice(1) %>%
  ungroup()

cat("Top BLAST hits: ", nrow(top_hits), " OTUs with hits\n")

# Extract organism name and check for Mucoromycota from subject title
# GenBank format is typically: "description [organism]"
# Example: "Rhizophagus irregularis 18S ribosomal RNA gene [Rhizophagus irregularis]"
top_hits <- top_hits %>%
  mutate(
    # Extract organism name from brackets at end
    organism = str_extract(stitle, "\\[([^\\]]+)\\]$", group = 1),
    # Check if description contains Mucoromycota, Glomeromycota, or other AMF clues
    is_mucoromycota = case_when(
      str_detect(stitle, regex("Mucoromycota|Mucoromycotina|mucoromycota", ignore_case = TRUE)) ~ TRUE,
      str_detect(organism, regex("Mucoromycota|Mucoromycotina|mucoromycota", ignore_case = TRUE)) ~ TRUE,
      # Also check organism names known to be Mucoromycota AMF
      str_detect(organism, regex("Rhizophagus|Funneliformis|Septoglomus|Claroideoglomus|Diversispora|Paraglomus|Scutellospora|Gigaspora|Racocetra", ignore_case = TRUE)) ~ TRUE,
      TRUE ~ FALSE
    ),
    # Check for uncultured/environmental
    is_uncultured = str_detect(stitle, regex("uncultured|environmental|metagenome|clone|sample", ignore_case = TRUE)),
    # Determine filter decision
    filter_decision = case_when(
      is.na(organism) | organism == "" ~ "REMOVED",
      is_uncultured ~ "REMOVED",
      is_mucoromycota ~ "KEPT",
      TRUE ~ "REMOVED"
    ),
    filter_reason = case_when(
      is.na(organism) | organism == "" ~ "No_organism_identified",
      is_uncultured ~ "Uncultured_or_environmental",
      is_mucoromycota ~ "Mucoromycota_confirmed",
      TRUE ~ paste0("Non_Mucoromycota: ", organism)
    )
  )

# OTUs with no BLAST hits
otus_no_hit <- setdiff(unique(otu_table$OTUId), top_hits$OTUId)
if (length(otus_no_hit) > 0) {
  cat("OTUs with no BLAST hits: ", length(otus_no_hit), "\n")
  # Add these as REMOVED
  no_hit_df <- data.frame(
    OTUId = otus_no_hit,
    subject_id = NA,
    pident = NA,
    length = NA,
    mismatch = NA,
    gapopen = NA,
    qstart = NA,
    qend = NA,
    sstart = NA,
    send = NA,
    evalue = NA,
    bitscore = NA,
    stitle = NA,
    organism = NA,
    is_mucoromycota = FALSE,
    is_uncultured = FALSE,
    filter_decision = "REMOVED",
    filter_reason = "No_BLAST_hit",
    stringsAsFactors = FALSE
  )
  top_hits <- bind_rows(top_hits, no_hit_df)
}

# Summary statistics
kept_count <- sum(top_hits$filter_decision == "KEPT")
removed_count <- sum(top_hits$filter_decision == "REMOVED")

cat("\nFiltering Summary:\n")
cat("  Total OTUs processed:        ", nrow(top_hits), "\n")
cat("  OTUs KEPT (Mucoromycota):    ", kept_count, "\n")
cat("  OTUs REMOVED:                ", removed_count, "\n")

# Breakdown of removals
if (removed_count > 0) {
  cat("\nRemoval Reason Breakdown:\n")
  removal_breakdown <- top_hits %>%
    filter(filter_decision == "REMOVED") %>%
    count(filter_reason) %>%
    arrange(desc(n))
  for (i in seq_len(nrow(removal_breakdown))) {
    cat("  ", removal_breakdown$filter_reason[i], ": ", removal_breakdown$n[i], "\n")
  }
}

# Summary of KEPT sequences
if (kept_count > 0) {
  cat("\nKEPT Sequences - Top organisms:\n")
  kept_organisms <- top_hits %>%
    filter(filter_decision == "KEPT") %>%
    count(organism) %>%
    arrange(desc(n)) %>%
    head(10)
  for (i in seq_len(nrow(kept_organisms))) {
    cat("  ", kept_organisms$organism[i], ": ", kept_organisms$n[i], "\n")
  }
}

# Write BLAST summary to TSV
cat("\nWriting BLAST filtering summary...\n")
write_delim(top_hits, SUMMARY_PATH, delim = "\t")

# ------------------------------------------------------------------------------
# Process OTU tables: Create filtered and unfiltered versions
# ------------------------------------------------------------------------------

cat("\nProcessing OTU abundance tables...\n")

# Get OTU IDs to keep
kept_otus <- top_hits$OTUId[top_hits$filter_decision == "KEPT"]

# Create unfiltered version (all OTUs with MaarjAM taxonomy)
otu_unfiltered <- otu_table

# Create filtered version (only kept OTUs)
otu_filtered <- otu_table %>%
  filter(OTUId %in% kept_otus)

cat("OTU counts:\n")
cat("  Unfiltered (all OTUs):        ", nrow(otu_unfiltered), "\n")
cat("  Filtered (Mucoromycota only): ", nrow(otu_filtered), "\n")
cat("  OTUs removed by filter:       ", nrow(otu_unfiltered) - nrow(otu_filtered), "\n")

# Get corresponding taxonomy rows
maarjam_unfiltered <- maarjam_tax
maarjam_filtered <- maarjam_tax %>%
  filter(OTUId %in% kept_otus)

# Combine OTU tables with MaarjAM taxonomy for both versions
# Match sample columns from original OTU table
sample_cols <- setdiff(colnames(otu_table), "OTUId")

final_unfiltered <- otu_unfiltered %>%
  left_join(maarjam_unfiltered, by = "OTUId") %>%
  select(OTUId, all_of(sample_cols), everything())

final_filtered <- otu_filtered %>%
  left_join(maarjam_filtered, by = "OTUId") %>%
  select(OTUId, all_of(sample_cols), everything())

# Write final output files
cat("\nWriting final OTU+taxonomy tables...\n")
write.table(final_unfiltered,
            OUT_UNFILTERED,
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(final_filtered,
            OUT_FILTERED,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nFiltering complete!\n")
cat("Outputs written to:\n")
cat("  Unfiltered OTU table: ", OUT_UNFILTERED, "\n")
cat("  Filtered OTU table:   ", OUT_FILTERED, "\n")
cat("  BLAST summary:        ", SUMMARY_PATH, "\n")
