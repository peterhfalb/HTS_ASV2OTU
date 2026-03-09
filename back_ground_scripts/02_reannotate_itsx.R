# 02_reannotate_itsx.R
# Re-adds size= abundance annotations to ITSx output using original ASV abundances
# Usage: Rscript 02_reannotate_itsx.R <itsx_fasta> <asv_abundance_table> <output_fasta>
# Use personal R library (required on clusters)

options(repos = c(CRAN = "https://cloud.r-project.org"))  # set CRAN mirror

# ------------------------------------------------------------------------------
# Library paths — must be set before any library() calls
# ------------------------------------------------------------------------------

# Conda env library (highest priority — has dada2, tidyverse etc.)
conda_lib <- file.path(Sys.getenv("HOME"), ".conda/envs/itsx_env/lib/R/library")

# Personal R library (fallback for any manually installed packages)
personal_lib <- Sys.getenv("R_LIBS_USER")
if (personal_lib == "") {
  personal_lib <- file.path(path.expand("~"), "R",
                            paste0(R.version$platform, "-library"),
                            paste0(R.version$major, ".",
                                   substr(R.version$minor, 1, 1)))
}
dir.create(personal_lib, recursive = TRUE, showWarnings = FALSE)

# Set library paths — conda first, then personal, then system defaults
if (file.exists(conda_lib)) {
  .libPaths(c(conda_lib, personal_lib))
} else {
  .libPaths(c(personal_lib, .libPaths()))
}

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(purrr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
ITSX_FASTA  <- args[1]
ASV_ABUND   <- args[2]
OUTPUT      <- args[3]

# Load abundance table and compute total reads per ASV
abund <- read_tsv(ASV_ABUND, col_types = cols(.default = "d", SeqID = "c"))
totals <- abund %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  select(SeqID, Total)

# Read ITSx fasta — ITSx may modify headers slightly, so strip back to SeqXXXXXX
lines <- readLines(ITSX_FASTA)
out_lines <- character(length(lines))

for (i in seq_along(lines)) {
  if (startsWith(lines[i], ">")) {
    # Extract SeqXXXXXX ID — take first field before any space or pipe
    seq_id <- sub("^>([^ |]+).*", "\\1", lines[i])
    # Strip any existing size annotation
    seq_id <- sub(";size=.*", "", seq_id)
    # Look up total abundance
    total <- totals$Total[totals$SeqID == seq_id]
    if (length(total) == 0 || is.na(total)) {
      total <- 1  # fallback if not found
      warning("SeqID not found in abundance table: ", seq_id)
    }
    out_lines[i] <- paste0(">", seq_id, ";size=", round(total))
  } else {
    out_lines[i] <- lines[i]
  }
}

writeLines(out_lines, OUTPUT)
cat("Re-annotated", sum(startsWith(lines, ">")), "sequences with abundance tags\n")
cat("Output written to:", OUTPUT, "\n")