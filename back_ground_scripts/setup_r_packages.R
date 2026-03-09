# setup_r_packages.R
# Run this once to install all required R packages
# Handles both local and cluster environments
# Usage: Rscript setup_r_packages.R

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ------------------------------------------------------------------------------
# Set up personal R library
# ------------------------------------------------------------------------------

personal_lib <- Sys.getenv("R_LIBS_USER")
if (personal_lib == "") {
  personal_lib <- file.path(path.expand("~"), "R",
                            paste0(R.version$platform, "-library"),
                            paste0(R.version$major, ".",
                                   substr(R.version$minor, 1, 1)))
}

cat("Personal R library path:", personal_lib, "\n")
dir.create(personal_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(personal_lib, .libPaths()))
cat("Library paths:\n")
print(.libPaths())

# ------------------------------------------------------------------------------
# Helper function — install with error handling
# ------------------------------------------------------------------------------

install_if_missing <- function(pkg, bioc = FALSE, skip_on_fail = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    tryCatch({
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", lib = personal_lib)
        }
        BiocManager::install(pkg, lib = personal_lib, ask = FALSE,
                             update = FALSE, dependencies = FALSE)
      } else {
        install.packages(pkg, lib = personal_lib)
      }
    }, error = function(e) {
      cat("ERROR during installation of", pkg, ":", conditionMessage(e), "\n")
    })
    
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (skip_on_fail) {
        cat("WARNING: Could not install", pkg, "— skipping\n")
        return(invisible(FALSE))
      } else {
        stop("Failed to install package: ", pkg)
      }
    } else {
      cat("Successfully installed:", pkg, "\n")
      return(invisible(TRUE))
    }
  } else {
    cat("Already installed:", pkg, "\n")
    return(invisible(TRUE))
  }
}

# ------------------------------------------------------------------------------
# CRAN packages
# ------------------------------------------------------------------------------

cat("\n--- Installing CRAN packages ---\n")

# Install tidyverse core packages individually
# (avoids ragg dependency which requires system graphics libs not on all clusters)
tidyverse_core <- c("dplyr", "tidyr", "readr", "tibble", "stringr", "purrr")
for (pkg in tidyverse_core) {
  install_if_missing(pkg)
}

# ragg is optional — skip if system graphics libs not available
install_if_missing("ragg", skip_on_fail = TRUE)

# ------------------------------------------------------------------------------
# Bioconductor packages — installed in explicit dependency order
# This is necessary on clusters with old conda/mamba that can't solve
# complex Bioconductor dependency chains automatically
# ------------------------------------------------------------------------------

cat("\n--- Installing Bioconductor packages (explicit dependency order) ---\n")

install_if_missing("BiocManager")
library(BiocManager)

# Tier 1 — base Bioconductor infrastructure
bioc_tier1 <- c("BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb",
                "XVector", "Biostrings")
for (pkg in bioc_tier1) {
  install_if_missing(pkg, bioc = TRUE)
}

# Tier 2 — depends on tier 1
bioc_tier2 <- c("MatrixGenerics", "DelayedArray", "GenomicRanges")
for (pkg in bioc_tier2) {
  install_if_missing(pkg, bioc = TRUE)
}

# Tier 3 — depends on tier 2
bioc_tier3 <- c("SummarizedExperiment", "Biobase", "Rsamtools")
for (pkg in bioc_tier3) {
  install_if_missing(pkg, bioc = TRUE)
}

# Tier 4 — depends on tier 3
bioc_tier4 <- c("GenomicAlignments", "BiocParallel")
for (pkg in bioc_tier4) {
  install_if_missing(pkg, bioc = TRUE)
}

# Tier 5 — ShortRead depends on tier 4
# latticeExtra is a CRAN package required by ShortRead
install_if_missing("latticeExtra", skip_on_fail = TRUE)
install_if_missing("ShortRead", bioc = TRUE)

# Final — dada2
install_if_missing("dada2", bioc = TRUE)

# ------------------------------------------------------------------------------
# Final verification
# ------------------------------------------------------------------------------

cat("\n--- Package verification ---\n")

required <- list(
  cran = c("dplyr", "tidyr", "readr", "tibble", "stringr", "purrr"),
  bioc = c("Biostrings", "dada2")
)

all_ok <- TRUE
for (pkg in c(required$cran, required$bioc)) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %-20s OK  (%s)\n", pkg, packageVersion(pkg)))
  } else {
    cat(sprintf("  %-20s MISSING\n", pkg))
    all_ok <- FALSE
  }
}

if (!all_ok) {
  stop("Some required packages are missing. Check errors above.")
} else {
  cat("\nAll required packages installed successfully.\n")
  cat("You can now run the pipeline.\n")
}