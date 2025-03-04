#!/usr/bin/env Rscript
# Install required R packages for RNA-Seq pipeline
# Author: Ajuni Sohota
# Usage: Rscript install_r_packages.R

# Function to install packages if not already installed
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    message(paste0("Installing ", package, "..."))
    
    # BiocManager for Bioconductor packages
    if (package %in% c("DESeq2", "edgeR", "limma", "clusterProfiler", 
                      "GO.db", "org.Hs.eg.db", "org.Mm.eg.db", 
                      "KEGG.db", "pathview", "tximport", "apeglm",
                      "GenomicFeatures", "rtracklayer")) {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cran.rstudio.com/")
      }
      BiocManager::install(package, update = FALSE, ask = FALSE)
    } else {
      # CRAN packages
      install.packages(package, repos = "https://cran.rstudio.com/")
    }
    
    # Check if installation was successful
    if (require(package, character.only = TRUE, quietly = TRUE)) {
      message(paste0(package, " installed successfully!"))
    } else {
      warning(paste0("Failed to install ", package))
    }
  } else {
    message(paste0(package, " is already installed."))
  }
}

# Main function
main <- function() {
  message("Installing required R packages for RNA-Seq pipeline...")
  
  # CRAN packages
  cran_packages <- c(
    "ggplot2",       # Plotting
    "pheatmap",      # Heatmaps
    "RColorBrewer",  # Color palettes
    "dplyr",         # Data manipulation
    "tidyr",         # Data tidying
    "readr",         # File reading
    "tibble",        # Modern data frames
    "gplots",        # Additional plotting
    "ggrepel",       # Text repelling for plots
    "reshape2",      # Data reshaping
    "yaml",          # YAML parsing
    "jsonlite",      # JSON manipulation
    "glue",          # String interpolation
    "optparse",      # Command line options
    "knitr",         # Report generation
    "rmarkdown",     # R Markdown support
    "plotly",        # Interactive plots
    "viridis",       # Color scales
    "DT",            # Interactive tables
    "heatmaply",     # Interactive heatmaps
    "stringr",       # String manipulation
    "scales"         # Scale functions for visualization
  )
  
  # Bioconductor packages
  bioc_packages <- c(
    "DESeq2",        # Differential expression
    "edgeR",         # Differential expression alternative
    "limma",         # Linear models for microarray/RNA-Seq data
    "clusterProfiler", # Functional enrichment
    "org.Hs.eg.db",  # Human genome annotation
    "org.Mm.eg.db",  # Mouse genome annotation
    "pathview",      # Pathway visualization
    "tximport",      # Transcript import for Salmon/Kallisto
    "apeglm",        # Shrinkage estimators for DESeq2
    "EnhancedVolcano", # Better volcano plots
    "ComplexHeatmap", # Advanced heatmaps
    "GenomicFeatures", # For working with genomic annotations
    "rtracklayer"    # Interface to genome browsers
  )
  
  # Install CRAN packages
  message("\nInstalling CRAN packages...")
  for (pkg in cran_packages) {
    install_if_missing(pkg)
  }
  
  # Install BiocManager if not already installed
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cran.rstudio.com/")
  }
  
  # Install Bioconductor packages
  message("\nInstalling Bioconductor packages...")
  for (pkg in bioc_packages) {
    install_if_missing(pkg)
  }
  
  message("\nPackage installation completed!")
  
  # Print version information
  message("\nR version information:")
  print(version)
  
  message("\nBioconductor version:")
  print(BiocManager::version())
  
  message("\nInstalled packages:")
  installed_packages <- installed.packages()[, c("Package", "Version")]
  all_packages <- c(cran_packages, bioc_packages)
  installed_required <- installed_packages[installed_packages[, "Package"] %in% all_packages, ]
  print(installed_required)
  
  # Check if any packages failed to install
  missing_packages <- all_packages[!all_packages %in% installed_packages[, "Package"]]
  if (length(missing_packages) > 0) {
    warning("The following packages failed to install: ", paste(missing_packages, collapse = ", "))
  } else {
    message("\nAll required packages were successfully installed!")
  }
}

# Run the main function
main()
