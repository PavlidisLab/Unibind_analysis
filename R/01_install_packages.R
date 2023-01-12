## Install packages used throughout analysis
## -----------------------------------------------------------------------------


packages <- c(
  "tidyverse",
  "BiocManager",
  "parallel",
  "pheatmap",
  "cowplot",
  "RColorBrewer",
  "preprocessCore",
  "GenomicRanges",
  "limma",
  "edgeR"
)


installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("biomaRt")
