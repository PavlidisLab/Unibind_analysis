## Creating a gene by experiment matrix of binding scores for mouse experiments
## from the Unibind database. Because all robust experiments are contained
## within the permissive set, score all permissive experiments.
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("R/00_config.R")
source("R/Utils/functions.R")

# Load metadata, protein coding genes, and ENCODE blacklisted regions
pc <- read.delim(ref_mm, stringsAsFactors = FALSE)
bl <- read.delim(bl_path_mm, stringsAsFactors = FALSE)
meta_l <- readRDS(meta_outfile)

# Convert range tables to GR objects for scoring
bl_gr <- bl_to_gr(bl)
pc_gr <- pc_to_gr(pc)

# Remove experiments under peak cutoff
meta <- filter(meta_l$Permissive_mm, N_peaks >= min_peaks)


# Gene scoring: Each experiment from meta is loaded, converted to a GRanges
# object, and a binding score for each gene is calculated. These vectors of
# binding scores are then bound in a matrix.
# NOTE: SLOW! Took ~ 30 hours on frink with 8 cores
# -----------------------------------------------------------------------------


if (!file.exists(bmat_path_mm)) {
  
  bscore_l <- load_and_score(
    dir = perm_path_mm,
    input_df = meta,
    bl_gr = bl_gr,
    pc_gr = pc_gr,
    ncores = cores
  )
  
  names(bscore_l) <- meta$File
  
  bscore_mat <- do.call(cbind, bscore_l)
  
  saveRDS(bscore_mat, file = bmat_path_mm)
  
}
