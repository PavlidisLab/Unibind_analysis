## Saving Unibind experiments as GenomicRanges list object
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(parallel)
source("~/regnetR/R/utils/bscore_unibind_functions.R")

#
gr_hg <- readRDS("~/scratch/R_objects/unibind_grlist_human.RDS")
gr_mm <- readRDS("~/scratch/R_objects/unibind_grlist_mouse.RDS")


#
pc_hg <- read.delim("~/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)
dupl <- pc_hg$Symbol[duplicated(pc_hg$Symbol)]
pc_hg <- filter(pc_hg, !(Symbol %in% dupl & Chromosome == "Y"))
pc_gr_hg <- pc_to_gr(pc_hg)

pc_mm <- read.delim("~/Data/Metadata/refseq_select_mm10.tsv", stringsAsFactors = FALSE)
pc_gr_mm <- pc_to_gr(pc_mm)


# 
# ------------------------------------------------------------------------------



nearest_dist_mat <- function(pc_gr, peak_gr_l, cores = 8) {
  
  dist_l <- mclapply(peak_gr_l, function(x) {
    dist_vec <- rep(NA, length(pc_gr))
    dist_gr <- distanceToNearest(pc_gr, x)
    dist_vec[dist_gr@from] <- dist_gr@elementMetadata$distance
    return(dist_vec)
  }, mc.cores = cores)
  
  dist_mat <- do.call(rbind, dist_l)
  colnames(dist_mat) <- pc_gr$Symbol
  return(dist_mat)
}


# examples of genes with high and low aggregate binding scores

hi_hg <- c("IRF2BP2")
low_hg <- c("OR4F29")

mb_hg <- pc_gr_hg[pc_gr_hg$Symbol %in% mb_hg]

dist_hg <- nearest_dist_mat(
  pc_gr = pc_gr_hg[pc_gr_hg$Symbol %in% c(hi_hg, low_hg)], 
  peak_gr_l = gr_hg$Robust
)

summary(dist_hg)
view(data.frame(dist_hg))


plot(density(log10(dist_hg[, "IRF2BP2"]+1), na.rm = TRUE), ylim = c(0, 2.5))
lines(density(log10(dist_hg[, "OR4F29"]+1), na.rm = TRUE), col = "red")

sum(dist_hg[, "IRF2BP2"] < 10e3, na.rm = TRUE)/nrow(tt2)
sum(dist_hg[, "OR4F29"] < 10e3, na.rm = TRUE)/nrow(tt2)


plot(sort(log10(dist_hg[, "IRF2BP2"]+1)))
abline(h = log10(10e3+1))
     