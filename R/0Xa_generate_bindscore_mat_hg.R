## Creating a gene by experiment matrix of binding scores for human experiments
## from the Unibind database.
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("R/00_config.R")

path_hg <- "/home/amorin/Data/Peak_files/Unibind/Hg38/Permissive/damo_hg38_TFBS_per_TF/"
pc_hg <- read.delim("~/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)
bl_hg <- read.delim("~/Data/Chromosome_info/blacklist_hg38.tsv", stringsAsFactors = FALSE)
meta_hg <- read.delim("~/Data/Metadata/Chipseq/unibind_human.tsv", stringsAsFactors = FALSE)
meta_rob_hg <- read.delim("~/Data/Metadata/Chipseq/unibind_human.tsv", stringsAsFactors = FALSE)

# Remove duplicated pseudoautosomal genes (keep X copy)
dupl <- pc_hg$Symbol[duplicated(pc_hg$Symbol)]
pc_hg <- filter(pc_hg, !(Symbol %in% dupl & Chromosome == "Y"))

# Range tables -> GR objects
bl_gr_hg <- bl_to_gr(bl_hg)
pc_gr_hg <- pc_to_gr(pc_hg)


# Remove experiments under peak cutoff
meta_hg <- filter(meta_hg, N >= 100)


# Scoring all is slow - took ~ 30 hours on frink with 8 cores

if(!file.exists("~/scratch/R_objects/unibind_all_scores_hg.RDS")) {
  bscore_hg <- load_and_score(path_hg, meta_hg, bl_gr_hg, pc_gr_hg)
  names(bscore_hg) <- meta_hg$File
  saveRDS(bscore_hg, file = "~/scratch/R_objects/unibind_all_scores_hg.RDS")
} else {
  bscore_hg <- readRDS("~/scratch/R_objects/unibind_all_scores_hg.RDS")
}


bmat_hg <- do.call(cbind, bscore_hg)