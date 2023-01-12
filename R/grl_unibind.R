## Saving Unibind experiments as GenomicRanges list object
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
source("~/regnetR/R/utils/bscore_unibind_functions.R")

# Note only need permissive path as it is a superset of robust experiments
# and the corresponding experiments are identical
path_perm_hg <- "/home/amorin/Data/Peak_files/Unibind/Hg38/Permissive/damo_hg38_TFBS_per_TF/"
path_perm_mm <- "/home/amorin/Data/Peak_files/Unibind/Mm10/Permissive/damo_mm10_TFBS_per_TF/"

# Meta pointing to paths of robust and permissive collections
meta_perm_hg <- read.delim("~/Data/Metadata/Chipseq/unibind_human.tsv", stringsAsFactors = FALSE)
meta_rob_hg <- read.delim("~/Data/Metadata/Chipseq/unibind_human.tsv", stringsAsFactors = FALSE)
meta_perm_mm <- read.delim("~/Data/Metadata/Chipseq/unibind_mouse.tsv", stringsAsFactors = FALSE)
meta_rob_mm <- read.delim("~/Data/Metadata/Chipseq/unibind_mouse.tsv", stringsAsFactors = FALSE)

# Blacklisted regions to be removed
bl_hg <- read.delim("~/Data/Chromosome_info/blacklist_hg38.tsv", stringsAsFactors = FALSE)
bl_mm <- read.delim("~/Data/Chromosome_info/blacklist_mm10.tsv", stringsAsFactors = FALSE)
bl_gr_hg <- bl_to_gr(bl_hg)
bl_gr_mm <- bl_to_gr(bl_mm)


# Load each experiment of the input/meta df into a list of GR objects


load_gr <- function(dir, input_df, bl_gr) {
  
  gr_l <- lapply(1:nrow(input_df), function(x) {
    
    path <- paste0(dir, "/", input_df$Symbol[x], "/", input_df$File[x])
    
    tryCatch({
      
      peak_gr <- load_unibind(path) %>% 
        unibind_to_gr(bl_gr)
      
      message(input_df$File[x], " complete ", Sys.time())
      
      return(peak_gr)
      
    }, error = function(e) NULL)
    
  })
  
  names(gr_l) <- input_df$File
  return(gr_l)
}


# Load human permissive and create subset list of robust collection
gr_perm_hg <- load_gr(path_perm_hg, meta_perm_hg, bl_gr_hg)
gr_rob_hg <- gr_perm_hg[intersect(meta_rob_hg$File, meta_perm_hg$File)]

# Ditto for mouse
gr_perm_mm <- load_gr(path_perm_mm, meta_perm_mm, bl_gr_mm)
gr_rob_mm <- gr_perm_mm[intersect(meta_rob_mm$File, meta_perm_mm$File)]


saveRDS(gr_perm_hg,
        file = "~/scratch/R_objects/unibind_grlist_perm_human.RDS")

saveRDS(gr_rob_hg,
        file = "~/scratch/R_objects/unibind_grlist_rob_human.RDS")

saveRDS(gr_perm_mm,
        file = "~/scratch/R_objects/unibind_grlist_perm_mouse.RDS")

saveRDS(gr_rob_mm,
        file = "~/scratch/R_objects/unibind_grlist_rob_mouse.RDS")
