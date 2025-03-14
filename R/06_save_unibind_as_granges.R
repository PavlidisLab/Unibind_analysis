## Saving Unibind experiments as GenomicRanges list object
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
source("R/utils/functions.R")
source("R/00_config.R")

# ENCODE blacklisted regions to be filtered
bl_hg <- bl_to_gr(read.delim(bl_path_hg, stringsAsFactors = FALSE))
bl_mm <- bl_to_gr(read.delim(bl_path_mm, stringsAsFactors = FALSE))

# Meta data of experiments to be loaded
meta_l <- readRDS(meta_path)


# Load each experiment of the input/meta df into a list of GR objects
# ------------------------------------------------------------------------------


load_gr <- function(dir, input_df, bl_gr) {
  
  gr_l <- lapply(1:nrow(input_df), function(x) {
    
    message(paste(input_df$File[x], Sys.time()))
    path <- file.path(dir, input_df$Symbol[x], input_df$File[x])
    
    tryCatch({
      
      gr <- load_unibind(path) %>% unibind_to_gr(bl_gr)

    }, error = function(e) NULL)
    
  })
  
  names(gr_l) <- input_df$File
  return(gr_l)
}


# Load human permissive and create subset list of robust collection

gr_perm_hg <- load_gr(dir = perm_path_hg, 
                      input_df = meta_l$Permissive_hg, 
                      bl_gr = bl_hg)

gr_rob_hg <- gr_perm_hg[intersect(meta_l$Robust_hg$File, meta_l$Permissive_hg$File)]

# Ditto for mouse

gr_perm_mm <- load_gr(dir = perm_path_mm, 
                      input_df = meta_l$Permissive_mm, 
                      bl_gr = bl_mm)

gr_rob_mm <- gr_perm_mm[intersect(meta_l$Robust_mm$File, meta_l$Permissive_mm$File)]


# Save out
# ------------------------------------------------------------------------------


saveRDS(gr_perm_hg, file = gr_perm_path_hg)
saveRDS(gr_rob_hg, file = gr_rob_path_hg)
saveRDS(gr_perm_mm, file = gr_perm_path_mm)
saveRDS(gr_rob_mm, file = gr_rob_path_mm)
