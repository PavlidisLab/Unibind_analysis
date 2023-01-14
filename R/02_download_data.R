## Download Unibind data, ...
## TODO: port over pcoding/ccre/etc table calls from TR agg OR make package
## TODO: robust/permissive def
## TODO: look at compressed https://unibind.uio.no/static/data/20220914/bulk_Robust/Homo_sapiens/hg38_compressed_TFBSs.bed.gz
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")

options(timeout = 1000)  # default 60 was insufficient


# Unibind data: Robust versus Permissive collection of peaks are each saved
# in a single directory for all experiments.
# Permissive:
# Robust: 
# ------------------------------------------------------------------------------


rob_link_hg <- "https://unibind.uio.no/static/data/20220914/bulk_Robust/Homo_sapiens/damo_hg38_TFBS_per_TF.tar.gz"
rob_link_mm <- "https://unibind.uio.no/static/data/20220914/bulk_Robust/Mus_musculus/damo_mm10_TFBS_per_TF.tar.gz"
perm_link_hg <- "https://unibind.uio.no/static/data/20220914/bulk_Permissive/Homo_sapiens/damo_hg38_TFBS_per_TF.tar.gz"
perm_link_mm <- "https://unibind.uio.no/static/data/20220914/bulk_Permissive/Mus_musculus/damo_mm10_TFBS_per_TF.tar.gz"




input_df <- data.frame(
  link = c(rob_link_mm, perm_link_mm, rob_link_hg, perm_link_hg),
  path = c(rob_path_mm, perm_path_mm, rob_path_hg, perm_path_hg),
  stringsAsFactors = FALSE
)


input_df$path <- paste0(input_df$path, ".tar.gz")


# NOTE: The untar creates a dir with another dir one level down with identical 
# name, which then contains the individual TF dirs. I ended up just manually 
# moving these dirs a level up. 


for (i in 1:nrow(input_df)) {
  
  if (!file.exists(input_df$path[i])) {
    
    download.file(input_df$link[i], input_df$path[i])
    
    untar(tarfile = input_df$path[i], 
          exdir = str_replace(input_df$path[i], ".tar.gz", ""))
  }
}



# UK Biobank depletion ranks

dr_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9329122/bin/41586_2022_4965_MOESM3_ESM.gz"

if (!file.exists(dr_path)) {
  download.file(dr_url, dr_path)  
}
