## Export a tsv of basic Unibind metadata, adding count of regions per experiment
## and creating a version that collapses samples duplicated over different motifs.
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("R/Utils/functions.R")
source("R/00_config.R")


# Init metadata df from files for each TF. Assumes that path points to a main 
# dir which contains only subdirs named with each TF. Note that the permissive
# collection is a superset of the robust collection - the individual experiments
# common to each are identical.
# ------------------------------------------------------------------------------


meta_df <- function(path) {
  
  tfs <- list.files(path)
  
  stopifnot(length(tfs) > 0)
  
  files <- lapply(tfs, function(x) list.files(file.path(path, x)))
  names(files) <- tfs
  
  df <- data.frame(
    Symbol = unlist(lapply(tfs, function(x) rep(x, length(files[[x]])))),
    File = unlist(files))
  
  df$ID <- str_replace(df$File, "\\.MA.*", "")
  df$ID <- str_replace_all(df$ID, "\\.", "_")

  return(df)
}


meta_l <- list(
  Permissive_hg = meta_df(perm_path_hg),
  Robust_hg = meta_df(rob_path_hg),
  Permissive_mm = meta_df(perm_path_mm),
  Robust_mm = meta_df(rob_path_mm)
)


stopifnot(all(meta_l$Robust_hg$File %in% meta_l$Permissive_hg$File))
stopifnot(all(meta_l$Robust_mm$File %in% meta_l$Permissive_mm$File))


# Some TF/experiments have "duplicated" samples where the same experiment is
# scored using different motifs. Identify all such cases.
# ------------------------------------------------------------------------------


is_dup <- function(df) {
  df$Duplicate <- df$ID %in% df$ID[duplicated(df$ID)]
  return(df)
}


meta_l <- lapply(meta_l, is_dup)


# Load permissive files to get count of peaks for each experiment, then join 
# with robust (as the matched experiments are identical). 
# Note that this step is redundant with the following scoring steps, since they 
# also load the data. For time being wanted to separate meta construction and 
# scoring. (Also really should just use a wc -l > output.tsv command...)
# ------------------------------------------------------------------------------


count_regions <- function(dir, input_df, ncores = 1) {
  
  counts <- mclapply(1:nrow(input_df), function(x) {
    
    path <- paste0(dir, "/", input_df$Symbol[x], "/", input_df$File[x])
    
    tryCatch({
      n <- nrow(load_unibind(path))
      message(input_df$ID[x], " complete ", Sys.time())
      return(n)
    }, error = function(e) NULL)
    
  }, mc.cores = ncores)
  
  return(counts)
}



meta_l$Permissive_hg$N_peaks <-
  unlist(count_regions(perm_path_hg, meta_l$Permissive_hg))

meta_l$Permissive_mm$N_peaks <-
  unlist(count_regions(perm_path_mm, meta_l$Permissive_mm))


meta_l$Robust_hg <- left_join(meta_l$Robust_hg,
                              meta_l$Permissive_hg[, c("File", "N_peaks")],
                              by = "File")

meta_l$Robust_mm <- left_join(meta_l$Robust_mm,
                              meta_l$Permissive_mm[, c("File", "N_peaks")],
                              by = "File")


# Tallying TFs and seeing which are depleted in the robust set. In human, 59 TFs
# have data only in the permissive set. In mouse 50 TFs have data only in
# permissive set
# ------------------------------------------------------------------------------


n_tf <- lapply(meta_l, function(x) arrange(count(x, Symbol), desc(n)))


diff_hg <- left_join(n_tf$Permissive_hg, 
                     n_tf$Robust_hg, 
                     by = "Symbol",
                     suffix = c("Permissive", "Robust"))

diff_mm <- left_join(n_tf$Permissive_mm, 
                     n_tf$Robust_mm, 
                     by = "Symbol",
                     suffix = c("Permissive", "Robust"))

perm_only_hg <- filter(diff_hg, is.na(nRobust))

perm_only_mm <- filter(diff_mm, is.na(nRobust))


# Save out
# ------------------------------------------------------------------------------

saveRDS(meta_l, meta_outfile)
