## Export a tsv of basic unibind metadata, adding count of regions per experiment
## and creating a version that collapses samples duplicated over different motifs
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("/home/amorin/regnetR/R/utils/bscore_unibind_functions.R")

# Note that permissive collection (lax motif filterin) contains all robust
path_hg <- "/home/amorin/Data/Peak_files/Unibind/Hg38/Permissive/damo_hg38_TFBS_per_TF/"
path_rob_hg <- "/home/amorin/Data/Peak_files/Unibind/Hg38/Robust/damo_hg38_TFBS_per_TF/"
path_mm <- "/home/amorin/Data/Peak_files/Unibind/Mm10/Permissive/damo_mm10_TFBS_per_TF/"
path_rob_mm <- "/home/amorin/Data/Peak_files/Unibind/Mm10/Robust/damo_mm10_TFBS_per_TF/"


# Init metadata df from files for each TF. Assumes main dir only contains only
# subdirs named with TF.


meta_df <- function(path) {
  
  tfs <- list.files(path)
  
  files <- lapply(tfs, function(x) {
    list.files(paste0(path, x))
  })
  names(files) <- tfs
  
  df <- data.frame(
    Symbol = unlist(lapply(tfs, function(x) rep(x, length(files[[x]])))),
    File = unlist(files))
  
  df$ID <- str_replace(df$File, "\\.MA.*", "")
  
  return(df)
}


df_hg <- meta_df(path_hg)
df_rob_hg <- meta_df(path_rob_hg)
df_mm <- meta_df(path_rob_mm)
df_rob_mm <- meta_df(path_mm)


stopifnot(all(df_rob_hg$File) %in% df_hg$File)
stopifnot(all(df_rob_mm$File) %in% df_mm$File)


# Load permissive files to get npeak, then join with robust. 
# NOTE: This step is redundant with the following scoring steps, since they also
# load the data. For time being wanted to separate meta construction and scoring 
# (Also really should just use a wc -l > output.tsv command...)


count_regions <- function(dir, input_df) {
  
  counts <- lapply(1:nrow(input_df), function(x) {
    
    path <- paste0(dir, "/", input_df$Symbol[x], "/", input_df$File[x])
    
    tryCatch({
      n <- nrow(load_unibind(path))
      message(input_df$ID[x], " complete ", Sys.time())
      return(n)
    }, error = function(e) NULL)
    
  })
  
}


df_hg$N <- unlist(count_regions(path_hg, df_hg))
df_mm$N <- unlist(count_regions(path_mm, df_mm))

df_rob_hg <- left_join(df_rob_hg, df_hg[, c("File", "N")], by = "File")
df_rob_mm <- left_join(df_rob_mm, df_mm[, c("File", "N")], by = "File")


write.table(df_hg,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            file = "~/Data/Metadata/Chipseq/unibind_human.tsv")


write.table(df_rob_hg,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            file = "~/Data/Metadata/Chipseq/unibind_robust_human.tsv")


write.table(df_mm,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            file = "~/Data/Metadata/Chipseq/unibind_mouse.tsv")


write.table(df_rob_mm,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            file = "~/Data/Metadata/Chipseq/unibind_robust_mouse.tsv")
