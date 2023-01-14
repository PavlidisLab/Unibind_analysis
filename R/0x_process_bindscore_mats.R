## Read in the raw bind score matrices, process them, and save out as a list.
## -----------------------------------------------------------------------------


library(tidyverse)
library(preprocessCore)
source("R/00_config.R")
# source("R/Utils/functions.R")


# TODO: remove when regen
bmat_path_hg <- "/space/scratch/amorin/R_objects/cp_unibind_all_scores_hg.RDS"
bmat_path_mm <- "/space/scratch/amorin/R_objects/cp_unibind_all_scores_mm.RDS"

# Unibind metadata
meta_l <- readRDS(meta_outfile)

# Loading of matrices of raw scores
bmat_hg <- readRDS(bmat_path_hg) 
bmat_mm <- readRDS(bmat_path_mm) 

# TODO: remove when regen
bmat_hg <- do.call(cbind, bmat_hg)
bmat_mm <- do.call(cbind, bmat_mm)


# Only considering experiments that have a minimum count of peaks
meta_l <- lapply(meta_l, filter, N_peaks >= min_peaks)

# All matrices to a list
bmat_l <- list(
  Permissive_hg = bmat_hg,
  Robust_hg = bmat_hg[, meta_l$Robust_hg$File],
  Permissive_mm = bmat_mm,
  Robust_mm = bmat_mm[, meta_l$Robust_mm$File]
)


# Unibind includes "duplicated" experiments that have been scored with different
# motifs. Calculate the correlation of binding scores between these experiments
# ------------------------------------------------------------------------------


# Hacky: removing experiments under peak filter can also remove some (but not
# all) experiments from duplicated sets. Need to remove these cases before 
# calculating cor between duplicated experiments.


filter_dupl <- function(meta_df) {
  meta_df <- filter(meta_df, Duplicate)
  n_id <- table(meta_df$ID)
  meta_df <- filter(meta_df, ID %in% names(n_id[n_id > 1]))
}


# Return the unique correlations between duplicate bind scores as a df

get_cor <- function(meta_df, bmat) {
  
  dup_cor <- lapply(unique(meta_df$ID), function(x) {
    
    exps <- filter(meta_df, ID == x)
    mat <- bmat[, exps$File]
    cormat <- cor(mat)
    cors <- unique(cormat[cormat != 1])
    data.frame(ID = x,
               Cor = cors,
               Symbol = exps$Symbol[1])
  })
  
  dup_cor <- do.call(rbind, dup_cor)
  
}


# Return a data frame of summary for experiment correlation for each TF

cor_summary <- function(cor_df) {
  
  tfs <- unique(cor_df$Symbol)
  
  summ_l <- lapply(tfs, function(x) {
    summary(filter(cor_df, Symbol == x)$Cor)
  })
  
  summ_df <- data.frame(Symbol = tfs,
                        do.call(rbind, summ_l))
  
  return(summ_df)
}


# Keeping only duplicated experiment metadata
meta_dup_l <- lapply(meta_l, filter_dupl)


# List of data frames of the correlation between duplicates
cor_l <- lapply(1:length(bmat_l), function(i) {
  get_cor(meta_dup_l[[i]], bmat_l[[i]])
})
names(cor_l) <- names(bmat_l)


# MEIS1, SMAD3, STAT1 examples of TFs with low cor: suggests that they have
# multiple and very distinct motifs on JASPAR. FOSL2 meanwhile most similar
# in all comparisons.

summ_l <- lapply(cor_l, cor_summary)
cor_min <- lapply(summ_l, filter, Median == min(Median))
cor_max <- lapply(summ_l, filter, Median == max(Median))


# Average duplicates in matrix and collapse info in metadata


dup_out <- bmat_hg[, !colnames(bmat_hg) %in% meta_dedup_hg$File]

dup_avg <- lapply(unique(meta_dedup_hg$ID), function(x) {
  cols <- filter(meta_dedup_hg, ID == x)
  mat <- bmat_hg[, cols$File]
  rowMeans(mat)
})

dup_avg <- do.call(cbind, dup_avg)

meta_dedup_hg <- meta_dedup_hg %>% 
  mutate(File = str_replace(File, "\\.MA.*", "\\.avg")) %>% 
  distinct(File, .keep_all = TRUE)

stopifnot(identical(nrow(meta_dedup_hg), ncol(dup_avg)))

colnames(dup_avg) <- meta_dedup_hg$File


# Combine average and dup out to get de-duplicated matrix and meta


bmat_dedup_hg <- cbind(dup_out, dup_avg)

meta_final_hg <- filter(meta_hg, !(ID %in% meta_dedup_hg$ID)) %>% 
  rbind(meta_dedup_hg) %>% 
  arrange(Symbol)


bmat_dedup_hg <- bmat_dedup_hg[, meta_final_hg$File]


stopifnot(identical(
  bmat_dedup_hg[, "EXP036852.renal_tubular_cells.ARNT.avg"],
  rowMeans(bmat_hg[, c("EXP036852.renal_tubular_cells.ARNT.MA0004.1.damo.bed",
                       "EXP036852.renal_tubular_cells.ARNT.MA0006.1.damo.bed",
                       "EXP036852.renal_tubular_cells.ARNT.MA0259.1.damo.bed")])
))


# Plots
# ------------------------------------------------------------------------------


# Boxplots of the binding score correlations of duplicates, grouped by TF

cor_boxplot <- function(cor_df) {
  
  ggplot(cor_df, aes(x = Symbol, y = Cor)) +
    geom_boxplot() +
    theme_classic() +
    ylab("Pcor of binding scores") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 25),
          axis.text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
}


plot_l1 <- lapply(cor_l, cor_boxplot)
