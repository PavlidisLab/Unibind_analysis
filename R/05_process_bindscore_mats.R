## Read in the raw bind score matrices, process them, and save out as a list.
## TODO: Need to square away discrepancy in meta counts between dat and de-dup'd
## -----------------------------------------------------------------------------

library(tidyverse)
library(preprocessCore)
library(parallel)
source("R/00_config.R")

# Load Unibind metadata and only consider experiments with a min count of peaks
meta_l <- readRDS(meta_outfile)
meta_l <- lapply(meta_l, filter, N_peaks >= min_peaks)

# Loading of matrices of raw scores
bmat_hg <- readRDS(bmat_path_hg) 
bmat_mm <- readRDS(bmat_path_mm) 

# All matrices to a list
bmat_l <- list(
  Permissive_hg = bmat_hg[, meta_l$Permissive_hg$File],
  Robust_hg = bmat_hg[, meta_l$Robust_hg$File],
  Permissive_mm = bmat_mm[, meta_l$Permissive_mm$File],
  Robust_mm = bmat_mm[, meta_l$Robust_mm$File]
)


# Unibind includes "duplicated" experiments that have been scored with different
# motifs. Removing experiments under peak filter can also remove some (but not
# all) experiments from duplicated sets. Split duplicated and non-dup meta.
# ------------------------------------------------------------------------------


filter_dupl <- function(meta_df) {
  meta_df <- filter(meta_df, Duplicate)
  n_id <- table(meta_df$ID)
  meta_df <- filter(meta_df, ID %in% names(n_id[n_id > 1]))
}


meta_dup_l <- lapply(meta_l, filter_dupl)
meta_nodup_l <- lapply(meta_l, filter, !Duplicate)


# Calculate the correlation of binding scores between dup experiments
# ------------------------------------------------------------------------------


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
# ------------------------------------------------------------------------------


# Return a matrix where the duplicated experiments are averaged

avg_dup <- function(meta_df, bmat) {
  
  stopifnot(all(colnames(bmat) %in% meta_df$File))
  
  dup_avg <- lapply(unique(meta_df$ID), function(x) {
    cols <- filter(meta_df, ID == x)
    bmat <- bmat[, cols$File]
    rowMeans(bmat)
  })
  
  dup_avg <- do.call(cbind, dup_avg)
  return(dup_avg)
}


# Return a df where the duplicated experiment info has been collapsed

dedup_meta <- function(meta_df) {
  
  meta_df <- meta_df %>% 
    mutate(File = str_replace(File, "\\.MA.*", "\\.avg")) %>% 
    distinct(File, .keep_all = TRUE)
  
  return(meta_df)
}


# Return a list of the de-duplicated meta and matrix

dedup_data <- function(dup_meta, nodup_meta, dup_mat, nodup_mat) {
  
  avg_dup_mat <- avg_dup(dup_meta, dup_mat)
  avg_dup_meta <- dedup_meta(dup_meta)
  colnames(avg_dup_mat) <- avg_dup_meta$File
  final_meta <- rbind(nodup_meta, avg_dup_meta) %>% arrange(Symbol)
  final_mat <- cbind(nodup_mat, avg_dup_mat)[, final_meta$File]
  
  return(list(Meta = final_meta, Mat_raw = final_mat))
  
}


dedup_l <- mclapply(names(bmat_l), function(x) {
  
  dedup_data(dup_meta = meta_dup_l[[x]], 
             nodup_meta = meta_nodup_l[[x]],
             dup_mat = bmat_l[[x]][,  meta_dup_l[[x]]$File],
             nodup_mat = bmat_l[[x]][, meta_nodup_l[[x]]$File])
  
}, mc.cores = 4)
names(dedup_l) <- names(bmat_l)



stopifnot(identical(
  dedup_l$Permissive_hg$Mat[, "EXP036852.renal_tubular_cells.ARNT.avg"],
  rowMeans(bmat_hg[, c("EXP036852.renal_tubular_cells.ARNT.MA0004.1.damo.bed",
                       "EXP036852.renal_tubular_cells.ARNT.MA0006.1.damo.bed",
                       "EXP036852.renal_tubular_cells.ARNT.MA0259.1.damo.bed")])
))


# Log transform + quantile norm the raw binding scores
# ------------------------------------------------------------------------------


# Returns dat_l with an additional list element of the quantile norm(log2(mat+1))

add_qnl <- function(dat_l) {
  
  dat_l <- lapply(dat_l, function(x) {
    x$Mat_qnl <- normalize.quantiles(log2(x$Mat_raw + 1), keep.names = TRUE)
    return(x)
  })
  
  return(dat_l)
}


dedup_l <- add_qnl(dedup_l)


test_cor <- cor(
  dedup_l$Permissive_hg$Mat_raw[, 1],
  dedup_l$Permissive_hg$Mat_qnl[, 1],
  method = "spearman"
)

stopifnot(round(test_cor, 3) == 1)


# Save out


saveRDS(dedup_l, bind_dat_path)



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
