## This script generates ranked summaries of TR-target binding evidence
## -----------------------------------------------------------------------------

library(tidyverse)
library(limma)
library(edgeR)
source("R/utils/functions.R")
source("R/00_config.R")

# Permissive or Robust data collection
collection <- "Permissive"
stopifnot(collection %in% c("Robust", "Permissive"))

# Load de-duplicated data
dat <- readRDS(bind_dat_path)

# Metadata
meta_l <- readRDS(meta_path)

# Isolating data to use (permissive or robust collection)
bind_hg <- dat[[paste0(collection, "_hg")]]$Mat_qnl
bind_mm <- dat[[paste0(collection, "_mm")]]$Mat_qnl
meta_hg <- dat[[paste0(collection, "_hg")]]$Meta
meta_mm <- dat[[paste0(collection, "_mm")]]$Meta

# Load protein coding genes and convert to GRanges
pc_hg <- pc_to_gr(read.delim(ref_path_hg, stringsAsFactors = FALSE))
pc_mm <- pc_to_gr(read.delim(ref_path_mm, stringsAsFactors = FALSE))

bind_summary_path <- file.path(dat_dir, paste0("unibind_", collection, "_bindscore_summary.RDS"))
bind_model_path <- file.path(dat_dir, paste0("unibind_", collection, "_bindscore_modelfit.RDS"))


# Generate average binding scores across experiments: 1) Across all experiments;
# 2) Grouped by TF; 3) Mean of means across TFs, to summarize across all 
# experiments but to reduce effect of imbalanced TF counts.
# ------------------------------------------------------------------------------


# Return a df of the average bind score per gene across all experiments. This
# will be influenced by TFs that have more data

get_all_mean <- function(mat) {
  data.frame(Symbol = rownames(mat), Mean = rowMeans(mat))
}


# Return a gene by TF matrix of the average binding scores for each TF's set 
# of experiments.

get_tf_mean <- function(mat, meta) {
  
  tfs <- unique(meta$Symbol)
  
  mean_l <- lapply(unique(tfs), function(x) {
    meta <- filter(meta, Symbol == x)
    rowMeans(mat[, meta$ID, drop = FALSE])
  })
  
  mean_mat <- as.matrix(do.call(cbind, mean_l))
  colnames(mean_mat) <- tfs
  
  return(mean_mat)
  
}


# Mean summaries into a list: include a groupwise average that gives each
# TF equal weight to the global average, regardless of its number of experiments

mean_l <- list(
  Human_all = get_all_mean(bind_hg),
  Human_TF = get_tf_mean(bind_hg, meta_hg),
  Human_MoM = get_all_mean(get_tf_mean(bind_hg, meta_hg)),
  Mouse_all = get_all_mean(bind_mm),
  Mouse_TF = get_tf_mean(bind_mm, meta_mm),
  Mouse_MoM = get_all_mean(get_tf_mean(bind_mm, meta_mm))
)


# Top bound genes, using mean of means.
# ------------------------------------------------------------------------------


filter_top <- function(mean_df, qtl = 0.99) {
  filter(mean_df, Mean >= quantile(Mean, qtl)) %>% arrange(desc(Mean))
}



# Human
top_bound_hg <- filter_top(mean_l$Human_MoM)

# Mouse
top_bound_mm <- filter_top(mean_l$Mouse_MoM)


# Rarely/never bound genes
# ------------------------------------------------------------------------------


filter_bottom <- function(mean_df, qtl = 0.01) {
  filter(mean_df, Mean <= quantile(Mean, qtl)) %>% arrange(Mean)
}


# Human
btm_bound_hg <- filter_bottom(mean_l$Human_MoM)

# Mouse
btm_bound_mm <- filter_bottom(mean_l$Mouse_MoM)


# Binding specificity model: Use limma voom framework to get the expected
# binding score of a TR-gene interaction, using a "one-versus rest" contrast
# to compare a group of TR experiments against all others, controlling for the
# count of peak for the experiment.
# --
# First prepare raw bind score matrices for limma voom, as it performs quantile
# norm and log. Examining hist to filter genes by raw counts across experiments.
# Mouse more zero-bound genes compared to human (likely due to bad annotations).
# ------------------------------------------------------------------------------


# Human

mat_raw_hg <- dat[[paste0(collection, "_hg")]]$Mat_raw
sum_hg <- rowSums(mat_raw_hg)
# hist(sum_hg, breaks = 1000)
min_hg <- 15
# abline(v = min_hg, col = "red")
# view(sort(sum_hg))

keep_hg <- sum_hg > min_hg
mat_hg <- mat_raw_hg[keep_hg, meta_hg$ID]
# hist(rowSums(mat_hg), breaks = 1000)

# Mouse 

mat_raw_mm <- dat[[paste0(collection, "_mm")]]$Mat_raw
sum_mm <- rowSums(mat_raw_mm)
# hist(sum_mm, breaks = 1000)
min_mm <- 120
# abline(v = min_mm, col = "red")
# view(sort(sum_mm))

keep_mm <- sum_mm > min_mm
mat_mm <- mat_raw_mm[keep_mm, meta_mm$ID]
# hist(rowSums(mat_mm), breaks = 1000)


# Design matrices, voom, limma model
# ~ 0 + Symbol + log(N_peaks)
# Means model (~0) so each TF has a coef corresponding to group mean with no intercept
#-------------------------------------------------------------------------------


design_mat <- function(meta) {
  
  mat <- model.matrix(
    ~ 0 + Symbol + log(N_peaks), data = meta)
  
  rownames(mat) <- meta$ID
  colnames(mat) <- str_replace(colnames(mat), "Symbol", "")
  
  return(mat)
}


design_hg <- design_mat(meta_hg)
design_mm <- design_mat(meta_mm)


# voom data and fit model
# --
# voom is a transformation for count matrices with a mean-variance relationship,
# as seen in the binding matrices. While the bind matrices are not integer
# counts like in RNA-seq, G. Smyth says it is appropriate to use on numeric
# matrices https://support.bioconductor.org/p/45695/. 
#-------------------------------------------------------------------------------


get_fit <- function(dat_mat, design_mat) {
  
  v <- voom(dat_mat, design = design_mat, normalize.method = "quantile")
  fit <- lmFit(v, design = design_mat)
  
  return(fit)
}

message("Fitting human model ", Sys.time())
fit_hg <- get_fit(mat_hg, design_hg)

message("Fitting mouse model ", Sys.time())
fit_mm <- get_fit(mat_mm, design_mm)


# Interested in each TR vs the rest, so must construct the appropriate contrast
# vectors. This is setting 1 for the TR of interest, and -1/(nTR-1) for the 
# other nTR-1 TRs, averaging their contributions as the "out" group. Nuisance
# variables are kept as 0 (the output still corrects for them)
# Law et al., 2020 used as reference 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/ 
#-------------------------------------------------------------------------------


contr_list <- function(meta, fit) {
  
  tfs <- unique(meta$Symbol)
  
  contr_vec <- rep(0, length(colnames(coef(fit))))
  names(contr_vec) <- colnames(coef(fit))
  
  clist <- lapply(tfs, function(x) {
    contr_vec[x] <- 1
    contr_vec[names(contr_vec) %in% setdiff(tfs, x)] <- -(1/(length(tfs) - 1))
    return(contr_vec)
  })
  names(clist) <- tfs
  
  return(clist)
}


clist_hg <- contr_list(meta_hg, fit_hg)
clist_mm <- contr_list(meta_mm, fit_mm)


# For each symbol, get the model estimates for the symbol vs all contrast
#-------------------------------------------------------------------------------

# Apply eBayes + contrast fit to the model and extract toptable for each contr

top_fit <- function(contr_list, fit) {
  lapply(contr_list, function(x) {
    contr_fit <- eBayes(contrasts.fit(fit, x))
    topTable(contr_fit, adjust.method = "BH", number = Inf)
  })
}


message("Top fit human model ", Sys.time())
top_hg <- top_fit(clist_hg, fit_hg)

message("Top fit mouse model ", Sys.time())
top_mm <- top_fit(clist_mm, fit_mm)


# Save out
# ------------------------------------------------------------------------------


saveRDS(mean_l, bind_summary_path)


fit_l <- list(
  Human_fit = fit_hg,
  Human_top = top_hg,
  Mouse_fit = fit_mm,
  Mouse_top = top_mm
)


saveRDS(fit_l, bind_model_path)


# Plots
# ------------------------------------------------------------------------------


# Histogram of mean binding scores

plot_hist <- function(mean_df, species, qtl_upper = 0.99) {
  
  if (species == "Human") {
    fill <- "royalblue"
  } else if (species == "Mouse") {
    fill <- "goldenrod"
  }
  
  ggplot(mean_df, aes(x = Mean)) +
    geom_histogram(alpha = 0.6, fill = fill, bins = 100) +
    geom_vline(xintercept = quantile(mean_df$Mean, qtl_upper), 
               linewidth = 1.2,
               colour = "firebrick3") +
    theme_classic() +
    ylab("Density") +
    xlab("Mean binding score") +
    ggtitle(species) +
    theme(
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 25),
      plot.title = element_text(hjust = 0.5, size = 30),
      plot.margin = margin(10, 20, 10, 10) 
    )
}


p1a <- plot_hist(mean_l$Human_MoM, "Human")
p1b <- plot_hist(mean_l$Mouse_MoM, "Mouse")
