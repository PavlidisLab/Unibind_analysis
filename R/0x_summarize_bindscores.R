##
## TODO: replace meta, bmat, GR, DHS with regen
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")
source("R/Utils/functions.R")

bind_hg <- readRDS("~/scratch/R_objects/unibind_bindscore_human.RDS")
bind_mm <- readRDS("~/scratch/R_objects/unibind_bindscore_mouse.RDS")
# bind_l <- list(Human = bind_hg$Mat_QNL, Mouse = bind_mm$Mat_QNL)
bind_l <- list(Human = bind_hg$Mat_raw, Mouse = bind_mm$Mat_raw)

# Load protein coding genes and convert to GRanges
pc_hg <- pc_to_gr(read.delim(ref_hg, stringsAsFactors = FALSE))
pc_mm <- pc_to_gr(read.delim(ref_mm, stringsAsFactors = FALSE))

#
gr_hg <- readRDS("~/scratch/R_objects/unibind_grlist_perm_human.RDS")
gr_mm <- readRDS("~/scratch/R_objects/unibind_grlist_perm_mouse.RDS")

#
dhs_hg <- read.delim("~/Data/Chromosome_info/ENCODE_human_DHS.tsv", stringsAsFactors = FALSE)
dhs_hg <- mutate(dhs_hg, seqname = str_replace(seqname, "chr", ""))
dhs_hg <- makeGRangesFromDataFrame(dhs_hg, keep.extra.columns = TRUE)


# Use GR to generate a gene x experiment matrix where elements are the distance
# in bp between the gene TSS and the nearest peak in the given experiment. This
# is used to explore peak distances for commonly and rarely bound genes.
# ------------------------------------------------------------------------------


nearest_dist_mat <- function(pc, gr_l, ncores = cores) {
  
  dist_vec <- rep(NA, length(pc$Symbol))
  
  dist_l <- mclapply(gr_l, function(x) {
    dist <- distanceToNearest(pc, x)
    dist_vec[dist@from] <- dist@elementMetadata$distance
    return(dist_vec)
  })
  
  mat <- do.call(cbind, dist_l)
  colnames(mat) <- names(gr_l)
  rownames(mat) <- pc$Symbol
  return(mat)
}


# dist_hg <- nearest_dist_mat(pc = pc_hg, gr_l = gr_hg)
# dist_mm <- nearest_dist_mat(pc = pc_mm, gr_l = gr_mm)
# 
# saveRDS(list(Human = dist_hg, Mouse = dist_mm),
#         "~/scratch/R_objects/unibind_dist_to_nearest_mats.RDS")


# Generate average binding scores across experiments: 1) Across all experiments;
# 2) Grouped by TF; 3) Mean of means across TFs, to summarize across all 
# experiments but to reduce effect of imbalanced TF counts.
# ------------------------------------------------------------------------------


# Return a df of the average bind score per gene across all experiments

get_all_mean <- function(mat) {
  data.frame(Symbol = rownames(mat), Mean = rowMeans(mat))
}


# Return a gene by TF matrix of the average binding scores

get_tf_mean <- function(mat, meta) {
  
  tfs <- unique(meta$Symbol)
  
  mean_l <- lapply(unique(tfs), function(x) {
    meta <- filter(meta, Symbol == x)
    rowMeans(mat[, meta$File, drop = FALSE])
  })
  
  mean_mat <- as.matrix(do.call(cbind, mean_l))
  colnames(mean_mat) <- tfs
  
  return(mean_mat)
  
}


# Mean summaries into a list

mean_l <- list(
  Human_all = get_all_mean(bind_l$Human),
  Human_TF = get_tf_mean(bind_l$Human, bind_hg$Meta),
  Human_mom = get_all_mean(get_tf_mean(bind_l$Human, bind_hg$Meta)),
  Mouse_all = get_all_mean(bind_l$Mouse),
  Mouse_TF = get_tf_mean(bind_l$Mouse, bind_mm$Meta),
  Mouse_mom = get_all_mean(get_tf_mean(bind_l$Mouse, bind_mm$Meta))
)



# plot(mean_l$Human_all$Mean, mean_l$Human_mom$Mean)
# cor(mean_l$Human_all$Mean, mean_l$Human_mom$Mean)
# plot(mean_l$Mouse_all$Mean, mean_l$Mouse_mom$Mean)
# cor(mean_l$Mouse_all$Mean, mean_l$Mouse_mom$Mean)


# Top bound genes, using mean of means.
# ------------------------------------------------------------------------------


filter_top <- function(mean_df, qtl = 0.99) {
  filter(mean_df, Mean >= quantile(Mean, qtl)) %>% arrange(desc(Mean))
}


topbound_hg <- filter_top(mean_l$Human_mom)




# Rarely/never bound genes
# ------------------------------------------------------------------------------


mat_raw <- bind_hg$Mat_raw
mat_qnl <- bind_hg$Mat_QNL

mean_raw <- get_all_mean(get_tf_mean(mat_raw, bind_hg$Meta))
mean_qnl <- get_all_mean(get_tf_mean(mat_qnl, bind_hg$Meta))

mean_raw %>% head
mean_qnl %>% head

plot(mean_raw$Mean, mean_qnl$Mean)
cor(mean_raw$Mean, mean_qnl$Mean)


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


p1a <- plot_hist(mean_l$Human_mom, "Human")
p1b <- plot_hist(mean_l$Mouse_mom, "Mouse")
