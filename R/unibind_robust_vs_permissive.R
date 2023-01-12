#

library(tidyverse)
library(reshape2)
library(parallel)
library(pheatmap)
library(RColorBrewer)
source("~/regnetR/R/utils/bscore_unibind_functions.R")

#
# tfs_hg <- read.delim("~/Data/Metadata/human_tfs_lambert2018.tsv", stringsAsFactors = FALSE)
pc_hg <- read.delim("~/Data/Metadata/ensembl_human_protein_coding_105.tsv", stringsAsFactors = FALSE)

#
prot <- read.delim("~/Data/Expression_files/HPA/normal_tissue.tsv", stringsAsFactors = FALSE)

#
allrank <- readRDS("~/scratch/R_objects/Apr2022_ranked_target_list.RDS")

#
bind_hg <- readRDS("~/scratch/R_objects/unibind_bindscore_robust_human.RDS")
bind_perm_hg <- readRDS("~/scratch/R_objects/unibind_bindscore_human.RDS")
bind_mm <- readRDS("~/scratch/R_objects/unibind_bindscore_robust_mouse.RDS")
bind_perm_mm <- readRDS("~/scratch/R_objects/unibind_bindscore_mouse.RDS")

#
lt <- read.delim("~/Data/Metadata/Curated_targets_all_Oct2022.tsv")

# lt <- mutate(lt,
#              TF_Symbol = str_to_upper(TF_Symbol),
#              Target_Symbol = str_to_upper(Target_Symbol))



# TFs in permissive but not robust
diff_hg <- setdiff(names(bind_perm_hg$Group_mean), names(bind_hg$Group_mean))
diff_mm <- setdiff(names(bind_perm_mm$Group_mean), names(bind_mm$Group_mean))



# First looking at the concordance of the summarized results from the unibind
# permissive and robust (required presence of canonical motif) collection
# ------------------------------------------------------------------------------


# Compare mean score from robust and permissive
# Human PAX6 example of weaker cor (0.6658) between robust and permissive


bindscore_cor <- function(robust_l, permissive_l) {
  
  cor_l <- lapply(names(robust_l$Group_mean), function(x) {
    rob <- dplyr::rename(robust_l$Group_mean[[x]], Mean_robust = Mean)
    perm <- dplyr::rename(permissive_l$Group_mean[[x]], Mean_permissive = Mean)
    df <- left_join(rob, perm, by = "Symbol")
    cor(df$Mean_robust, df$Mean_permissive)
  })
  
  names(cor_l) <- names(robust_l$Group_mean)
  sort(unlist(cor_l))
}


uni_cor_hg <- bindscore_cor(bind_hg, bind_perm_hg)
uni_cor_mm <- bindscore_cor(bind_mm, bind_perm_mm)


# Count of topn genes in mean bind score between the two collections


topn <- 500


bindscore_topn <- function(robust_l, permissive_l) {
  
  topn_l <- lapply(names(robust_l$Group_mean), function(x) {
    rob <- arrange(robust_l$Group_mean[[x]], desc(Mean))$Symbol[1:500]
    perm <- arrange(permissive_l$Group_mean[[x]], desc(Mean))$Symbol[1:500]
    length(intersect(rob, perm))
  })
  
  names(topn_l) <- names(robust_l$Group_mean)
  sort(unlist(topn_l))
}


uni_topn_hg <- bindscore_topn(bind_hg, bind_perm_hg)
uni_topn_mm <- bindscore_topn(bind_mm, bind_perm_mm)


# Jaccard of diffbind genes for robust and permissive unibind


jacc_diffbind <- function(robust_l, permissive_l) {
  
  jacc_l <- lapply(names(robust_l$Group_mean), function(x) {
    rob <- rownames(filter(robust_l$Fit[[x]], logFC > 0 & adj.P.Val < 0.05))
    perm <- rownames(filter(permissive_l$Fit[[x]], logFC > 0 & adj.P.Val < 0.05))
    length(intersect(rob, perm)) / length(union(rob, perm))
  })
  
  names(jacc_l) <- names(robust_l$Group_mean)
  sort(unlist(jacc_l))
}


uni_jacc_hg <- jacc_diffbind(bind_hg, bind_perm_hg)
uni_jacc_mm <- jacc_diffbind(bind_mm, bind_perm_mm)


# Plots
# ------------------------------------------------------------------------------


# Barchart of count of diff bound genes


plot_barchart <- function(df) {
  
  ggplot(df, aes(x =  reorder(Symbol, N), y = N)) +
    geom_bar(stat = "identity") +
    ylab("Count of differentially bound genes") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1),
          axis.ticks.x = element_blank())
  
}


plot_l <- lapply(ndiff_l, plot_barchart)


# Count of diff bound genes when using robust vs permissive

plot(count_hg$N.x, count_hg$N.y)
plot(count_mm$N.x, count_mm$N.y)
