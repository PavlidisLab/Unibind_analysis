## This script compares the aggregated rankings from the Encpipe output to that
## obtained from the robust+permissive collection from unibind
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)

# The final high-throughput TR rankings
allrank <- readRDS("~/scratch/R_objects/Apr2022_ranked_target_list.RDS")

# list objects contain the data matrices [robust/permissive]
bind_hg <- readRDS("~/scratch/R_objects/unibind_bindscore_robust_human.RDS") 
bind_mm <- readRDS("~/scratch/R_objects/unibind_bindscore_robust_mouse.RDS")
bind_perm_hg <- readRDS("~/scratch/R_objects/unibind_bindscore_human.RDS")
bind_perm_mm <- readRDS("~/scratch/R_objects/unibind_bindscore_mouse.RDS")


topn <- 500

# Merging allrank with unibind

tf_rank_hg <- intersect(names(allrank$Human), names(bind_perm_hg$Group_mean))
tf_rank_mm <- intersect(str_to_upper(names(allrank$Mouse)), names(bind_perm_mm$Group_mean))


keep_cols <- c(
  "Symbol",
  "Count_DE",
  "Avg_abs_FC",
  "Mean_bind",
  "Bind_logFC",
  "Bind_adj_Pval",
  "Curated_target",
  "Rank_perturbation",
  "Rank_binding",
  "Rank_integrated"
)



allrank_hg <- lapply(tf_rank_hg, function(x) {
  
  perm_df <- bind_perm_hg$Fit[[x]][, c("logFC", "adj.P.Val")] %>% 
    rownames_to_column(var = "Symbol") %>% 
    left_join(bind_perm_hg$Group_mean[[x]], by = "Symbol")
  colnames(perm_df) <- c("Symbol", "logFC_permissive", "FDR_permissive", "Mean_permissive")
  
  if (!is.null(bind_hg$Fit[[x]])) {
    
    rob_df <- bind_hg$Fit[[x]][, c("logFC", "adj.P.Val")] %>% 
      rownames_to_column(var = "Symbol") %>% 
      left_join(bind_hg$Group_mean[[x]], by = "Symbol")
    colnames(rob_df) <- c("Symbol", "logFC_robust", "FDR_robust", "Mean_robust")
    
  } else {
    rob_df <- perm_df[, "Symbol", drop = FALSE]
  }
  
  rank_df <- allrank$Human[[x]][, keep_cols]
  
  rank_df <- dplyr::rename(rank_df, 
                           Mean_encpipe = Mean_bind,
                           logFC_encpipe = Bind_logFC,
                           FDR_encpipe = Bind_adj_Pval)
  
  
  all_df <- left_join(rank_df, rob_df, by = "Symbol") %>% 
    left_join(perm_df, by = "Symbol")
  
})
names(allrank_hg) <- tf_rank_hg



allrank_mm <- lapply(tf_rank_mm, function(x) {
  
  perm_df <- bind_perm_mm$Fit[[x]][, c("logFC", "adj.P.Val")] %>% 
    rownames_to_column(var = "Symbol") %>% 
    left_join(bind_perm_mm$Group_mean[[x]], by = "Symbol")
  colnames(perm_df) <- c("Symbol", "logFC_permissive", "FDR_permissive", "Mean_permissive")
  
  if (!is.null(bind_mm$Fit[[x]])) {
    
    rob_df <- bind_mm$Fit[[x]][, c("logFC", "adj.P.Val")] %>% 
      rownames_to_column(var = "Symbol") %>% 
      left_join(bind_mm$Group_mean[[x]], by = "Symbol")
    colnames(rob_df) <- c("Symbol", "logFC_robust", "FDR_robust", "Mean_robust")
    
  } else {
    rob_df <- perm_df[, "Symbol", drop = FALSE]
  }
  
  rank_df <- allrank$Mouse[[str_to_title(x)]][, keep_cols]
  
  rank_df <- dplyr::rename(rank_df, 
                           Mean_encpipe = Mean_bind,
                           logFC_encpipe = Bind_logFC,
                           FDR_encpipe = Bind_adj_Pval)
  
  
  all_df <- left_join(rank_df, rob_df, by = "Symbol") %>% 
    left_join(perm_df, by = "Symbol")
  
})
names(allrank_mm) <- tf_rank_mm


# Get cor and topn of encpipe vs robust/permissive


compare_hg <- lapply(tf_rank_hg, function(x) {
  
  df <- allrank_hg[[x]]
  
  cor_perm <- cor(df$Mean_permissive, df$Mean_encpipe, use = "pairwise.complete.obs")
  
  topn_perm <- length(intersect(
    arrange(df, desc(Mean_encpipe))$Symbol[1:topn],
    arrange(df, desc(Mean_permissive))$Symbol[1:topn]
  ))
  
  if (x %in% names(bind_hg$Group_mean)) {
    
    cor_rob <- cor(df$Mean_robust, df$Mean_encpipe, use = "pairwise.complete.obs")
    
    topn_rob <- length(intersect(
      arrange(df, desc(Mean_encpipe))$Symbol[1:topn],
      arrange(df, desc(Mean_robust))$Symbol[1:topn]
    ))
    
  } else {
    
    cor_rob <- topn_rob <- NA
  }
  
  data.frame(
    Cor_robust = cor_rob,
    Cor_permissive = cor_perm,
    Topn_robust = topn_rob,
    Topn_permissive = topn_perm
  )
  
})
names(compare_hg) <- tf_rank_hg

compare_hg <- do.call(rbind, compare_hg)
rownames(compare_hg) <- paste0(rownames(compare_hg), "_Human")



compare_mm <- lapply(tf_rank_mm, function(x) {
  
  df <- allrank_mm[[x]]
  
  cor_perm <- cor(df$Mean_permissive, df$Mean_encpipe, use = "pairwise.complete.obs")
  
  topn_perm <- length(intersect(
    arrange(df, desc(Mean_encpipe))$Symbol[1:topn],
    arrange(df, desc(Mean_permissive))$Symbol[1:topn]
  ))
  
  if (x %in% names(bind_mm$Group_mean)) {
    
    cor_rob <- cor(df$Mean_robust, df$Mean_encpipe, use = "pairwise.complete.obs")
    
    topn_rob <- length(intersect(
      arrange(df, desc(Mean_encpipe))$Symbol[1:topn],
      arrange(df, desc(Mean_robust))$Symbol[1:topn]
    ))
    
  } else {
    
    cor_rob <- topn_rob <- NA
  }
  
  data.frame(
    Cor_robust = cor_rob,
    Cor_permissive = cor_perm,
    Topn_robust = topn_rob,
    Topn_permissive = topn_perm
  )
  
})
names(compare_mm) <- tf_rank_mm

compare_mm <- do.call(rbind, compare_mm)
rownames(compare_mm) <- paste0(rownames(compare_mm), "_Mouse")


# Heatmap/table of ENCODE vs unibind


mat_cor <- round(rbind(compare_hg, compare_mm)[, 1:2], 3)
mat_cor <- mat_cor[sort(rownames(mat_cor)),]
colnames(mat_cor) <- c("Robust", "Permissive")

mat_topn <- round(rbind(compare_hg, compare_mm)[, 3:4], 3)
mat_topn <- mat_topn[sort(rownames(mat_topn)),]
colnames(mat_topn) <- c("Robust", "Permissive")


pal_length <- 100
bluered_pal <- colorRampPalette(c("#0571b0", "white", "#ca0020"))(pal_length)
color_breaks <- seq(0, 1, length.out = pal_length)
green_pal <- c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45')



pheatmap(mat_cor,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = bluered_pal,
         breaks = color_breaks,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 30,
         cellwidth = 50,
         cellheight = 50,
         angle_col = 90,
         gaps_col = 1)


pheatmap(mat_topn,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = green_pal,
         display_numbers = TRUE,
         number_format = "%.0f",
         number_color = "black",
         fontsize = 30,
         cellwidth = 50,
         cellheight = 50,
         angle_col = 90,
         gaps_col = 1)
