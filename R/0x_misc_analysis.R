## This script is a placeholder for work that was started but is incomplete, or
## for quick hits that were written for an earlier draft of this project. It
## therefore contains (messy) code that currently will not run. 
## -----------------------------------------------------------------------------


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# The following was written to begin exploring the Unibind bind score matrices,
# getting summaries of binding and exploring the similarity of experiments.
# ------------------------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(vegan)
library(ggrepel)
library(umap)
library(factoextra)
# library(dendextend)
# library(cluster)
source("~/regnetR/R/utils/bscore_unibind_functions.R")


# list objects contain the data matrices [robust/permissive]
bind_hg <- readRDS("~/scratch/R_objects/unibind_bindscore_robust_human.RDS")
bind_mm <- readRDS("~/scratch/R_objects/unibind_bindscore_robust_mouse.RDS")
# bind_hg <- readRDS("~/scratch/R_objects/unibind_bindscore_human.RDS")
# bind_mm <- readRDS("~/scratch/R_objects/unibind_bindscore_mouse.RDS")

# ortho genes
pc_ortho <- read.delim("~/Data/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv", stringsAsFactors = FALSE)

# TFs common to both species
tf_common <- intersect(names(bind_hg$Group_mean), names(bind_mm$Group_mean))


# Get matrix of mean TR binding scores
# TODO: basic summary
# ------------------------------------------------------------------------------


get_mean_mat <- function(bind_list) {
  mean_mat <- do.call(cbind, lapply(bind_list$Group_mean, `[[`, "Mean"))
  rownames(mean_mat) <- bind_list$Group_mean[[1]]$Symbol
  return(mean_mat)
}



mmat_l <- list(
  Human = get_mean_mat(bind_hg),
  Mouse = get_mean_mat(bind_mm)
)


# Get the average bind score for each gene across experiments. Want both the
# flat average, which may be influenced by uneven sample sizes, as well as mean
# of means. Don't want a typical weighted mean, as this would give highest 
# weights to CTCF. Instead want all TF groups treated equally
# ------------------------------------------------------------------------------


allmean_l <- list(
  Human_flat = bind_hg$All_mean,
  Human_mofm = rownames_to_column(data.frame(Mean = rowMeans(mmat_l$Human)), var = "Symbol"),
  Mouse_flat = bind_mm$All_mean,
  Mouse_mofm = rownames_to_column(data.frame(Mean = rowMeans(mmat_l$Mouse)), var = "Symbol")
)


identical(allmean_l$Human_flat$Symbol, allmean_l$Human_mofm$Symbol)
plot(allmean_l$Human_flat$Mean, allmean_l$Human_mofm$Mean)
cor(allmean_l$Human_flat$Mean, allmean_l$Human_mofm$Mean)


# Inspect top bound across TFs

head(arrange(allmean_l$Human_mofm, desc(Mean)), 30)
head(arrange(allmean_l$Mouse_mofm, desc(Mean)), 30)


# How many never bound? NOTE - binding score accepting far distances + quantile
# norm makes this a bit fuzzy. using a simple and arbitrary cutoff based on min
# value to grab genes that are "essentially never bound." may be desirable to 
# instead pick a window (like the cutoff for bind score) and look for genes that
# never have a discrete peak in that window

fuzz <- 1.1
filter(allmean_l$Human_mofm, Mean <= min(Mean) * fuzz)
filter(allmean_l$Mouse_mofm, Mean <= min(Mean) * fuzz)




# Compare the mean bind score from all experiments for ortho genes between
# mouse and human. Find strong correlation: ortho genes frequently found to have
# signal across human ChIP-seq experiments also tend to have signal in across
# mouse experiments
# ------------------------------------------------------------------------------


allmean_hg <- allmean_l$Human_mofm %>% 
  filter(Symbol %in% pc_ortho$Symbol_hg) %>% 
  arrange(match(Symbol, pc_ortho$Symbol_hg)) %>% 
  dplyr::rename(Mean_hg = Mean) %>% 
  mutate(ID = pc_ortho$ID)

allmean_mm <- allmean_l$Mouse_mofm %>% 
  filter(Symbol %in% pc_ortho$Symbol_mm) %>% 
  arrange(match(Symbol, pc_ortho$Symbol_mm)) %>% 
  dplyr::rename(Mean_mm = Mean) %>% 
  mutate(ID = pc_ortho$ID)

allmean_ortho <- left_join(allmean_hg[, c("Mean_hg", "ID")],
                           allmean_mm[, c("Mean_mm", "ID")],
                           by = "ID") %>% 
  mutate(Rank_hg = rank(-Mean_hg),
         Rank_mm = rank(-Mean_mm),
         Diff_rank = Rank_hg - Rank_mm)


cor(allmean_ortho$Mean_hg, allmean_ortho$Mean_mm)



# Get correlation of mean binding scores between TFs, the pairwise df of this mat,
# and a hclust objects of the mean binding mats
# ------------------------------------------------------------------------------


# matrix for clustering: rows/observations/TFs x cols/features/genes

hclust_mat <- function(mat,
                       scale = FALSE,
                       dis = "pearson", 
                       linkage = "complete") {
  
  if (scale) {
    mat <- scale(mat)
  }
  
  # dist: distance between rows/observations/TFs
  dist <- factoextra::get_dist(mat, method = dis)
  hc <- hclust(dist, method = linkage)
  return(hc)
}



cmat_l <- lapply(mmat_l, cor)
cdf_l <- lapply(cmat_l, mat_to_df)
# tcmat_l <- lapply(cmat_l, tri_cormat)



# For hclust focus on pearson distance, no scaling of mean mat, and complete
# linkage (most extreme dissimilar pair as cluster distance)

hc_l <- list(Human = hclust_mat(t(mmat_l$Human)),
             Human_common = hclust_mat(t(mmat_l$Human[pc_ortho$Symbol_hg, tf_common])),
             Mouse = hclust_mat(t(mmat_l$Mouse)),
             Mouse_common = hclust_mat(t(mmat_l$Mouse[pc_ortho$Symbol_mm, tf_common])))


plot(hc_l$Human, hang = -1)
plot(hc_l$Human_common, hang = -1)
plot(hc_l$Mouse, hang = -1)
plot(hc_l$Mouse_common, hang = -1)




# tt <- lapply(2:(length(tf_common) - 1), function(x) {
tt <- lapply(2:30, function(x) {
  a <- cutree(hc_l$Human_common, k = x)
  b <- cutree(hc_l$Mouse_common, k = x)
  ri <- mclust::adjustedRandIndex(a, b)
  # tab <- table(a, b)
  # which(tab == 1, arr.ind = TRUE)
})

which.max(unlist(tt))

a <- rep(1:3, 4)
a
b <- rep(c("A", "B", "C", "D"), 3)
b
mclust::adjustedRandIndex(a, b)


# Summarize the range of inter-TF cor for each TF. Order by median cor
# ------------------------------------------------------------------------------


# Return col vector named str from mat, excluding the row element named str.
# Eg, removing self entry from a correlation matrix

excl_vec <- function(str, mat) {
  vec <- mat[, str]
  vec <- vec[names(vec) != str]
  return(vec)
}


# Return a dataframe of correlation ranges for each TF in cmat

range_cor <- function(cmat) {
  
  cor_l <- lapply(colnames(cmat), function(x) {
    
    vec <- excl_vec(x, cmat)
    summ <- summary(vec)
    
    data.frame(
      Symbol = x,
      Group = c("Min", "Q1", "Med", "Q3", "Max"),
      Value = c(summ[[1]], summ[[2]], summ[[3]], summ[[5]], summ[[6]])
    )
  })
  
  cor_df <- do.call(rbind, cor_l)
  
  med_df <- cor_df %>% 
    filter(Group == "Med") %>% 
    group_by(Symbol) %>%
    arrange(Value)
  
  cor_df <- cor_df[order(match(cor_df$Symbol, med_df$Symbol)), ]
  cor_df$Symbol <- factor(cor_df$Symbol, levels = unique(cor_df$Symbol))
  
  return(cor_df)
}



crange_l <- lapply(cmat_l, range_cor)


# TFs with the highest median correlation. Note that the median cor is itself
# correlated between species. IE, TFs that tend to have similar binding profiles
# with other TFs in human tend to have also have this quality in mouse


crange_l$Human %>% 
  filter(Group == "Med") %>% 
  arrange(desc(Value)) %>% 
  head(10)

crange_l$Mouse %>% 
  filter(Group == "Med") %>% 
  arrange(desc(Value)) %>% 
  head(10)


medcor_ortho <- left_join(
  filter(crange_l$Human, Group == "Med"),
  filter(crange_l$Mouse, Group == "Med"),
  by = "Symbol")

cor(medcor_ortho$Value.x, medcor_ortho$Value.y, use = "pairwise.complete.obs")


# For each TF, get the TF for which it has the highest cor
# ------------------------------------------------------------------------------


get_topcor <- function(cmat) {
  
  topcor <- lapply(colnames(cmat), function(x) {
    vec <- excl_vec(x, cmat)
    topcor <- head(sort(vec, decreasing = TRUE), 1)
    
    data.frame(TF1 = x,
               TF2 = names(topcor),
               Cor = topcor,
               row.names = NULL)
  })
  do.call(rbind, topcor)
}

topc_l <- lapply(cmat_l, get_topcor)


# For top pairs, how highly ranked is that pair in the other species?


tf_ortho <- intersect(names(bind_hg$Group_mean), names(bind_mm$Group_mean))
stopifnot(all(tf_ortho %in% pc_ortho$Symbol_hg))



top_xspecies <- lapply(tf_ortho, function(x) {
  
  top_hg <- filter(topc_l$Human, TF1 == x)$TF2
  top_mm <- filter(topc_l$Mouse, TF1 == x)$TF2
  
  vec_hg <- sort(excl_vec(x, cmat_l$Human), decreasing = TRUE)
  hg_rank <- which(names(vec_hg) == top_mm)
  
  vec_mm <- sort(excl_vec(x, cmat_l$Mouse), decreasing = TRUE)
  mm_rank <- which(names(vec_mm) == top_hg)
  
  data.frame(
    Symbol = x,
    Top_cor_hg = top_hg,
    Top_cor_mm = top_mm,
    Top_mouse_rank_in_human = ifelse(length(hg_rank) > 0, hg_rank, NA),
    Top_human_rank_in_mouse = ifelse(length(mm_rank) > 0, mm_rank, NA)
  )
})


top_xspecies <- do.call(rbind, top_xspecies)

sum(top_xspecies$Top_mouse_rank_in_human == 1, na.rm = TRUE)/nrow(top_xspecies)
sum(top_xspecies$Top_human_rank_in_mouse == 1, na.rm = TRUE)/nrow(top_xspecies)


# Get gene x TF binary matrix of differential binding status. Compare with
# and without a logFC filter. Describe the DB counts across genes and TFs
# TODO: organize comments on LFC
# Find dramatic drop in DB counts for factors like CTCF, AR, ESRG1... uncertain
# if this is data-driven (these have lots of samples) or biological (these
# bind many regions)
# ------------------------------------------------------------------------------


get_db_mat <- function(fit_list, fdr_cutoff = 0.05, lfc_cutoff = 1) {
  
  gene_vec <- rep(0, length(rownames(fit_list[[1]])))
  names(gene_vec) <- sort(rownames(fit_list[[1]]))
  
  db_list <- lapply(names(fit_list), function(x) {
    
    diffbound <- fit_list[[x]] %>% 
      filter(adj.P.Val < fdr_cutoff & logFC >= lfc_cutoff) %>% 
      rownames()
    
    gene_vec[diffbound] <- 1
    return(gene_vec)
  })
  
  db_mat <- do.call(cbind, db_list)
  colnames(db_mat) <- names(fit_list)
  return(db_mat)
  
}



count_df <- function(mat1, mat2) {
  
  df1 <- data.frame(Count = rowSums(mat1)) %>%
    rownames_to_column(var = "Symbol")
  
  df2 <- data.frame(Count = rowSums(mat2)) %>%
    rownames_to_column(var = "Symbol")
  
  count_df <- left_join(df1, df2, by = "Symbol")
  colnames(count_df) <- c("Symbol", "LFC0", "LFC1")
  
  return(count_df)
}



dbm_l <- list(
  Human0 = get_db_mat(bind_hg$Fit, lfc_cutoff = 0),
  Human1 = get_db_mat(bind_hg$Fit, lfc_cutoff = 1),
  Mouse0 = get_db_mat(bind_mm$Fit, lfc_cutoff = 0),
  Mouse1 = get_db_mat(bind_mm$Fit, lfc_cutoff = 1)
)


count_l <- list(
  Human_gene = count_df(dbm_l$Human0, dbm_l$Human1),
  Human_TF = count_df(t(dbm_l$Human0), t(dbm_l$Human1)),
  Mouse_gene = count_df(dbm_l$Mouse0, dbm_l$Mouse1),
  Mouse_TF = count_df(t(dbm_l$Mouse0), t(dbm_l$Mouse1))
)


# Jaccard of DB genes for every TR pair
# TODO: describe
# ------------------------------------------------------------------------------


get_jac_mat <- function(db_mat) {
  
  all0_tf <- which(colSums(db_mat) == 0)
  all0_gene <- which(rowSums(db_mat) == 0)
  
  if (length(all0_tf) > 0) {
    db_mat <- db_mat[, -all0_tf]
  }
  
  if (length(all0_gene) > 0) {
    db_mat <- db_mat[-all0_gene,]
  }
  
  db_mat <- t(db_mat)  # vegdist works on rows, want cols/TFs compared
  
  jac_mat <- as.matrix(vegdist(
    db_mat,
    method = "jaccard",
    upper = FALSE,
    diag = TRUE,
    binary = "TRUE"
  ))
  
  return(1 - jac_mat)
  
}


jmat_l <- lapply(dbm_l, get_jac_mat)
jdf_l <- lapply(jmat_l, mat_to_df)


# Looking at PCA and UMAP of TF profiles
# ------------------------------------------------------------------------------


remove_nonvar <- function(mat) {
  # Remove rows of matrix that show no variation
  no_sd <- which(apply(mat, 1, sd) == 0)
  if (length(no_sd) > 1) {
    mat <- mat[-no_sd, ]
  }
  return(mat)
}


pca_and_var <- function(mat) {
  
  # Performs PCA with prcomp and returns list of the resulting 
  # object as well as the variance explained
  
  # prcomp expects samples as rows, features (genes) as columns so transpose
  pcmat <- prcomp(t(mat), scale = TRUE)
  
  # variance explained by the PCs
  prc_var <- pcmat$sdev ^ 2
  var_explained <- round(prc_var / sum(prc_var) * 100, 2)
  cumvar_explained <- cumsum(var_explained)/sum(var_explained)
  return(list(PC = pcmat, 
              Var_explained = var_explained, 
              Cumvar_explained = cumvar_explained))
}


pc_l <- lapply(mmat_l, function(x) pca_and_var(remove_nonvar(x)))

umap_l <- lapply(mmat_l, function(x) umap(t(x)))


# Plots
# ------------------------------------------------------------------------------


# Scatter of all mean for ortho genes between mouse and human

cor(allmean_ortho$Mean_hg, allmean_ortho$Mean_mm)

p1 <- 
  ggplot(allmean_ortho, aes(x = Mean_hg, y = Mean_mm)) +
  # geom_point(size = 3, alpha = 0.4, shape = 21, fill = "slategrey") +
  geom_point(size = 3, alpha = 0.4, shape = 21, fill = "#bcbddc") +
  geom_text_repel(
    data = filter(allmean_ortho, Rank_hg < 50 & Rank_mm < 50),
    aes(x = Mean_hg, y = Mean_mm, label = str_replace(ID, "_.*", "")),
    force = 1,
    force_pull = 2,
    size = 5,
    max.overlaps = 20
  ) +
  xlab("Human aggregated binding score") +
  ylab("Mouse aggregated binding score") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))




# Boxplot of cor ranges


p2a <- 
  ggplot(crange_l$Human, aes(x = Symbol, y = Value)) +
  geom_boxplot() +
  ylab("Pearson correlation") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20))

p2b <- 
  ggplot(crange_l$Mouse, aes(x = Symbol, y = Value)) +
  geom_boxplot() +
  ylab("Pearson correlation") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20))


p2 <- plot_grid(p2a, p2b, ncol = 1)


# scatter plot of median correlation each TF has with the rest of the TFs


p3 <- ggplot(medcor_ortho, aes(x = Value.x, y = medcor_ortho$Value.y)) +
  geom_point() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))



# heatmap of cor
pal_length <- 100
bluered_pal <- colorRampPalette(c("#0571b0", "white", "#ca0020"))(pal_length)
# color_breaks <- seq(min(cmat_l$Human, na.rm = TRUE), max(cmat_l$Human, na.rm = TRUE), length.out = pal_length)
color_breaks <- seq(-1, 1, length.out = pal_length)


pheatmap(tcmat_l$Human, 
         cluster_col = FALSE, 
         cluster_row = FALSE,
         show_colnames = FALSE,
         color = bluered_pal,
         breaks = color_breaks,
         na_col = "white",
         border_color = NA,
         fontsize = 5)

pheatmap(cmat_l$Human,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = bluered_pal,
         breaks = color_breaks,
         fontsize = 5)


# Plots comparing +/- LFC

p1 <- 
  ggplot(count_l$Human_TF, aes(x = LFC0, y = LFC1)) +
  geom_jitter(shape = 21, colour = "black", fill = "royalblue", size = 3) +
  theme_classic() +
  ggtitle("Count of diff bound per TF") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

p2 <- 
  ggplot(count_l$Human_gene, aes(x = LFC0, y = LFC1)) +
  geom_jitter(shape = 21, alpha = 0.3, fill = "goldenrod", size = 3) +
  theme_classic() +
  ggtitle("Count of diff bound per gene") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

plot_grid(p1, p2, ncol = 1)


hist(count_l$Human_gene$LFC0, breaks = 100)
hist(count_l$Human_gene$LFC1, breaks = 100)
hist(count_l$Human_TF$LFC0, breaks = 100)
hist(count_l$Human_TF$LFC1, breaks = 100)

hist(count_l$Mouse_gene$LFC0, breaks = 100)
hist(count_l$Mouse_gene$LFC1, breaks = 100)
hist(count_l$Mouse_TF$LFC0, breaks = 100)
hist(count_l$Mouse_TF$LFC1, breaks = 100)


# 
ctcf_df <- bind_hg$Fit$CTCF %>% 
  rownames_to_column(var = "Symbol")

min_apval <- sort(ctcf_df$adj.P.Val)
min_apval <- head(min_apval[min_apval != 0], 1)

ctcf_df$adj.P.Val <- ifelse(ctcf_df$adj.P.Val == 0, min_apval, ctcf_df$adj.P.Val)

ggplot(ctcf_df, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(shape = 21, alpha = 0.4, size = 3, fill = "slategrey") +
  geom_hline(yintercept = -log10(0.05), colour = "red") +
  geom_vline(xintercept = 1, colour = "blue") +
  ylab("-log10(adj Pval)") +
  ggtitle("Human CTCF") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))


gene <- "GBX1"  # example of gene with LFC=0.96 (so would lose) and highly sig 
boxplot(bind_hg$Mat_QNL[gene, ] ~ bind_hg$Meta$Symbol)
boxplot(bind_hg$Mat_QNL[gene, ] ~ bind_hg$Meta$Symbol == "CTCF", cex.lab = 1.5)


# plot of PCs

data.frame(pc_l$Human$PC$x) %>% 
  rownames_to_column(var = "Symbol") %>%
  mutate(Group = Symbol %in% c("CEBPA", "CEBPB")) %>% 
  ggplot(., aes(x = PC1, y = PC2, fill = Group)) +
  geom_point(shape = 21, colour = "black", size = 3) +
  theme_classic()

# Plot of UMAP

udf_hg <- data.frame(umap_l$Human$layout) %>% 
  rownames_to_column(var = "Symbol") 

udf_mm <- data.frame(umap_l$Mouse$layout) %>% 
  rownames_to_column(var = "Symbol") 


ggplot(udf_mm, aes(x = X1, y = X2)) +
  geom_point(size = 3, shape = 21, fill = "#31a354") +
  geom_text_repel(
    # data = filter(udf_hg, X1 < -2),
    aes(x = X1, y = X2, label = Symbol),
    force = 1,
    force_pull = 2,
    size = 5,
    max.overlaps = 20
  ) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

# PC1 vs UMAP 1
plot(umap_l$Mouse$layout[, 1], pc_l$Mouse$PC$x[, 1])



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# The following was written to compare results obtained using the permissive
# versus robust collection of Unibind data. It is outdated and will not run.
# ------------------------------------------------------------------------------


library(tidyverse)
library(reshape2)
library(parallel)
library(pheatmap)
library(RColorBrewer)
source("~/regnetR/R/utils/bscore_unibind_functions.R")


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



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# The following was written to compare the aggregated rankings from the ENCODE
# pipeline to Unibind. This is outdated and will not run.
# ------------------------------------------------------------------------------

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
