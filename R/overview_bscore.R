## Analysis of Unibind bind matrices
## TODO: organize functions
# TODO: explicit clustering here?
# TODO: describe range
# TODO: better description of top cor pairs across species
## -----------------------------------------------------------------------------

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
