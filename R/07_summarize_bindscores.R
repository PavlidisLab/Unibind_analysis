## This script generates ranked summaries of TR-target binding evidence
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")
source("R/Utils/functions.R")

# Load de-duplicated data
dat <- readRDS(bind_dat_path)

# Isolating data to use (permissive or robust collection) 
bind_l <- list(Human = dat$Permissive_hg$Mat_qnl, Mouse = dat$Permissive_mm$Mat_qnl)
# bind_l <- list(Human = dat$Permissive_hg$Mat_raw, Mouse = dat$Permissive_mm$Mat_raw)

# TODO: consolidate metas across GR/mats
# meta_l <- list(Human = dat$Permissive_hg$Meta, Mouse = dat$Permissive_mm$Meta)
meta_l <- readRDS(meta_outfile)
meta_hg <- meta_l$Permissive_hg
meta_mm <- meta_l$Permissive_mm

# Load protein coding genes and convert to GRanges
pc_hg <- pc_to_gr(read.delim(ref_hg, stringsAsFactors = FALSE))
pc_mm <- pc_to_gr(read.delim(ref_mm, stringsAsFactors = FALSE))

        
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
  Human_TF = get_tf_mean(bind_l$Human, dat$Permissive_hg$Meta),
  Human_mom = get_all_mean(get_tf_mean(bind_l$Human, dat$Permissive_hg$Meta)),
  Mouse_all = get_all_mean(bind_l$Mouse),
  Mouse_TF = get_tf_mean(bind_l$Mouse, dat$Permissive_mm$Meta),
  Mouse_mom = get_all_mean(get_tf_mean(bind_l$Mouse, dat$Permissive_mm$Meta))
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



# Human
top_bound_hg <- filter_top(mean_l$Human_mom)
top_dist_hg <- dist_hg[top_bound_hg$Symbol, ]
top_summ_hg <- dist_summ_hg[top_bound_hg$Symbol, ]

# Mouse
top_bound_mm <- filter_top(mean_l$Mouse_mom)
top_dist_mm <- dist_mm[top_bound_mm$Symbol, ]
top_summ_mm <- dist_summ_mm[top_bound_mm$Symbol, ]


# Relationship between count of DHSs around genes and mean binding


tt <- left_join(mean_l$Human_mom, pc_dhs_count, by = "Symbol")

ggplot(tt, aes(x = Count, y = Mean)) +
  geom_point(alpha = 0.4, colour = "royalblue") +
  xlab("Count of DHS elements around TSS") +
  ylab("Mean binding score") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.margin = margin(10, 20, 10, 10))

cor(tt$Mean, tt$Count, method = "spearman")


# Binding but no DHS

tt %>% 
  filter(Count == 0) %>% 
  arrange(desc(Mean)) %>% 
  head(20)

pc_dhs_count_top <- mutate(pc_dhs_count, Top = Symbol %in% top_bound_hg$Symbol)
boxplot(pc_dhs_count_top$Count ~ pc_dhs_count_top$Top)
wilcox.test(pc_dhs_count_top$Count ~ pc_dhs_count_top$Top)



# Rarely/never bound genes
# ------------------------------------------------------------------------------


filter_bottom <- function(mean_df, qtl = 0.01) {
  filter(mean_df, Mean <= quantile(Mean, qtl)) %>% arrange(Mean)
}


# TODO: finalize mat for infreq calc
mat_raw <- dat$Permissive_hg$Mat_raw
mat_qnl <- dat$Permissive_hg$Mat_qnl
mean_raw <- get_all_mean(get_tf_mean(mat_raw, dat$Permissive_hg$Meta))
mean_qnl <- get_all_mean(get_tf_mean(mat_qnl, dat$Permissive_hg$Meta))
# plot(mean_raw$Mean, mean_qnl$Mean)
# cor(mean_raw$Mean, mean_qnl$Mean)

btm_bound_hg1 <- filter_bottom(mean_raw)
btm_bound_hg2 <- filter_bottom(mean_qnl)
length(intersect(btm_bound_hg1$Symbol, btm_bound_hg2$Symbol))


dist_summ_hg[btm_bound_hg1$Symbol, ]


# Example of infreq bound gene that still has a proximal peak
dist_summ_hg["CT45A6", ]
summary(mat_raw["CT45A6", ])
head(sort(mat_raw["CT45A6", ], decreasing = TRUE))
summary(mat_qnl["CT45A6", ])
head(sort(mat_qnl["CT45A6", ], decreasing = TRUE))


# Example of infreq bound for DHS

pc_dhs_count2 <- mutate(pc_dhs_count,
                        Raw = Symbol %in% btm_bound_hg1$Symbol,
                        QNL = Symbol %in% btm_bound_hg2$Symbol)

pc_dhs_count2 %>% head
boxplot(pc_dhs_count2$Count ~ pc_dhs_count2$Raw)
boxplot(pc_dhs_count2$Count ~ pc_dhs_count2$QNL)


# For frequent bound, are peaks dispersed?

pc_sub <- pc_hg[pc_hg$Symbol == "MIDN"]

downstream <- dhs_gr[dhs_gr@seqnames == as.numeric(pc_sub@seqnames) & dhs_gr@ranges < pc_sub@ranges]
upstream <- dhs_gr[dhs_gr@seqnames == as.numeric(pc_sub@seqnames) & dhs_gr@ranges > pc_sub@ranges]
downstream$Distance <- distance(pc_sub, downstream)
upstream$Distance <- distance(pc_sub, upstream)





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



##

# limma voom for in vs out


# prepare data matrices - min count filter on raw/no-norm score matrix (will be 
# QN+log in voom). for now currently just from examining distn of row/gene sums

sum_hg <- rowSums(bmat_dedup_hg)
hist(sum_hg, breaks = 100)
min_hg <- 50
abline(v = min_hg, col = "red")

# Subset matrix to only minimum raw counts across experiments
keep_hg <- sum_hg > min_hg
mat_hg <- bmat_dedup_hg[keep_hg, meta_final_hg$File]


# Design matrices, voom, limma model
# means model (~0) so each TF has a coef corresponding to group mean with no intercept
#-------------------------------------------------------------------------------


design_hg <- model.matrix(
  ~ 0 + Symbol + log10(N), 
  data = meta_final_hg)

rownames(design_hg) <- meta_final_hg$File
colnames(design_hg) <- str_replace(colnames(design_hg), "Symbol", "")

voom_hg <- voom(mat_hg, 
                design = design_hg, 
                normalize.method = "quantile")

fit_hg <- lmFit(voom_hg, design = design_hg)



# interested in each TR vs the rest, so must iteratively construct the 
# appropriate contrast vector. (Law et al., 2020) used as reference
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/ 
#-------------------------------------------------------------------------------


contr_list <- function(meta, contr_vec) {
  # Make a contrast vector for each TF vs the rest: 1 for TF of interest, 
  # -(1/(#TF - 1)) for the rest, leaving nuisance variables as 0
  
  tfs <- unique(meta$Symbol)
  
  clist <- lapply(tfs, function(x) {
    contr_vec[x] <- 1
    contr_vec[names(contr_vec) %in% setdiff(tfs, x)] <- -(1/(length(tfs)-1))
    return(contr_vec)
  })
  names(clist) <- tfs
  return(clist)
}


contr_hg <- rep(0, length(colnames(coef(fit_hg))))
names(contr_hg) <- colnames(coef(fit_hg))

clist_hg <- contr_list(meta_final_hg, contr_hg)


# For each symbol, get the top results for the symbol vs all contrast
#-------------------------------------------------------------------------------


top_fit <- function(contr_list, fit) {
  # Apply eBayes + contrast fit to the model and extract toptable for each contr
  lapply(contr_list, function(x) {
    contr_fit <- eBayes(contrasts.fit(fit, x))
    topTable(contr_fit, n = Inf)
  })
}


top_hg <- top_fit(clist_hg, fit_hg) 



saveRDS(
  list(
    Mat_raw = bmat_dedup_hg,
    Mat_QNL = qmat_hg,
    Meta = meta_final_hg,
    All_mean = all_mean,
    Group_mean = group_mean,
    Fit = top_hg
  ),
  file = "~/scratch/R_objects/unibind_bindscore_human.RDS"
)





# count of pos FC sig genes at FDR05
n_hg <- sapply(top_hg, function(x) nrow(filter(x, logFC > 0 & adj.P.Val < 0.05)))

n_df <- data.frame(Count_diff = n_hg) %>% 
  rownames_to_column(var = "Symbol") %>% 
  left_join(count(meta_final_hg, Symbol), by = "Symbol")

ggplot(n_df, aes(x =  reorder(Symbol, Count_diff), y = Count_diff)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust=1))


# head(top_hg$RUNX1)
# view(top_hg$RUNX1)

# tf <- "RUNX1"
# gene <- "GTF2A2"

# boxplot(qmat_hg[gene, ] ~ meta_final_hg$Symbol == tf)
# boxplot(voom_hg$E[gene, ] ~ meta_final_hg$Symbol == tf)
# 
# plot(density(qmat_hg[gene, meta_final_hg$Symbol == tf]), col = "red")
# lines(density(qmat_hg[gene, meta_final_hg$Symbol != tf]), col = "black")

