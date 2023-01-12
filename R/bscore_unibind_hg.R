## TODO: DOC

library(tidyverse)
library(parallel)
library(WGCNA)
library(limma)
source("~/regnetR/R/utils/bscore_unibind_functions.R")

path_hg <- "/home/amorin/Data/Peak_files/Unibind/Hg38/Permissive/damo_hg38_TFBS_per_TF/"
pc_hg <- read.delim("~/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)
bl_hg <- read.delim("~/Data/Chromosome_info/blacklist_hg38.tsv", stringsAsFactors = FALSE)
meta_hg <- read.delim("~/Data/Metadata/Chipseq/unibind_human.tsv", stringsAsFactors = FALSE)
meta_rob_hg <- read.delim("~/Data/Metadata/Chipseq/unibind_human.tsv", stringsAsFactors = FALSE)

# Remove duplicated pseudoautosomal genes (keep X copy)
dupl <- pc_hg$Symbol[duplicated(pc_hg$Symbol)]
pc_hg <- filter(pc_hg, !(Symbol %in% dupl & Chromosome == "Y"))

# Range tables -> GR objects
bl_gr_hg <- bl_to_gr(bl_hg)
pc_gr_hg <- pc_to_gr(pc_hg)


# Remove experiments under peak cutoff
meta_hg <- filter(meta_hg, N >= 100)


# Scoring all is slow - took ~ 30 hours on frink with 8 cores

if(!file.exists("~/scratch/R_objects/unibind_all_scores_hg.RDS")) {
  bscore_hg <- load_and_score(path_hg, meta_hg, bl_gr_hg, pc_gr_hg)
  names(bscore_hg) <- meta_hg$File
  saveRDS(bscore_hg, file = "~/scratch/R_objects/unibind_all_scores_hg.RDS")
} else {
  bscore_hg <- readRDS("~/scratch/R_objects/unibind_all_scores_hg.RDS")
}


bmat_hg <- do.call(cbind, bscore_hg)


# Looking for duplicate samples with multiple motifs


meta_dedup_hg <- meta_hg[meta_hg$ID %in% meta_hg$ID[duplicated(meta_hg$ID)], ]


# Get correlation of duplicates


dup_cor <- lapply(unique(meta_dedup_hg$ID), function(x) {
  cols <- filter(meta_dedup_hg, ID == x)
  mat <- bmat_hg[, cols$File]
  res <- cor(mat)
  res <- unique(res[res != 1])
  data.frame(ID = x,
             Cor = res,
             Symbol = cols$Symbol[1])
})

dup_cor <- do.call(rbind, dup_cor)


p1 <- 
  ggplot(dup_cor, aes(x  = Symbol, y = Cor)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Pcor of binding scores") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# summary(dup_cor$Cor)


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


# Quantile norm then describe in vs out cor



qmat_hg <- preprocessCore::normalize.quantiles(
  log10(bmat_dedup_hg+1), keep.names = TRUE)

cor_mat <- WGCNA::cor(qmat_hg)

cor_df <- mat_to_df(cor_mat) %>% 
  mutate(TF1 = str_split(Row, "\\.", simplify = TRUE)[, 3],
         TF2 = str_split(Col, "\\.", simplify = TRUE)[, 3],
         Group = (TF1 == TF2))


boxplot(cor_df$Value ~ cor_df$Group, 
        xlab = "Same TF", ylab = "Pcor",
        cex.lab = 1.5)


cor_sum <- mclapply(unique(meta_final_hg$Symbol), function(x) {
  df <- filter(cor_df, TF1 == x | TF2 == x)
  m_out <- mean(filter(df, !Group)$Value)
  m_in <- ifelse(nrow(filter(df, Group)) == 0, NA, mean(filter(df, Group)$Value))
  data.frame(Symbol = x, In = m_in, Out = m_out)
}, mc.cores = 8)


cor_sum <- bind_rows(cor_sum) %>% 
  left_join(count(meta_final_hg, Symbol), by = "Symbol")


plot(y = cor_sum$In, x = cor_sum$n,
     ylab = "Intra-TR Pcor", xlab = "Count of samples",
     cex.lab = 1.5)


# describe mean/var across all, and get group means


all_mean <- data.frame(Symbol = rownames(qmat_hg),
                       Mean = rowMeans(qmat_hg))


group_mean <- lapply(unique(meta_final_hg$Symbol), function(x) {
  
  meta <- filter(meta_final_hg, Symbol == x)
  data.frame(Symbol = rownames(qmat_hg),
             Mean = rowMeans(qmat_hg[, meta$File, drop = FALSE]))
})
names(group_mean) <- unique(meta_final_hg$Symbol)


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
