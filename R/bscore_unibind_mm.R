## TODO: DOC

library(tidyverse)
library(parallel)
library(WGCNA)
library(limma)
source("~/regnetR/R/utils/bscore_unibind_functions.R")

path_mm <- "/home/amorin/Data/Peak_files/Unibind/Mm10/Permissive/damo_mm10_TFBS_per_TF/"
pc_mm <- read.delim("~/Data/Metadata/refseq_select_mm10.tsv", stringsAsFactors = FALSE)
bl_mm <- read.delim("~/Data/Chromosome_info/blacklist_mm10.tsv", stringsAsFactors = FALSE)
meta_mm <- read.delim("~/Data/Metadata/Chipseq/unibind_mouse.tsv", stringsAsFactors = FALSE)


# Range tables -> GR objects
bl_gr_mm <- bl_to_gr(bl_mm)
pc_gr_mm <- pc_to_gr(pc_mm)


# Remove experiments under peak cutoff
meta_mm <- filter(meta_mm, N >= 100)


# Scoring all is slow - took ~ 30 hours on frink with 8 cores

if(!file.exists("~/scratch/R_objects/unibind_all_scores_mm.RDS")) {
  bscore_mm <- load_and_score(path_mm, meta_mm, bl_gr_mm, pc_gr_mm)
  names(bscore_mm) <- meta_mm$File
  saveRDS(bscore_mm, file = "~/scratch/R_objects/unibind_all_scores_mm.RDS")
} else {
  bscore_mm <- readRDS("~/scratch/R_objects/unibind_all_scores_mm.RDS")
}


bmat_mm <- do.call(cbind, bscore_mm)


# Looking for duplicate samples with multiple motifs


meta_dedup_mm <- meta_mm[meta_mm$ID %in% meta_mm$ID[duplicated(meta_mm$ID)], ]


# Get correlation of duplicates


dup_cor <- lapply(unique(meta_dedup_mm$ID), function(x) {
  cols <- filter(meta_dedup_mm, ID == x)
  mat <- bmat_mm[, cols$File]
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


dup_out <- bmat_mm[, !colnames(bmat_mm) %in% meta_dedup_mm$File]

dup_avg <- lapply(unique(meta_dedup_mm$ID), function(x) {
  cols <- filter(meta_dedup_mm, ID == x)
  mat <- bmat_mm[, cols$File]
  rowMeans(mat)
})

dup_avg <- do.call(cbind, dup_avg)

meta_dedup_mm <- meta_dedup_mm %>% 
  mutate(File = str_replace(File, "\\.MA.*", "\\.avg")) %>% 
  distinct(File, .keep_all = TRUE)

stopifnot(identical(nrow(meta_dedup_mm), ncol(dup_avg)))

colnames(dup_avg) <- meta_dedup_mm$File


# Combine average and dup out to get de-duplicated matrix and meta


bmat_dedup_mm <- cbind(dup_out, dup_avg)

meta_final_mm <- filter(meta_mm, !(ID %in% meta_dedup_mm$ID)) %>% 
  rbind(meta_dedup_mm) %>% 
  arrange(Symbol)


bmat_dedup_mm <- bmat_dedup_mm[, meta_final_mm$File]


stopifnot(identical(
  bmat_dedup_mm[, "EXP030565.embryonic_fibroblasts.ASCL1.avg"],
  rowMeans(bmat_mm[, c("EXP030565.embryonic_fibroblasts.ASCL1.MA1100.2.damo.bed",
                       "EXP030565.embryonic_fibroblasts.ASCL1.MA1631.1.damo.bed")])
))


# Quantile norm then describe in vs out cor


####### TODO: remove when re-ran
# meta_final_mm <- filter(meta_final_mm, N >= 100)
# bmat_dedup_mm <- bmat_dedup_mm[, meta_final_mm$File]
########


qmat_mm <- preprocessCore::normalize.quantiles(
  log10(bmat_dedup_mm+1), keep.names = TRUE)

cor_mat <- WGCNA::cor(qmat_mm)

cor_df <- mat_to_df(cor_mat) %>% 
  mutate(TF1 = str_split(Row, "\\.", simplify = TRUE)[, 3],
         TF2 = str_split(Col, "\\.", simplify = TRUE)[, 3],
         Group = (TF1 == TF2))


boxplot(cor_df$Value ~ cor_df$Group)


cor_sum <- mclapply(unique(meta_final_mm$Symbol), function(x) {
  df <- filter(cor_df, TF1 == x | TF2 == x)
  m_out <- mean(filter(df, !Group)$Value)
  m_in <- ifelse(nrow(filter(df, Group)) == 0, NA, mean(filter(df, Group)$Value))
  data.frame(Symbol = x, In = m_in, Out = m_out)
}, mc.cores = 8)


cor_sum <- bind_rows(cor_sum) %>% 
  left_join(count(meta_final_mm, Symbol), by = "Symbol")


# describe mean/var across all, and get group means


all_mean <- data.frame(Symbol = rownames(qmat_mm),
                       Mean = rowMeans(qmat_mm))


group_mean <- lapply(unique(meta_final_mm$Symbol), function(x) {
  
  meta <- filter(meta_final_mm, Symbol == x)
  data.frame(Symbol = rownames(qmat_mm),
             Mean = rowMeans(qmat_mm[, meta$File, drop = FALSE]))
})
names(group_mean) <- unique(meta_final_mm$Symbol)


# limma voom for in vs out


# prepare data matrices - min count filter on raw/no-norm score matrix (will be 
# QN+log in voom). for now currently just from examining distn of row/gene sums

sum_mm <- rowSums(bmat_dedup_mm)
hist(sum_mm, breaks = 100)
min_mm <- 60
abline(v = min_mm, col = "red")

# Subset matrix to only minimum raw counts across experiments
keep_mm <- sum_mm > min_mm
mat_mm <- bmat_dedup_mm[keep_mm, meta_final_mm$File]


# Design matrices, voom, limma model
# means model (~0) so each TF has a coef corresponding to group mean with no intercept
#-------------------------------------------------------------------------------


design_mm <- model.matrix(
  ~ 0 + Symbol + log10(N), 
  data = meta_final_mm)

rownames(design_mm) <- meta_final_mm$File
colnames(design_mm) <- str_replace(colnames(design_mm), "Symbol", "")

voom_mm <- voom(mat_mm, 
                design = design_mm, 
                normalize.method = "quantile")

fit_mm <- lmFit(voom_mm, design = design_mm)



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


contr_mm <- rep(0, length(colnames(coef(fit_mm))))
names(contr_mm) <- colnames(coef(fit_mm))

clist_mm <- contr_list(meta_final_mm, contr_mm)


# For each symbol, get the top results for the symbol vs all contrast
#-------------------------------------------------------------------------------


top_fit <- function(contr_list, fit) {
  # Apply eBayes + contrast fit to the model and extract toptable for each contr
  lapply(contr_list, function(x) {
    contr_fit <- eBayes(contrasts.fit(fit, x))
    topTable(contr_fit, n = Inf)
  })
}


top_mm <- top_fit(clist_mm, fit_mm) 


# count of pos FC sig genes at FDR05
n_mm <- sapply(top_mm, function(x) nrow(filter(x, logFC > 0 & adj.P.Val < 0.05)))


saveRDS(list(
  Mat_raw = bmat_dedup_mm,
  Mat_QNL = qmat_mm,
  Meta = meta_final_mm,
  All_mean = all_mean,
  Group_mean = group_mean,
  Fit = top_mm),
  file = "~/scratch/R_objects/unibind_bindscore_mouse.RDS"
)
