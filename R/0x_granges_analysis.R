## This script is a placeholder for work in progress analysis using Granges
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

# Load bound regions as GR objects
gr_hg <- readRDS(gr_path_perm_hg)
gr_mm <- readRDS(gr_path_perm_mm)

# Load ENCODE DHS regions and convert to GR object
dhs_hg <- read.delim("~/Data/Chromosome_info/ENCODE_human_DHS.tsv", stringsAsFactors = FALSE)
dhs_hg <- mutate(dhs_hg, seqname = str_replace(seqname, "chr", ""))
dhs_hg <- makeGRangesFromDataFrame(dhs_hg, keep.extra.columns = TRUE)

# TODO: keep?
dist_mat_out <- "~/scratch/R_objects/unibind_dist_to_nearest_mats.RDS" 
dhs_ol_out <- "~/scratch/R_objects/unibind_dhs_proportion_overlap.RDS" 


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


nearest_dist_summary <- function(mat) {
  df <- data.frame(t(apply(mat, 1, summary)))
  return(df)
}


if (!file.exists(dist_mat_out)) {
  dist_hg <- nearest_dist_mat(pc = pc_hg, gr_l = gr_hg)
  dist_mm <- nearest_dist_mat(pc = pc_mm, gr_l = gr_mm)
  saveRDS(list(Human = dist_hg, Mouse = dist_mm), dist_mat_out)
} else {
  temp <- readRDS(dist_mat_out)
  dist_hg <- temp$Human
  dist_mm <- temp$Mouse
  rm(temp)
}


dist_summ_hg <- nearest_dist_summary(dist_hg)
dist_summ_mm <- nearest_dist_summary(dist_mm)


# Use ENCODE DHS index (human only) to count open chromatin sites around genes
# ------------------------------------------------------------------------------


# the window in bp +/- from the TSS to consider overlaps with DHSs 
# TODO: to function

dhs_gap <- 50e3

pc_dhs_count <- data.frame(
  Symbol = pc_hg$Symbol,
  Count = countOverlaps(pc_hg, dhs_hg, maxgap = dhs_gap)
)

pc_dhs_count %>% head
summary(pc_dhs_count$Count)
hist(pc_dhs_count$Count, breaks = 100)


# count DHS around Pax6 over range of distances
# TODO: review function calls


dhs_distance <- function(gene, pc_gr, dhs_gr) {
  
  pc_sub <- pc_gr[pc_gr$Symbol == gene]
  
  downstream <- dhs_gr[dhs_gr@seqnames == as.numeric(pc_sub@seqnames) & dhs_gr@ranges < pc_sub@ranges]
  upstream <- dhs_gr[dhs_gr@seqnames == as.numeric(pc_sub@seqnames) & dhs_gr@ranges > pc_sub@ranges]
  downstream$Distance <- distance(pc_sub, downstream)
  upstream$Distance <- distance(pc_sub, upstream)
  
  return(list(Distance_upstream = upstream,
              Distance_downstream = downstream))
  
}

# PAX6 top target "DNAJB6"
dhs_dist <- dhs_distance("PAX6", pc_hg, dhs_hg)


# TODO: function
dist_l <- lapply(seq(from = 1e3, to = 1e5, by = 1e3), function(dist) {
  data.frame(
    Distance = dist,
    Count_upstream = length(dhs_dist[[1]][dhs_dist[[1]]$Distance <= dist]),
    Count_downstream = length(dhs_dist[[2]][dhs_dist[[2]]$Distance <= dist])
  )
})



dist_df <- do.call(rbind, dist_l)


ggplot(dist_df) +
  geom_line(aes(x = Distance/1e3, y = Count_upstream, col = "Count_upstream"), 
            linewidth = 1.4) +
  geom_line(aes(x = Distance/1e3, y = Count_downstream, col = "Count_downstream"), 
            linewidth = 1.4) +
  ylab("Count of DHS elements") +
  xlab("Distance from TSS (kb)") +
  scale_colour_manual(values = c("Count_upstream" = "red", 
                                 "Count_downstream" = "black")) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.title = element_blank(),
        legend.text = element_text(size = 15))



# Proportion of peaks in each human experiment that overlap a DHS
# ------------------------------------------------------------------------------


dhs_overlap <- function(gr_l, dhs, ncores) {
  
  mclapply(gr_l, function(x) {
    n_peaks <- length(x)
    ol <- findOverlaps(x, dhs_hg)
    n_distinct(ol@from) / n_peaks
  }, mc.cores = ncores)
}


if (!file.exists(dhs_ol_out)) {
  prop_ol <- dhs_overlap(gr_hg, dhs_hg, cores)
  saveRDS(prop_ol, dhs_ol_out)
} else {
  prop_ol <- readRDS(dhs_ol_out)
}


# TODO: consolidated meta
stopifnot(identical(names(prop_ol), meta_hg$File))


prop_ol_df <- left_join(
  meta_hg,
  rownames_to_column(data.frame(Proportion = unlist(prop_ol)), var = "File"),
  by = "File"
)


prop_summ <- prop_ol_df %>% 
  group_by(Symbol) %>% 
  dplyr::summarise(Mean = mean(Proportion, na.rm = TRUE), 
                   n = length(File)) %>% 
  arrange(desc(Mean))

head(prop_summ)
tail(prop_summ)

summary(prop_ol_df$Proportion)


plot(prop_summ$Mean, prop_summ$n)

ggplot(prop_ol_df, aes(x = reorder(Symbol, Proportion, FUN = median), y = Proportion)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Proportion of peaks overlapping a DHS") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))



# Relationship between count of DHSs around genes and mean binding


dhs_mean_df <- left_join(mean_l$Human_mom, pc_dhs_count, by = "Symbol")

ggplot(dhs_mean_df, aes(x = Count, y = Mean)) +
  geom_point(alpha = 0.4, colour = "royalblue") +
  xlab("Count of DHS elements around TSS") +
  ylab("Mean binding score") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.margin = margin(10, 20, 10, 10))

cor(dhs_mean_df$Mean, dhs_mean_df$Count, method = "spearman")


# Binding but no DHS

dhs_mean_df %>% 
  filter(Count == 0) %>% 
  arrange(desc(Mean)) %>% 
  head(20)

pc_dhs_count_top <- mutate(pc_dhs_count, Top = Symbol %in% top_bound_hg$Symbol)
boxplot(pc_dhs_count_top$Count ~ pc_dhs_count_top$Top)
wilcox.test(pc_dhs_count_top$Count ~ pc_dhs_count_top$Top)
